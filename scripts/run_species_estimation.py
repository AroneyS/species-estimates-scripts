#!/usr/bin/env python3

###############################################################################
#
#    Collate OTU sequences per marker across archive OTU tables, then cluster
#    with smafa per marker.
#
#    Steps:
#      1. Read archive OTU table(s) (JSON / JSON.gz)
#      2. Filter off-target OTUs (domain mismatch for the marker)
#      3. Order by GlobDB-known taxonomy first, then num_hits descending
#      4. Per marker: smafa cluster --max-divergence 2
#
#    With --run-through-mqsub:
#      - Steps 1-3 run locally, writing per-marker FASTA files
#      - Steps 4-5 are submitted as one mqsub job per marker
#
###############################################################################

import argparse
import concurrent.futures
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from multiprocessing import Pool

import extern
import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from singlem.metapackage import Metapackage

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Subdirectory within output_dir where per-sample per-marker FASTAs land
# after extraction (one subdir per sample, sharded by first 5 chars).
SAMPLE_FASTA_SUBDIR = "sample_fastas"

# Subdirectory within output_dir where per-sample FASTAs are cat'd into one
# file per marker before clustering.
FASTA_SUBDIR = "collated_fastas"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DOMAIN_PREFIXES = {
    "Bacteria":  "d__Bacteria",
    "Archaea":   "d__Archaea",
    "Eukaryota": "d__Eukaryota",
}


def open_archive(path):
    """Return a file-like for a plain or gzip-compressed archive OTU table."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_archive_otus(path):
    """
    Parse a SingleM archive OTU table (v4 JSON format).

    Returns a list of dicts, one per OTU, with keys matching the 'fields'
    array plus 'source_file'.
    """
    with open_archive(path) as fh:
        data = json.load(fh)

    fields = data["fields"]
    records = []
    for otu in data["otus"]:
        # otus array entries are parallel to fields
        rec = dict(zip(fields, otu))
        rec["source_file"] = path
        records.append(rec)
    return records


def build_marker_domain_map(metapackage_path):
    """
    Return a dict: marker_name -> set of domains (e.g. {"Bacteria", "Archaea"})
    using the singlem metapackage.
    """
    mp = Metapackage.acquire(metapackage_path)
    marker_domains = {}
    for spkg in mp.singlem_packages:
        name = spkg.graftm_package_basename()
        marker_domains[name] = set(spkg.target_domains())
    return marker_domains


def taxonomy_domain(taxonomy_str):
    """
    Extract the domain from a taxonomy string such as
    'Root; d__Bacteria; p__Firmicutes; ...'
    Returns one of 'Bacteria', 'Archaea', 'Eukaryota', or None.
    """
    if not taxonomy_str:
        return None
    for domain, prefix in DOMAIN_PREFIXES.items():
        if prefix in taxonomy_str:
            return domain
    return None


def is_on_target(otu, marker_domains):
    """
    Return True if the OTU's taxonomy domain is compatible with the marker's
    target domains.  OTUs with no taxonomy (unclassified) are kept.
    """
    allowed = marker_domains.get(otu["gene"], set())
    if not allowed:
        # Unknown marker — keep to be safe
        return True
    dom = taxonomy_domain(otu.get("taxonomy", ""))
    if dom is None:
        # Unclassified; retain
        return True
    return dom in allowed


def sort_key(otu):
    """
    Sort so that query-assigned OTUs come first, then by num_hits descending.

    'taxonomy_assignment_method' == 'singlem_query_based'  ->  known (lower sort value = first).
    """
    known_rank = 0 if otu.get("taxonomy_assignment_method") == "singlem_query_based" else 1
    hits = otu.get("num_hits", 0)
    try:
        hits = int(hits)
    except (TypeError, ValueError):
        hits = 0
    return (known_rank, -hits)


def otus_to_fasta(otus):
    """
    Convert a list of OTU dicts to a FASTA string.

    Sort fields are embedded in the header so that after per-sample FASTAs are
    cat'd together the global sort can be reconstructed without re-reading the
    original archive.  Format:
        >otu{i}|{sample}|{gene}|unknown={0or1}|hits={n}

    unknown=0: taxonomy_assignment_method is singlem_query_based; unknown=1: everything else.
    """
    lines = []
    for i, otu in enumerate(otus):
        seq = otu.get("sequence", "")
        if not seq or seq == "-":
            continue
        known_rank = 0 if otu.get("taxonomy_assignment_method") == "singlem_query_based" else 1
        hits = 0
        try:
            hits = int(otu.get("num_hits", 0))
        except (TypeError, ValueError):
            pass
        header = f">otu{i}|{otu['sample']}|{otu['gene']}|unknown={known_rank}|hits={hits}"
        lines.append(header)
        lines.append(seq)
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# smafa wrappers and clustering helpers
# ---------------------------------------------------------------------------

def _parse_header_stats(header):
    """Extract (unknown_rank, hits, n_occurrences, sum_hits) from a FASTA header.

    n_occurrences and sum_hits are written by _dedup_fasta_with_counts; for
    original per-sample headers they default to 1 and hits respectively.
    """
    try:
        parts = header.lstrip(">").split("|")
        unknown_rank = int(next((p.split("=")[1] for p in parts if p.startswith("unknown=")), "1"))
        hits = int(next((p.split("=")[1] for p in parts if p.startswith("hits=")), "0"))
        n_occurrences = int(next((p.split("=")[1] for p in parts if p.startswith("n_occurrences=")), "1"))
        sum_hits = int(next((p.split("=")[1] for p in parts if p.startswith("sum_hits=")), str(hits)))
    except (ValueError, IndexError):
        unknown_rank, hits, n_occurrences, sum_hits = 1, 0, 1, 0
    return unknown_rank, hits, n_occurrences, sum_hits


def _dedup_fasta_with_counts(collated_path, dedup_path):
    """
    Deduplicate a large FASTA by sequence, aggregating occurrence counts and
    summed hits across all occurrences.  Uses disk-based sort (GNU sort) to
    keep orchestrator memory at O(1) regardless of input size.

    Output sequences are in the same order as their FIRST occurrence in the
    collated FASTA (original sort-by-taxonomy/hits order is preserved), which
    is important because smafa's greedy algorithm is order-dependent.

    Writes a deduplicated FASTA where each header encodes:
      |unknown={min_unknown}|hits={max_hits}|n_occurrences={count}|sum_hits={sum}

    _parse_header_stats reads these fields back so that cluster_size and
    total_hits in the final representatives.fasta reflect the full dataset.
    """
    import subprocess
    tmp_tsv         = dedup_path + ".counts_tmp.tsv"
    tmp_seq_sorted  = dedup_path + ".counts_seq_sorted.tsv"
    tmp_agg         = dedup_path + ".counts_agg.tsv"
    tmp_pos_sorted  = dedup_path + ".counts_pos_sorted.tsv"
    try:
        # Pass 1: stream FASTA -> write one TSV row per sequence:
        #   position TAB sequence TAB hits TAB unknown_rank  (O(1) memory).
        logging.info("Building sequence-counts TSV ...")
        pos = 0
        with open(collated_path) as in_fh, open(tmp_tsv, "w") as out_fh:
            cur_hdr   = None
            seq_parts = []
            for line in in_fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if cur_hdr is not None and seq_parts:
                        seq = "".join(seq_parts)
                        unk, hits, _, _ = _parse_header_stats(cur_hdr)
                        out_fh.write(f"{pos}\t{seq}\t{hits}\t{unk}\n")
                        pos += 1
                    cur_hdr   = line
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if cur_hdr is not None and seq_parts:
                seq = "".join(seq_parts)
                unk, hits, _, _ = _parse_header_stats(cur_hdr)
                out_fh.write(f"{pos}\t{seq}\t{hits}\t{unk}\n")

        # Sort by sequence to group identical sequences for aggregation.
        logging.info("Sorting by sequence for aggregation ...")
        subprocess.run(
            ["sort", "-k2,2", "--buffer-size=1G", tmp_tsv, "-o", tmp_seq_sorted],
            check=True,
        )

        # Pass 2: aggregate counts per unique sequence; retain the minimum position
        # (= first occurrence) so we can restore original order afterwards.
        logging.info("Aggregating counts ...")
        cur_seq      = None
        first_pos    = None
        count        = 0
        sum_hits     = 0
        min_unknown  = 1
        max_hits_cur = 0

        with open(tmp_seq_sorted) as in_fh, open(tmp_agg, "w") as out_fh:
            def _flush_agg(out_fh):
                if cur_seq is None:
                    return
                out_fh.write(
                    f"{first_pos}\t{cur_seq}\t{min_unknown}\t{max_hits_cur}"
                    f"\t{count}\t{sum_hits}\n"
                )

            for line in in_fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                pos_str, seq, hits_str, unk_str = parts[0], parts[1], parts[2], parts[3]
                pos  = int(pos_str)
                hits = int(hits_str)
                unk  = int(unk_str)
                if seq == cur_seq:
                    count        += 1
                    sum_hits     += hits
                    min_unknown   = min(min_unknown, unk)
                    max_hits_cur  = max(max_hits_cur, hits)
                    first_pos     = min(first_pos, pos)
                else:
                    _flush_agg(out_fh)
                    cur_seq       = seq
                    first_pos     = pos
                    count         = 1
                    sum_hits      = hits
                    min_unknown   = unk
                    max_hits_cur  = hits
            _flush_agg(out_fh)

        # Sort aggregated rows by first_pos to restore original input order.
        logging.info("Restoring original sequence order ...")
        subprocess.run(
            ["sort", "-k1,1n", "--buffer-size=1G", tmp_agg, "-o", tmp_pos_sorted],
            check=True,
        )

        # Pass 3: write deduplicated FASTA in original order.
        logging.info("Writing deduplicated FASTA ...")
        seq_num = 0
        with open(tmp_pos_sorted) as in_fh, open(dedup_path, "w") as out_fh:
            for line in in_fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                _, seq, min_unk, max_hits, n_occ, s_hits = parts
                out_fh.write(
                    f">seq{seq_num}"
                    f"|unknown={min_unk}"
                    f"|hits={max_hits}"
                    f"|n_occurrences={n_occ}"
                    f"|sum_hits={s_hits}\n{seq}\n"
                )
                seq_num += 1

        logging.info(f"Deduplicated FASTA written with {seq_num} unique sequence(s).")
    finally:
        for f in [tmp_tsv, tmp_seq_sorted, tmp_agg, tmp_pos_sorted]:
            if os.path.exists(f):
                os.unlink(f)


def _run_smafa_cluster_file(fasta_path, max_divergence, marker_name, threads=1):
    """
    Run smafa cluster directly on fasta_path, writing output to a sibling temp file.
    Returns the path to the cluster TSV file (caller is responsible for deletion).
    """
    cluster_tmp = fasta_path + ".smafa_clusters.tsv"
    cmd = [
        "cargo", "run", "--manifest-path", "/home/aroneys/src/smafa/Cargo.toml", "--", "cluster",
        # "smafa", "cluster",
        "--max-divergence", str(max_divergence),
        "--threads", str(threads),
        "--input", fasta_path,
    ]
    logging.debug(f"[{marker_name}] Running: {' '.join(cmd)}")
    with open(cluster_tmp, "w") as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        if os.path.exists(cluster_tmp):
            os.unlink(cluster_tmp)
        raise RuntimeError(
            f"smafa cluster failed for marker {marker_name}:\n{result.stderr}"
        )
    return cluster_tmp


def _read_fasta_sequences(fasta_path):
    """
    Read a FASTA file and return (sequences, seq_to_header, seq_to_stats).
      sequences:     [(header, seq), ...]  in input order
      seq_to_header: seq -> first header seen
      seq_to_stats:  seq -> (unknown_rank, hits, n_occurrences, sum_hits)
    """
    sequences = []
    seq_to_header = {}
    seq_to_stats = {}
    with open(fasta_path) as fh:
        current_header = None
        seq_parts = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None and seq_parts:
                    seq = "".join(seq_parts)
                    sequences.append((current_header, seq))
                    seq_to_header.setdefault(seq, current_header)
                    seq_to_stats.setdefault(seq, _parse_header_stats(current_header))
                current_header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_header is not None and seq_parts:
            seq = "".join(seq_parts)
            sequences.append((current_header, seq))
            seq_to_header.setdefault(seq, current_header)
            seq_to_stats.setdefault(seq, _parse_header_stats(current_header))
    return sequences, seq_to_header, seq_to_stats


def _propagate_chunk_results(rep_to_all_members, chunk_clusters_path):
    """
    Read a smafa clusters TSV (member TAB centroid) and propagate accumulated
    cluster membership from previous rounds.

    For each new centroid, collects all original sequences that were members of
    any previous-round representative now assigned to that centroid.

    Returns the updated rep_to_all_members dict.
    """
    chunk_rep_to_members = defaultdict(list)
    with open(chunk_clusters_path) as cf:
        for line in cf:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            member_seq = parts[0]   # col 1 = sequence being assigned
            true_rep   = parts[1]   # col 2 = centroid
            chunk_rep_to_members[true_rep].append(member_seq)

    new_rep_to_all_members = {}
    for new_rep, chunk_members in chunk_rep_to_members.items():
        all_members = []
        for member in chunk_members:
            if member in rep_to_all_members:
                # member was itself a rep in a previous round; expand to all its members
                all_members.extend(rep_to_all_members[member])
            else:
                all_members.append(member)
        new_rep_to_all_members[new_rep] = all_members
    return new_rep_to_all_members


def _build_cluster_outputs(rep_to_all_members, seq_to_header, seq_to_stats, marker_name):
    """
    Build cluster_tsv and enriched rep_fasta from accumulated cluster state.
    Returns (cluster_tsv_str, rep_fasta_str).
    """
    cluster_tsv_lines = []
    for rep_seq, members in rep_to_all_members.items():
        for member in members:
            cluster_tsv_lines.append(f"{member}\t{rep_seq}")
    cluster_tsv = "\n".join(cluster_tsv_lines) + "\n" if cluster_tsv_lines else ""

    rep_fasta_lines = []
    for rep_seq, members in rep_to_all_members.items():
        n_seqs = sum(seq_to_stats.get(m, (1, 0, 1, 0))[2] for m in members)
        n_query = sum(seq_to_stats.get(m, (1, 0, 1, 0))[2] for m in members
                      if seq_to_stats.get(m, (1, 0, 1, 0))[0] == 0)
        max_hits = max((seq_to_stats.get(m, (1, 0, 1, 0))[1] for m in members), default=0)
        total_hits = sum(seq_to_stats.get(m, (1, 0, 1, 0))[3] for m in members)
        orig_header = seq_to_header.get(rep_seq)
        if orig_header is None:
            logging.warning(f"[{marker_name}] Rep sequence not found in FASTA -- skipping.")
            continue
        enriched_header = (
            f"{orig_header}"
            f"|cluster_size={n_seqs}"
            f"|n_query={n_query}"
            f"|max_hits={max_hits}"
            f"|total_hits={total_hits}"
        )
        rep_fasta_lines.append(enriched_header)
        rep_fasta_lines.append(rep_seq)
    rep_fasta = "\n".join(rep_fasta_lines) + "\n" if rep_fasta_lines else ""
    return cluster_tsv, rep_fasta


def _cluster_marker_file_chunked(fasta_path, max_divergence, marker_name, chunk_size=2000000, threads=16):
    """
    Cluster sequences in sequential greedy chunks (local path, no mqsub).

    Processes sequences in batches of `chunk_size`.  Each round runs smafa on
    (current_reps + next_chunk).  Member lists are propagated so the final
    cluster_tsv maps every original sequence to its final representative.

    Returns (cluster_tsv_str, rep_fasta_str).
    """
    sequences, seq_to_header, seq_to_stats = _read_fasta_sequences(fasta_path)
    n_total = len(sequences)
    if n_total == 0:
        return "", ""

    chunks = [sequences[i:i + chunk_size] for i in range(0, n_total, chunk_size)]
    logging.info(
        f"[{marker_name}] Chunked clustering: {n_total} sequence(s) "
        f"in {len(chunks)} chunk(s) of up to {chunk_size}."
    )

    rep_to_all_members = {}
    current_rep_seqs = []

    for chunk_idx, chunk in enumerate(chunks):
        chunk_seqs = [seq for _, seq in chunk]
        combined_seqs = current_rep_seqs + chunk_seqs

        tmp_fasta_path = None
        cluster_tmp = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False,
                prefix=f"smafa_chunk{chunk_idx}_",
            ) as tmp_fasta:
                for seq in combined_seqs:
                    tmp_fasta.write(seq_to_header[seq] + "\n" + seq + "\n")
                tmp_fasta_path = tmp_fasta.name

            logging.info(
                f"[{marker_name}] Chunk {chunk_idx + 1}/{len(chunks)}: "
                f"{len(current_rep_seqs)} rep(s) + {len(chunk_seqs)} new = "
                f"{len(combined_seqs)} sequence(s)."
            )
            cluster_tmp = _run_smafa_cluster_file(
                tmp_fasta_path, max_divergence, marker_name, threads=threads
            )
            rep_to_all_members = _propagate_chunk_results(rep_to_all_members, cluster_tmp)
        finally:
            if tmp_fasta_path and os.path.exists(tmp_fasta_path):
                os.unlink(tmp_fasta_path)
            if cluster_tmp and os.path.exists(cluster_tmp):
                os.unlink(cluster_tmp)

        current_rep_seqs = list(rep_to_all_members.keys())
        logging.info(
            f"[{marker_name}] After chunk {chunk_idx + 1}: {len(current_rep_seqs)} cluster(s)."
        )

    return _build_cluster_outputs(rep_to_all_members, seq_to_header, seq_to_stats, marker_name)


def _process_cluster_file(fasta_path, cluster_path, marker_name):
    """
    Build the enriched representative FASTA and return (cluster_tsv_str, rep_fasta_str)
    using three streaming passes over the cluster TSV and FASTA files.

    Pass 1 (cluster file):  collect the set of representative sequences
                            (O(n_clusters) memory).
    Pass 2 (FASTA):         build rep_to_header for rep sequences only
                            (O(n_clusters)); build seq_to_stats for all unique
                            sequences as four ints per entry, no header strings
                            (O(n_unique_seqs)).
    Pass 3 (cluster file):  aggregate per-cluster stats (n_seqs, n_query,
                            max_hits, total_hits) using seq_to_stats.
    """
    # Pass 1: get the set of representative sequences.
    # smafa writes  sequence TAB centroid  so col 2 is the centroid/rep.
    rep_seqs = set()
    with open(cluster_path) as cf:
        for line in cf:
            line = line.rstrip("\n")
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    rep_seqs.add(parts[1])

    # Pass 2: stream the FASTA once.
    #   - seq_to_stats: seq -> (unknown_rank, hits, n_occurrences, sum_hits)
    #     for every unique sequence (values are four ints, not header strings)
    #   - rep_to_header: seq -> header  only for sequences that are cluster reps
    seq_to_stats = {}   # seq -> (unknown_rank, hits, n_occurrences, sum_hits)
    rep_to_header = {}  # seq -> original header line (reps only)
    current_header = None
    current_seq_parts = []

    def _fasta_flush(hdr, seq_parts):
        seq = "".join(seq_parts)
        if not seq:
            return
        if seq not in seq_to_stats:
            seq_to_stats[seq] = _parse_header_stats(hdr)
        if seq in rep_seqs and seq not in rep_to_header:
            rep_to_header[seq] = hdr

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    _fasta_flush(current_header, current_seq_parts)
                current_header = line
                current_seq_parts = []
            else:
                current_seq_parts.append(line)
    if current_header is not None:
        _fasta_flush(current_header, current_seq_parts)

    # Pass 3: stream the cluster file to aggregate per-cluster stats.
    rep_stats = {}  # rep_seq -> {n_seqs, n_query, max_hits, total_hits}
    with open(cluster_path) as cf:
        for line in cf:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep_seq = parts[1]           # col 2 = centroid/rep
            member_seq = parts[0]         # col 1 = sequence being assigned
            if rep_seq not in rep_stats:
                rep_stats[rep_seq] = {"n_seqs": 0, "n_query": 0, "max_hits": 0, "total_hits": 0}
            s = rep_stats[rep_seq]
            unk, hits, n_occ, sum_hits_all = seq_to_stats.get(member_seq, (1, 0, 1, 0))
            s["n_seqs"] += n_occ
            s["n_query"] += n_occ if unk == 0 else 0
            s["max_hits"] = max(s["max_hits"], hits)
            s["total_hits"] += sum_hits_all

    # Build enriched rep FASTA.
    rep_fasta_lines = []
    for rep_seq, stats in rep_stats.items():
        orig_header = rep_to_header.get(rep_seq)
        if orig_header is None:
            logging.warning(f"[{marker_name}] Rep sequence not found in FASTA -- skipping.")
            continue
        enriched_header = (
            f"{orig_header}"
            f"|cluster_size={stats['n_seqs']}"
            f"|n_query={stats['n_query']}"
            f"|max_hits={stats['max_hits']}"
            f"|total_hits={stats['total_hits']}"
        )
        rep_fasta_lines.append(enriched_header)
        rep_fasta_lines.append(rep_seq)
    rep_fasta = "\n".join(rep_fasta_lines) + "\n" if rep_fasta_lines else ""

    # Read cluster TSV so callers can write it to disk.
    with open(cluster_path) as cf:
        cluster_tsv = cf.read()

    return cluster_tsv, rep_fasta


# ---------------------------------------------------------------------------
# Phase 1 (local): read one archive, filter, write per-marker FASTAs into a
# sample-specific subdirectory.  This is also the mqsub worker entry point.
# ---------------------------------------------------------------------------

def extract_sample_fastas(archive_path, marker_domains, markers_of_interest, sample_fasta_dir):
    """
    Read a single archive OTU table, filter off-target OTUs, and write one
    FASTA file per marker into sample_fasta_dir/<marker>.fasta.

    Sequences are NOT globally sorted here -- sorting across samples requires
    the full dataset and is done after cat (phase 2).

    Returns the set of marker names that received >=1 sequence.
    """
    otus = read_archive_otus(archive_path)

    # Domain filter
    if marker_domains:
        otus = [o for o in otus if is_on_target(o, marker_domains)]

    # Within-sample sort (GlobDB-known first, then num_hits desc) so that if
    # any sample happens to be processed alone the order is still sensible.
    otus.sort(key=sort_key)

    by_marker = defaultdict(list)
    for otu in otus:
        marker = otu["gene"]
        if markers_of_interest and marker not in markers_of_interest:
            continue
        by_marker[marker].append(otu)

    os.makedirs(sample_fasta_dir, exist_ok=True)
    seen_markers = set()
    for marker_name, marker_otus in by_marker.items():
        fasta_str = otus_to_fasta(marker_otus)
        if not fasta_str.strip():
            continue
        fasta_path = os.path.join(sample_fasta_dir, f"{marker_name}.fasta")
        with open(fasta_path, "w") as fh:
            fh.write(fasta_str)
        seen_markers.add(marker_name)

    # Write sentinel so reruns can skip this archive even when it produced no FASTAs.
    open(os.path.join(sample_fasta_dir, ".done"), "w").close()

    return seen_markers


def extract_sample_fastas_local_worker(args_tuple):
    """multiprocessing.Pool worker wrapper for extract_sample_fastas."""
    (archive_path, marker_domains, markers_of_interest, sample_fasta_dir) = args_tuple
    try:
        seen = extract_sample_fastas(
            archive_path, marker_domains, markers_of_interest, sample_fasta_dir
        )
        logging.debug(f"Extracted {len(seen)} marker(s) from {archive_path}")
        return seen
    except Exception as e:
        logging.error(f"Failed to extract {archive_path}: {e}")
        raise


# ---------------------------------------------------------------------------
# Phase 1 (mqsub): submit one job per archive, wait, then cat locally
# ---------------------------------------------------------------------------

def _sample_fasta_dir(output_dir, archive_path):
    """Deterministic per-sample subdirectory derived from the archive path."""
    basename = os.path.basename(archive_path)
    # Strip common extensions to get a clean sample name
    for ext in (".json.gz", ".json", ".gz"):
        if basename.endswith(ext):
            basename = basename[: -len(ext)]
            break
    # Use first 5 chars as a sharding prefix (same convention as original script)
    return os.path.join(output_dir, SAMPLE_FASTA_SUBDIR, basename[:5], basename)


def _already_extracted(output_dir, archive_path):
    """Return True if extraction has been completed for this archive (done sentinel exists)."""
    sample_dir = _sample_fasta_dir(output_dir, archive_path)
    return os.path.exists(os.path.join(sample_dir, ".done"))


def submit_extraction_jobs(
    archive_paths, marker_domains, markers_of_interest,
    output_dir, this_script_path
):
    """
    Submit one mqsub job per archive OTU table.  Each job calls this script
    with --_extract-sample-archive to produce per-marker FASTAs in a
    sample-specific subdirectory.  Blocks until all jobs finish.
    """
    logging.info(f"Submitting {len(archive_paths)} extraction job(s) via mqsub ...")
    abs_output_dir = os.path.abspath(output_dir)

    # Write a marker->domains TSV that worker jobs can read; workers rebuild
    # marker_domains from this file rather than pickling singlem objects.
    domain_map_path = os.path.join(abs_output_dir, "marker_domains.tsv")
    os.makedirs(abs_output_dir, exist_ok=True)
    with open(domain_map_path, "w") as fh:
        for marker, domains in marker_domains.items():
            fh.write(f"{marker}\t{','.join(sorted(domains))}\n")

    markers_arg = ""
    if markers_of_interest:
        markers_arg = " --markers " + " ".join(sorted(markers_of_interest))

    with tempfile.NamedTemporaryFile(
        mode="w", prefix="extract_mqsub_", suffix=".cmds", delete=False
    ) as cmd_file:
        cmd_file_path = cmd_file.name
        for archive_path in archive_paths:
            sample_dir = _sample_fasta_dir(output_dir, archive_path)
            cmd = (
                f"python3 {this_script_path}"
                f" --_extract-sample-archive {os.path.abspath(archive_path)}"
                f" --_marker-domains-tsv {domain_map_path}"
                f" --_sample-fasta-dir {os.path.abspath(sample_dir)}"
                f" --output-directory {abs_output_dir}"
                + markers_arg
            )
            cmd_file.write(cmd + "\n")

    try:
        mqsub_cmd = (
            f"mqsub -m 4 --name otu_extract"
            f" --segregated-log-files"
            f" --hours 4"
            f" --command-file {cmd_file_path}"
            f" --chunk-size 5000 2>&1"
        )
        logging.info(f"Running: {mqsub_cmd}")
        mqsub_stdout = extern.run(mqsub_cmd)
        logging.info(f"Submitted {len(archive_paths)} extraction job(s).")
        _mqwait(mqsub_stdout)
    finally:
        os.unlink(cmd_file_path)

    logging.info("All extraction jobs finished.")


def _cat_one_marker(args_tuple):
    """Stream all per-sample FASTAs for one marker into the collated output file."""
    marker_name, sample_fastas, collated_path = args_tuple
    with open(collated_path, "w") as out_fh:
        for sfasta in sample_fastas:
            with open(sfasta) as in_fh:
                shutil.copyfileobj(in_fh, out_fh)
    logging.debug(
        f"[{marker_name}] Collated {len(sample_fastas)} sample FASTA(s) into {collated_path}."
    )
    return marker_name, collated_path


def cat_sample_fastas(archive_paths, output_dir, markers_of_interest, threads=1):
    """
    After per-sample extraction jobs have run, cat all per-sample per-marker
    FASTAs into one file per marker.

    If a collated FASTA already exists for a marker it is reused as-is.

    Returns a dict: marker_name -> collated_fasta_path.
    """
    fasta_dir = os.path.join(output_dir, FASTA_SUBDIR)
    os.makedirs(fasta_dir, exist_ok=True)

    # If all markers of interest already have collated FASTAs, skip scanning entirely.
    if markers_of_interest:
        existing = {
            m: os.path.join(fasta_dir, f"{m}.fasta")
            for m in markers_of_interest
            if os.path.exists(os.path.join(fasta_dir, f"{m}.fasta"))
               and os.path.getsize(os.path.join(fasta_dir, f"{m}.fasta")) > 0
        }
        if len(existing) == len(markers_of_interest):
            logging.info(
                f"All {len(existing)} collated FASTA(s) already exist -- skipping concatenation."
            )
            return existing

    # Discover which markers need (re-)collation.
    marker_to_sample_fastas = defaultdict(list)
    for archive_path in archive_paths:
        sample_dir = _sample_fasta_dir(output_dir, archive_path)
        if not os.path.isdir(sample_dir):
            logging.warning(f"No output directory for {archive_path} -- skipping in cat.")
            continue
        for fname in os.listdir(sample_dir):
            if fname.endswith(".fasta"):
                marker_name = fname[: -len(".fasta")]
                if markers_of_interest and marker_name not in markers_of_interest:
                    continue
                marker_to_sample_fastas[marker_name].append(
                    os.path.join(sample_dir, fname)
                )

    if not marker_to_sample_fastas:
        logging.warning("No per-sample FASTA files found after extraction jobs.")
        return {}

    # Skip markers whose collated FASTA already exists.
    marker_to_fasta = {}
    pending_args = []
    for marker_name, sample_fastas in marker_to_sample_fastas.items():
        collated_path = os.path.join(fasta_dir, f"{marker_name}.fasta")
        if os.path.exists(collated_path) and os.path.getsize(collated_path) > 0:
            logging.debug(f"[{marker_name}] Collated FASTA already exists -- skipping.")
            marker_to_fasta[marker_name] = collated_path
        else:
            pending_args.append((marker_name, sample_fastas, collated_path))

    if pending_args:
        logging.info(
            f"Concatenating sample FASTAs for {len(pending_args)} marker(s) "
            f"across {threads} thread(s) ..."
        )
        if threads > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                for marker_name, collated_path in executor.map(_cat_one_marker, pending_args):
                    marker_to_fasta[marker_name] = collated_path
        else:
            for marker_name, collated_path in map(_cat_one_marker, pending_args):
                marker_to_fasta[marker_name] = collated_path
        logging.info(f"Collation complete -- {len(pending_args)} marker(s) concatenated.")
    else:
        logging.info("All collated FASTA(s) already exist -- skipping concatenation.")

    logging.info(f"{len(marker_to_fasta)} marker(s) ready for clustering.")
    return marker_to_fasta


def _header_sort_key(header):
    """
    Reconstruct a sort key from the FASTA header written by otus_to_fasta.
    Header format: >otu{i}|{sample}|{gene}|unknown={0or1}|hits={n}
    Falls back gracefully if fields are missing.
    """
    # We embed sort fields in the header so we can recover them post-cat.
    # See otus_to_fasta() -- the enriched format is set there.
    try:
        parts = header.lstrip(">").split("|")
        known_rank = int(next((p.split("=")[1] for p in parts if p.startswith("unknown=")), "1"))
        hits = int(next((p.split("=")[1] for p in parts if p.startswith("hits=")), "0"))
        return (known_rank, -hits)
    except (ValueError, IndexError):
        return (1, 0)


# ---------------------------------------------------------------------------
# Phase 2: cluster locally (multiprocessing across markers)
# ---------------------------------------------------------------------------

def cluster_markers_locally(marker_to_fasta, output_dir, max_divergence, threads, chunk_size=2000000, cluster_threads=16):
    """Run smafa cluster for each marker in parallel (file-based, no fasta_str)."""

    def _cluster_one(args_tuple):
        (marker_name, fasta_path, max_divergence) = args_tuple
        cluster_tsv, rep_fasta = _cluster_marker_file_chunked(
            fasta_path, max_divergence, marker_name, chunk_size=chunk_size, threads=cluster_threads
        )
        n_clusters = len(set(
            ln.split("\t")[0] for ln in cluster_tsv.splitlines() if ln.strip()
        ))
        n_seqs = sum(1 for ln in cluster_tsv.splitlines() if ln.strip())
        logging.info(f"[{marker_name}] Done -- {n_clusters} cluster(s) from {n_seqs} OTU(s).")
        return (marker_name, cluster_tsv, rep_fasta)

    worker_args = [
        (marker, path, max_divergence)
        for marker, path in marker_to_fasta.items()
    ]

    logging.info(
        f"Clustering {len(worker_args)} marker(s) locally across {threads} thread(s)."
    )
    if threads > 1:
        with Pool(threads) as pool:
            results = pool.map(_cluster_one, worker_args)
    else:
        results = list(map(_cluster_one, worker_args))

    summary_rows = _write_cluster_results(results, marker_to_fasta, output_dir)
    _write_summary(summary_rows, output_dir)


def _write_cluster_results(results, marker_to_fasta, output_dir):
    """Write cluster TSV and rep FASTA; return summary rows (does not write summary.tsv)."""
    summary_rows = []
    for (marker_name, cluster_tsv, rep_fasta) in results:
        marker_dir = os.path.join(output_dir, marker_name)
        os.makedirs(marker_dir, exist_ok=True)

        with open(os.path.join(marker_dir, "clusters.tsv"), "w") as fh:
            fh.write(cluster_tsv)
        with open(os.path.join(marker_dir, "representatives.fasta"), "w") as fh:
            fh.write(rep_fasta)

        num_clusters = sum(1 for ln in cluster_tsv.splitlines() if ln.strip())
        fasta_path = marker_to_fasta.get(marker_name, "")
        num_otus = 0
        if fasta_path and os.path.exists(fasta_path):
            with open(fasta_path) as fh:
                num_otus = sum(1 for ln in fh if ln.startswith(">"))
        summary_rows.append({
            "marker": marker_name,
            "num_otus": num_otus,
            "num_clusters": num_clusters,
        })

    return summary_rows


def _write_summary(summary_rows, output_dir):
    """Write summary.tsv from a list of summary row dicts."""
    summary_path = os.path.join(output_dir, "summary.tsv")
    pd.DataFrame(summary_rows).sort_values("marker").to_csv(
        summary_path, sep="\t", index=False
    )
    logging.info(f"Summary written to {summary_path}")


# ---------------------------------------------------------------------------
# Phase 3 helpers: streaming FASTA splitter and per-chunk stats reader
# ---------------------------------------------------------------------------

def _split_fasta_to_chunk_files(fasta_path, chunk_dir, chunk_size):
    """
    Stream fasta_path and write chunk_{N:04d}.fasta files to chunk_dir.
    Returns the number of chunk files written (0 if fasta_path is empty).
    No sequences are kept in memory beyond a single entry at a time.
    """
    os.makedirs(chunk_dir, exist_ok=True)
    chunk_idx = 0
    seq_count_in_chunk = 0
    out_fh = None
    current_header = None
    seq_parts = []

    def _write_seq(header, parts):
        nonlocal out_fh, chunk_idx, seq_count_in_chunk
        seq = "".join(parts)
        if not seq:
            return
        if out_fh is None or seq_count_in_chunk >= chunk_size:
            if out_fh is not None:
                out_fh.close()
                chunk_idx += 1
            out_fh = open(os.path.join(chunk_dir, f"chunk_{chunk_idx:04d}.fasta"), "w")
            seq_count_in_chunk = 0
        out_fh.write(header + "\n" + seq + "\n")
        seq_count_in_chunk += 1

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None and seq_parts:
                    _write_seq(current_header, seq_parts)
                current_header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_header is not None and seq_parts:
        _write_seq(current_header, seq_parts)

    if out_fh is not None:
        out_fh.close()
        return chunk_idx + 1
    return 0


def _stream_chunk_seq_stats(chunk_fasta_path):
    """
    Stream a FASTA file and return {seq -> (unknown_rank, hits, n_occurrences, sum_hits)}
    for all sequences.  Only the four stat integers (not the header string) are
    stored per sequence.
    """
    seq_to_stats = {}
    current_header = None
    seq_parts = []
    with open(chunk_fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None and seq_parts:
                    seq = "".join(seq_parts)
                    seq_to_stats.setdefault(seq, _parse_header_stats(current_header))
                current_header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_header is not None and seq_parts:
        seq_to_stats.setdefault("".join(seq_parts), _parse_header_stats(current_header))
    return seq_to_stats


# ---------------------------------------------------------------------------
# Phase 3 (mqsub): submit one smafa cluster job per marker
# ---------------------------------------------------------------------------

def _already_clustered(output_dir, marker_name):
    """Return True if this marker's representatives.fasta already exists and is non-empty."""
    rep_path = os.path.join(output_dir, marker_name, "representatives.fasta")
    return os.path.exists(rep_path) and os.path.getsize(rep_path) > 0


def _write_chunk_at_offset(dedup_path, seq_start, n_seqs, out_path):
    """
    Stream dedup_path and write sequences [seq_start, seq_start+n_seqs) to out_path.
    Returns the actual number of sequences written.
    """
    written = 0
    seq_idx = 0
    cur_hdr = None
    seq_parts = []
    with open(dedup_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_hdr is not None and seq_parts:
                    if seq_start <= seq_idx < seq_start + n_seqs:
                        out_fh.write(cur_hdr + "\n" + "".join(seq_parts) + "\n")
                        written += 1
                    seq_idx += 1
                    if seq_idx >= seq_start + n_seqs:
                        return written
                cur_hdr = line
                seq_parts = []
            else:
                seq_parts.append(line)
        if cur_hdr is not None and seq_parts and seq_start <= seq_idx < seq_start + n_seqs:
            out_fh.write(cur_hdr + "\n" + "".join(seq_parts) + "\n")
            written += 1
    return written


def _rep_stats_path(output_dir, marker_name, round_idx):
    """Path for persisted rep-stats dict after completing round `round_idx`."""
    chunk_dir = os.path.join(output_dir, marker_name, "chunks")
    return os.path.join(chunk_dir, f"round_{round_idx:04d}_rep_stats.pkl")


def _save_rep_stats(path, rep_stats):
    """Atomically persist rep_stats dict to a pickle file."""
    import pickle
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(rep_stats, f, protocol=4)
    os.replace(tmp, path)


def _load_rep_stats(path):
    """Load a rep_stats dict from a pickle file."""
    import pickle
    with open(path, "rb") as f:
        return pickle.load(f)


def _cluster_file_status(round_clusters, round_input):
    """
    Returns:
      'complete' -- line count >= n_input_seqs and non-empty.
      'partial'  -- file exists and non-empty but fewer lines than expected.
      'missing'  -- file absent or empty.
    """
    if not os.path.exists(round_clusters) or not os.path.exists(round_input):
        return 'missing'
    r_c = subprocess.run(["wc", "-l", round_clusters], capture_output=True, text=True)
    r_i = subprocess.run(["grep", "-c", "^>", round_input], capture_output=True, text=True)
    if r_c.returncode != 0 or r_i.returncode not in (0, 1):
        return 'missing'
    try:
        n_cluster_lines = int(r_c.stdout.split()[0])
        n_input_seqs    = int(r_i.stdout.strip())
    except (ValueError, IndexError):
        return 'missing'
    if n_cluster_lines == 0:
        return 'missing'
    if n_cluster_lines >= n_input_seqs:
        return 'complete'
    return 'partial'


def _rescue_partial_cluster(round_clusters, round_input, initial_chunk, round_reps):
    """
    Rescue a partial cluster file (e.g. from a PBS job killed at the walltime limit):
      1. Strip the last (possibly truncated) line.
      2. Rebuild round_reps from the stripped cluster file + round_input sequences.
      3. Return the number of *new* sequences (from initial_chunk) that were rescued,
         or 0 if the file is not salvageable.
    """
    r_c  = subprocess.run(["wc", "-l", round_clusters],  capture_output=True, text=True)
    r_i  = subprocess.run(["grep", "-c", "^>", round_input],   capture_output=True, text=True)
    r_ch = subprocess.run(["grep", "-c", "^>", initial_chunk], capture_output=True, text=True)
    if r_c.returncode != 0 or r_i.returncode not in (0, 1) or r_ch.returncode not in (0, 1):
        return 0
    try:
        n_cluster_lines = int(r_c.stdout.split()[0])
        n_input_seqs    = int(r_i.stdout.strip())
        n_chunk_seqs    = int(r_ch.stdout.strip())
    except (ValueError, IndexError):
        return 0
    R_k = n_input_seqs - n_chunk_seqs   # prev-rep sequences at front of round_input
    rescued_lines = n_cluster_lines - 1  # drop last (possibly truncated) line
    if rescued_lines <= R_k:
        return 0   # not even the reps were fully processed
    # Overwrite cluster file with stripped version.
    tmp = round_clusters + ".rescue_tmp"
    ret = subprocess.run(
        ["head", "-n", str(rescued_lines), round_clusters],
        stdout=open(tmp, "wb"),
    )
    if ret.returncode != 0:
        if os.path.exists(tmp):
            os.unlink(tmp)
        return 0
    os.replace(tmp, round_clusters)
    # Rebuild reps.fasta from the stripped clusters + round_input.
    centroid_seqs = set()
    with open(round_clusters) as cf:
        for line in cf:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                centroid_seqs.add(parts[1])
    with open(round_reps, "w") as out_fh:
        cur_hdr, seq_parts = None, []
        with open(round_input) as fh:
            for raw in fh:
                ln = raw.rstrip("\n")
                if ln.startswith(">"):
                    if cur_hdr is not None and seq_parts:
                        seq = "".join(seq_parts)
                        if seq.replace("-", "N") in centroid_seqs:
                            out_fh.write(cur_hdr + "\n" + seq + "\n")
                    cur_hdr, seq_parts = ln, []
                else:
                    seq_parts.append(ln)
        if cur_hdr is not None and seq_parts:
            seq = "".join(seq_parts)
            if seq.replace("-", "N") in centroid_seqs:
                out_fh.write(cur_hdr + "\n" + seq + "\n")
    n_rescued = rescued_lines - R_k
    return n_rescued


def _adaptive_chunk_size(model_b, model_c, R_k, target_min, seqs_remaining):
    """
    Solve b*n^2 + c*R_k*n = target_min for n, where:
      n   = new sequences to add in this round,
      R_k = number of current representative sequences (cross-comparison term),
      b   = per-sequence quadratic cost coefficient (fitted from round 0),
      c   = per-rep per-new-seq cross-comparison cost coefficient (fitted from round 1+).

    Falls back to seqs_remaining if model_b is not yet fitted.
    When c is unknown but R_k > 0, uses a conservative estimate of c = 4 * model_b
    so the first cross-comparison round does not wildly overshoot the target.
    """
    if model_b is None:
        return seqs_remaining
    if R_k == 0:
        # No cross-comparison; pure quadratic: n = sqrt(target_min / b).
        n = int((target_min / model_b) ** 0.5)
    elif model_c is not None:
        # Full quadratic model: solve via quadratic formula.
        disc = (model_c * R_k) ** 2 + 4 * model_b * target_min
        n = int((-model_c * R_k + disc ** 0.5) / (2 * model_b))
    else:
        # c unknown but cross-comparison exists -- use conservative c = 4 * b.
        c_est = model_b * 4
        disc  = (c_est * R_k) ** 2 + 4 * model_b * target_min
        n     = int((-c_est * R_k + disc ** 0.5) / (2 * model_b))
    return max(1, min(n, seqs_remaining))


def cluster_markers_via_mqsub(
    marker_to_fasta, output_dir, max_divergence, this_script_path, cluster_memory=64,
    cluster_threads=16, chunk_size=2000000, target_walltime_hours=24.0,
):
    """
    Submit smafa cluster jobs via mqsub, one job per (marker, round) pair.

    All markers' round-N jobs are submitted together as a single mqsub batch,
    waited on, then round N+1 is submitted, etc.

    The orchestrator NEVER loads full sequence data into memory.  Instead:
      Phase A -- each marker's collated FASTA is deduplicated (disk sort).
      Phase B -- per round, prev_reps.fasta is catenated with a new chunk of
                 sequences; the chunk size adapts each round to keep per-job
                 walltime at target_walltime_hours (model fitted from PBS .OU logs).
      Phase C -- representatives.fasta is written from the final reps.fasta +
                 accumulated per-cluster stats loaded from disk.
    """
    abs_output_dir = os.path.abspath(output_dir)

    # Skip markers already clustered successfully.
    pending = {
        m: p for m, p in marker_to_fasta.items()
        if not _already_clustered(abs_output_dir, m)
    }
    skipped = len(marker_to_fasta) - len(pending)
    if skipped:
        logging.info(f"Skipping {skipped} already-clustered marker(s).")
    if not pending:
        logging.info("All markers already clustered.")
    else:
        target_min = target_walltime_hours * 60.0

        # ------------------------------------------------------------------
        # Phase A: deduplicate each collated FASTA (streaming, disk sort).
        # ------------------------------------------------------------------
        marker_dedup_fasta = {}   # marker -> path to deduplicated.fasta
        marker_total_seqs  = {}   # marker -> total unique seq count
        for marker_name, fasta_path in pending.items():
            marker_dir = os.path.join(abs_output_dir, marker_name)
            os.makedirs(marker_dir, exist_ok=True)
            dedup_path = os.path.join(marker_dir, "deduplicated.fasta")
            if not os.path.exists(dedup_path) or os.path.getsize(dedup_path) == 0:
                logging.info(
                    f"[{marker_name}] Deduplicating collated FASTA (aggregating counts) ..."
                )
                _dedup_fasta_with_counts(fasta_path, dedup_path)
                logging.info(f"[{marker_name}] Deduplication complete.")
            else:
                logging.info(
                    f"[{marker_name}] Deduplicated FASTA already exists -- skipping."
                )
            n_unique = sum(1 for ln in open(dedup_path) if ln.startswith(">"))
            marker_dedup_fasta[marker_name] = dedup_path
            marker_total_seqs[marker_name]  = n_unique
            logging.info(f"[{marker_name}] {n_unique} unique sequence(s) to cluster.")

        # Per-marker lightweight orchestrator state (no full stats dicts in memory
        # simultaneously).  Full stats are persisted to disk as pickle files and
        # loaded one marker at a time to keep peak memory to O(one_marker).
        marker_cluster_count      = {m: 0 for m in pending}
        marker_current_reps_fasta = {m: None for m in pending}

        # Track how many sequences have been consumed (clustered) per marker.
        marker_seqs_consumed = {m: 0 for m in pending}

        # Find the highest completed round index across all markers.
        max_existing_round = -1
        for marker_name in pending:
            chunk_dir = os.path.join(abs_output_dir, marker_name, "chunks")
            if not os.path.isdir(chunk_dir):
                continue
            for fname in os.listdir(chunk_dir):
                if fname.startswith("round_") and fname.endswith("_clusters.tsv"):
                    try:
                        idx = int(fname[len("round_"):len("round_") + 4])
                        max_existing_round = max(max_existing_round, idx)
                    except ValueError:
                        pass

        if max_existing_round >= 0:
            logging.info(
                f"Resuming -- replaying {max_existing_round + 1} completed round(s) "
                f"to rebuild orchestrator state ..."
            )
            # Process one marker at a time (all rounds) so only one marker's
            # rep-stats dict is in memory at once, reducing peak RAM dramatically.
            for marker_name in list(pending.keys()):
                persisted = _rep_stats_path(abs_output_dir, marker_name, max_existing_round)
                if os.path.exists(persisted):
                    # Fast path: load already-persisted stats from a previous run.
                    local_rep_stats = _load_rep_stats(persisted)
                    # Recover current_reps_fasta and seqs_consumed from disk files.
                    chunk_dir = os.path.join(abs_output_dir, marker_name, "chunks")
                    for ridx in range(max_existing_round + 1):
                        round_reps    = os.path.join(chunk_dir, f"round_{ridx:04d}_reps.fasta")
                        initial_chunk = os.path.join(chunk_dir, f"round_{ridx:04d}_chunk.fasta")
                        round_input_r = os.path.join(chunk_dir, f"round_{ridx:04d}_input.fasta")
                        round_clusters = os.path.join(chunk_dir, f"round_{ridx:04d}_clusters.tsv")
                        if not os.path.exists(round_reps) or not os.path.exists(initial_chunk):
                            continue
                        status = _cluster_file_status(round_clusters, round_input_r)
                        if status == 'missing':
                            continue
                        n_actual = sum(1 for ln in open(initial_chunk) if ln.startswith(">"))
                        marker_current_reps_fasta[marker_name] = round_reps
                        marker_seqs_consumed[marker_name] += n_actual
                    logging.info(
                        f"[{marker_name}] Loaded persisted stats "
                        f"({len(local_rep_stats)} cluster(s))."
                    )
                else:
                    # Slow path: replay all rounds for this marker to rebuild stats.
                    local_rep_stats = {}
                    chunk_dir = os.path.join(abs_output_dir, marker_name, "chunks")
                    for ridx in range(max_existing_round + 1):
                        initial_chunk  = os.path.join(chunk_dir, f"round_{ridx:04d}_chunk.fasta")
                        round_clusters = os.path.join(chunk_dir, f"round_{ridx:04d}_clusters.tsv")
                        round_reps     = os.path.join(chunk_dir, f"round_{ridx:04d}_reps.fasta")
                        round_input_r  = os.path.join(chunk_dir, f"round_{ridx:04d}_input.fasta")
                        status = _cluster_file_status(round_clusters, round_input_r)
                        if status == 'missing':
                            continue
                        if status == 'partial':
                            if not os.path.exists(initial_chunk):
                                continue
                            n_rescued = _rescue_partial_cluster(
                                round_clusters, round_input_r, initial_chunk, round_reps
                            )
                            if n_rescued == 0:
                                logging.warning(
                                    f"[{marker_name}] Replay round {ridx}: "
                                    f"partial cluster file, nothing rescuable -- skipping."
                                )
                                continue
                            logging.info(
                                f"[{marker_name}] Replay round {ridx}: "
                                f"rescued {n_rescued} seq(s) from partial cluster file."
                            )
                        if not os.path.exists(initial_chunk):
                            continue

                        old_rep_stats = local_rep_stats
                        chunk_seq_to_stats = _stream_chunk_seq_stats(initial_chunk)
                        new_rep_stats = {}
                        with open(round_clusters) as cf:
                            for line in cf:
                                parts = line.rstrip("\n").split("\t")
                                if len(parts) < 2:
                                    continue
                                member_seq, rep_seq = parts[0], parts[1]
                                if rep_seq not in new_rep_stats:
                                    new_rep_stats[rep_seq] = {
                                        "n_seqs": 0, "n_query": 0,
                                        "max_hits": 0, "total_hits": 0,
                                    }
                                s = new_rep_stats[rep_seq]
                                if member_seq in old_rep_stats:
                                    # member was a rep in a previous round; propagate its stats
                                    ms = old_rep_stats[member_seq]
                                    s["n_seqs"]     += ms["n_seqs"]
                                    s["n_query"]    += ms["n_query"]
                                    s["max_hits"]    = max(s["max_hits"], ms["max_hits"])
                                    s["total_hits"] += ms["total_hits"]
                                else:
                                    # member is a new sequence from this chunk
                                    unk, hits, n_occ, sum_hits_all = chunk_seq_to_stats.get(member_seq, (1, 0, 1, 0))
                                    s["n_seqs"]     += n_occ
                                    s["n_query"]    += (n_occ if unk == 0 else 0)
                                    s["max_hits"]    = max(s["max_hits"], hits)
                                    s["total_hits"] += sum_hits_all

                        del old_rep_stats, chunk_seq_to_stats  # free memory promptly
                        n_actual = sum(1 for ln in open(initial_chunk) if ln.startswith(">"))
                        if status == 'partial':
                            r_c2 = subprocess.run(["wc", "-l", round_clusters], capture_output=True, text=True)
                            r_i2 = subprocess.run(["grep", "-c", "^>", round_input_r], capture_output=True, text=True)
                            try:
                                n_actual = int(r_c2.stdout.split()[0]) - (int(r_i2.stdout.strip()) - n_actual)
                            except (ValueError, IndexError):
                                pass
                        local_rep_stats = new_rep_stats
                        del new_rep_stats  # keep only the latest round's stats
                        marker_current_reps_fasta[marker_name] = round_reps
                        marker_seqs_consumed[marker_name]      += n_actual

                    # Persist to disk so future runs skip replay for this marker.
                    if local_rep_stats:
                        _save_rep_stats(persisted, local_rep_stats)

                marker_cluster_count[marker_name] = len(local_rep_stats)
                del local_rep_stats  # free: will be reloaded from disk when needed in Phase B
                logging.info(
                    f"[{marker_name}] Resumed: {marker_seqs_consumed[marker_name]}/"
                    f"{marker_total_seqs[marker_name]} seq(s) consumed, "
                    f"{marker_cluster_count[marker_name]} cluster(s)."
                )

        # Adaptive timing model shared across markers; fitted from PBS .OU walltimes.
        # Models per-job walltime as: t ~= model_b * n^2 + model_c * R_k * n
        # where n = new sequences in this round and R_k = current rep count.
        model_file = os.path.join(abs_output_dir, "timing_model.json")
        model_b: float | None = None
        model_c: float | None = None
        if os.path.exists(model_file):
            with open(model_file) as _mf:
                _md = json.load(_mf)
            model_b = _md.get("model_b") or None
            model_c = _md.get("model_c") or None
            logging.info(
                f"Loaded timing model from {model_file}: "
                f"b={model_b}, c={model_c}"
            )

        # ------------------------------------------------------------------
        # Phase B: round-by-round mqsub submission with adaptive chunk sizes.
        # ------------------------------------------------------------------
        round_idx = max_existing_round + 1
        while any(marker_seqs_consumed[m] < marker_total_seqs[m] for m in pending):
            round_jobs = []  # (marker, initial_chunk, round_input, round_clusters, round_reps, n_k)
            for marker_name in list(pending.keys()):
                seqs_remaining = marker_total_seqs[marker_name] - marker_seqs_consumed[marker_name]
                if seqs_remaining <= 0:
                    continue

                R_k   = marker_cluster_count[marker_name]
                n_k   = _adaptive_chunk_size(model_b, model_c, R_k, target_min, seqs_remaining)
                # For round 0, honour the user-supplied chunk_size as initial cap.
                if round_idx == 0:
                    n_k = min(n_k, chunk_size)

                chunk_dir      = os.path.join(abs_output_dir, marker_name, "chunks")
                os.makedirs(chunk_dir, exist_ok=True)
                initial_chunk  = os.path.join(chunk_dir, f"round_{round_idx:04d}_chunk.fasta")
                round_input    = os.path.join(chunk_dir, f"round_{round_idx:04d}_input.fasta")
                round_clusters = os.path.join(chunk_dir, f"round_{round_idx:04d}_clusters.tsv")
                round_reps     = os.path.join(chunk_dir, f"round_{round_idx:04d}_reps.fasta")

                # Write the new-sequences chunk for this round.
                if not os.path.exists(initial_chunk) or os.path.getsize(initial_chunk) == 0:
                    n_written = _write_chunk_at_offset(
                        marker_dedup_fasta[marker_name],
                        marker_seqs_consumed[marker_name],
                        n_k, initial_chunk,
                    )
                    logging.info(
                        f"[{marker_name}] Round {round_idx + 1}: "
                        f"chunk of {n_written} seq(s) (R={R_k} reps, target={target_walltime_hours}h)."
                    )
                # Build round_input = prev_reps + initial_chunk (prev_reps omitted in round 0).
                prev_reps = marker_current_reps_fasta[marker_name]
                sources   = ([prev_reps] if prev_reps else []) + [initial_chunk]
                with open(round_input, "w") as out_fh:
                    for src in sources:
                        with open(src) as in_fh:
                            shutil.copyfileobj(in_fh, out_fh)

                # Count actual seqs in chunk (may differ from n_k at EOF).
                n_actual = sum(1 for ln in open(initial_chunk) if ln.startswith(">"))
                round_jobs.append(
                    (marker_name, initial_chunk, round_input, round_clusters, round_reps, n_actual)
                )

            if not round_jobs:
                break

            # Submit only jobs whose cluster output is missing (partial files are rescued below).
            jobs_to_submit = [
                (mn, ic, ri, rc, rr) for (mn, ic, ri, rc, rr, _) in round_jobs
                if _cluster_file_status(rc, ri) == 'missing'
            ]

            if jobs_to_submit:
                logging.info(
                    f"Round {round_idx + 1}: submitting {len(jobs_to_submit)} smafa job(s) via mqsub ..."
                )
                with tempfile.NamedTemporaryFile(
                    mode="w", prefix=f"smafa_chunk{round_idx:04d}_mqsub_",
                    suffix=".cmds", delete=False,
                ) as cmd_file:
                    cmd_file_path = cmd_file.name
                    for _mn, _ic, round_input, round_clusters, round_reps in jobs_to_submit:
                        cmd = (
                            f"python3 {this_script_path}"
                            f" --_cluster-chunk-fasta {round_input}"
                            f" --_chunk-clusters-tsv {round_clusters}"
                            f" --_chunk-reps-fasta {round_reps}"
                            f" --max-divergence {max_divergence}"
                            f" --_cluster-threads {cluster_threads}"
                        )
                        cmd_file.write(cmd + "\n")
                try:
                    mqsub_cmd = (
                        f"mqsub -m {cluster_memory} -t {cluster_threads}"
                        f" --name smafa_chunk_{round_idx:04d}"
                        f" --segregated-log-files --hours 48"
                        f" --command-file {cmd_file_path}"
                        f" --chunk-size 1 2>&1"
                    )
                    logging.info(f"Running: {mqsub_cmd}")
                    mqsub_stdout = extern.run(mqsub_cmd)
                    job_ids = _mqwait(mqsub_stdout)
                    # Build job_id -> marker_name map (submission order matches return order).
                    _marker_job_id = {mn: job_ids[i] for i, (mn, *_) in enumerate(jobs_to_submit)}
                finally:
                    os.unlink(cmd_file_path)
            else:
                _marker_job_id = {}

            logging.info(f"Round {round_idx + 1} mqsub jobs complete.")

            # Propagate state: stream clusters.tsv to update per-cluster stats.
            failed = []
            for marker_name, initial_chunk, round_input, round_clusters, round_reps, n_actual in round_jobs:
                status = _cluster_file_status(round_clusters, round_input)
                if status == 'partial':
                    n_rescued = _rescue_partial_cluster(
                        round_clusters, round_input, initial_chunk, round_reps
                    )
                    if n_rescued == 0:
                        logging.warning(
                            f"[{marker_name}] Round {round_idx + 1}: cluster file partial, "
                            f"nothing rescuable -- will retry."
                        )
                        failed.append(marker_name)
                        continue
                    logging.info(
                        f"[{marker_name}] Round {round_idx + 1}: rescued {n_rescued}/{n_actual} "
                        f"seq(s) from partial cluster file -- remainder deferred to next round."
                    )
                    n_actual = n_rescued  # propagate only the rescued portion
                elif status == 'missing':
                    continue

                # Load this marker's accumulated stats from disk (only one marker
                # in memory at a time to keep peak RAM at O(one_marker)).
                prev_stats_path = _rep_stats_path(abs_output_dir, marker_name, round_idx - 1)
                if os.path.exists(prev_stats_path):
                    old_rep_stats = _load_rep_stats(prev_stats_path)
                else:
                    old_rep_stats = {}

                # Stats for brand-new sequences in this chunk (O(chunk_size) memory).
                chunk_seq_to_stats = _stream_chunk_seq_stats(initial_chunk)

                new_rep_stats = {}
                with open(round_clusters) as cf:
                    for line in cf:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 2:
                            continue
                        member_seq = parts[0]
                        rep_seq    = parts[1]
                        if rep_seq not in new_rep_stats:
                            new_rep_stats[rep_seq] = {
                                "n_seqs": 0, "n_query": 0,
                                "max_hits": 0, "total_hits": 0,
                            }
                        s = new_rep_stats[rep_seq]
                        if member_seq in old_rep_stats:
                            # member was a rep in a previous round; propagate its stats
                            ms = old_rep_stats[member_seq]
                            s["n_seqs"]     += ms["n_seqs"]
                            s["n_query"]    += ms["n_query"]
                            s["max_hits"]    = max(s["max_hits"], ms["max_hits"])
                            s["total_hits"] += ms["total_hits"]
                        else:
                            # member is a new sequence from this chunk
                            unk, hits, n_occ, sum_hits_all = chunk_seq_to_stats.get(member_seq, (1, 0, 1, 0))
                            s["n_seqs"]     += n_occ
                            s["n_query"]    += (n_occ if unk == 0 else 0)
                            s["max_hits"]    = max(s["max_hits"], hits)
                            s["total_hits"] += sum_hits_all

                del old_rep_stats, chunk_seq_to_stats  # free before persisting

                # Persist updated stats to disk; free from memory.
                new_stats_path = _rep_stats_path(abs_output_dir, marker_name, round_idx)
                _save_rep_stats(new_stats_path, new_rep_stats)
                marker_cluster_count[marker_name]      = len(new_rep_stats)
                del new_rep_stats

                marker_current_reps_fasta[marker_name] = round_reps
                marker_seqs_consumed[marker_name]      += n_actual
                logging.info(
                    f"[{marker_name}] After round {round_idx + 1}: "
                    f"{marker_cluster_count[marker_name]} cluster(s), "
                    f"{marker_seqs_consumed[marker_name]}/{marker_total_seqs[marker_name]} seq(s) done."
                )

            if failed:
                logging.warning(
                    f"{len(failed)} smafa chunk job(s) produced incomplete output in round "
                    f"{round_idx + 1} -- will retry next round: "
                    + ", ".join(sorted(failed))
                )

            # Update timing model using walltimes only from jobs that succeeded.
            successful_markers = set(mn for mn, *_ in round_jobs) - set(failed)
            successful_job_ids = [
                _marker_job_id[mn]
                for mn in (mn for mn, *_ in jobs_to_submit if mn in _marker_job_id)
                if mn in successful_markers
            ]
            walltimes = _parse_walltimes_from_job_ids(successful_job_ids)
            if walltimes:
                wt_med = sorted(walltimes)[len(walltimes) // 2]
                n_vals = [n for mn, _, _, _, _, n in round_jobs if mn in successful_markers]
                n_med  = sorted(n_vals)[len(n_vals) // 2] if n_vals else chunk_size
                R_vals = [marker_cluster_count[mn] for mn, *_ in round_jobs
                          if mn in successful_markers]
                R_med  = sorted(R_vals)[len(R_vals) // 2] if R_vals else 0
                if round_idx == 0:
                    # Fit b from round-0 jobs: t0 = b * n^2  ->  b = t0 / n^2
                    model_b = wt_med / (n_med ** 2) if n_med > 0 else None
                    logging.info(
                        f"Timing model: b={model_b:.3e} min/seq^2 "
                        f"(from round 0 median {wt_med:.1f} min, n={n_med}, "
                        f"{len(walltimes)}/{len(round_jobs)} succeeded)"
                    )
                elif model_b is not None and R_med > 0:
                    # Fit c from round 1+ jobs: t1 = b*n^2 + c*R*n  ->  c = (t - b*n^2) / (R*n)
                    cross = wt_med - model_b * n_med ** 2
                    if cross > 0:
                        model_c = cross / (R_med * n_med)
                        logging.info(
                            f"Timing model: b={model_b:.3e}, c={model_c:.3e} min/rep/seq "
                            f"(from round {round_idx + 1} median {wt_med:.1f} min, "
                            f"n={n_med}, R={R_med}, {len(walltimes)}/{len(round_jobs)} succeeded)"
                        )
                # Persist updated model so resuming runs inherit it.
                with open(model_file, "w") as _mf:
                    json.dump({"model_b": model_b, "model_c": model_c}, _mf)

            round_idx += 1


        # ------------------------------------------------------------------
        # Phase C: write representatives.fasta per marker from accumulated stats.
        # ------------------------------------------------------------------
        for marker_name in list(pending.keys()):
            fasta_path = pending[marker_name]
            reps_fasta = marker_current_reps_fasta[marker_name]

            marker_dir = os.path.join(abs_output_dir, marker_name)
            os.makedirs(marker_dir, exist_ok=True)

            # Load this marker's accumulated stats from the latest persisted pickle.
            chunk_dir = os.path.join(abs_output_dir, marker_name, "chunks")
            stats_files = sorted(
                f for f in os.listdir(chunk_dir) if f.endswith("_rep_stats.pkl")
            ) if os.path.isdir(chunk_dir) else []
            if not stats_files:
                logging.warning(f"[{marker_name}] No persisted rep stats found -- skipping.")
                continue
            rep_stats = _load_rep_stats(os.path.join(chunk_dir, stats_files[-1]))

            # Read the final round's reps.fasta to get original headers for each rep seq.
            rep_to_header = {}
            if reps_fasta and os.path.exists(reps_fasta):
                cur_hdr = None
                seq_parts = []
                with open(reps_fasta) as fh:
                    for line in fh:
                        line = line.rstrip("\n")
                        if line.startswith(">"):
                            if cur_hdr is not None and seq_parts:
                                rep_to_header["".join(seq_parts).replace('-', 'N')] = cur_hdr
                            cur_hdr = line
                            seq_parts = []
                        else:
                            seq_parts.append(line)
                if cur_hdr is not None and seq_parts:
                    rep_to_header["".join(seq_parts).replace('-', 'N')] = cur_hdr

            # Write representatives.fasta with enriched headers.
            rep_fasta_path = os.path.join(marker_dir, "representatives.fasta")
            n_clusters = 0
            with open(rep_fasta_path, "w") as out_fh:
                for rep_seq, stats in rep_stats.items():
                    orig_header = rep_to_header.get(rep_seq)
                    if orig_header is None:
                        logging.warning(
                            f"[{marker_name}] Rep sequence not found in reps.fasta -- skipping."
                        )
                        continue
                    enriched = (
                        f"{orig_header}"
                        f"|cluster_size={stats['n_seqs']}"
                        f"|n_query={stats['n_query']}"
                        f"|max_hits={stats['max_hits']}"
                        f"|total_hits={stats['total_hits']}"
                    )
                    out_fh.write(enriched + "\n" + rep_seq + "\n")
                    n_clusters += 1

            del rep_stats  # free after writing

            n_otus = 0
            if os.path.exists(fasta_path):
                with open(fasta_path) as fh:
                    n_otus = sum(1 for ln in fh if ln.startswith(">"))
            logging.info(
                f"[{marker_name}] Done -- {n_clusters} cluster(s) from {n_otus} OTU(s)."
            )

        logging.info("All smafa cluster jobs finished.")

    # Aggregate summary across all marker dirs (read from representatives.fasta).
    summary_rows = []
    for marker_name, fasta_path in marker_to_fasta.items():
        marker_dir     = os.path.join(abs_output_dir, marker_name)
        rep_fasta_path = os.path.join(marker_dir, "representatives.fasta")
        num_clusters = 0
        if os.path.exists(rep_fasta_path):
            with open(rep_fasta_path) as fh:
                num_clusters = sum(1 for ln in fh if ln.startswith(">"))
        num_otus = 0
        if os.path.exists(fasta_path):
            with open(fasta_path) as fh:
                num_otus = sum(1 for ln in fh if ln.startswith(">"))
        summary_rows.append({
            "marker": marker_name,
            "num_otus": num_otus,
            "num_clusters": num_clusters,
        })
    _write_summary(summary_rows, abs_output_dir)


# ---------------------------------------------------------------------------
# Shared mqsub wait helper
# ---------------------------------------------------------------------------

def _mqwait(mqsub_log):
    """Wait for all job IDs parsed from mqsub stdout. Returns list of job ID strings."""
    r = re.compile(r'^qsub stdout: (\d+\.aqua)$')
    job_ids = []
    for line in mqsub_log.split('\n'):
        m = r.match(line)
        if m:
            job_ids.append(m.group(1))
    logging.info(f"Waiting for {len(job_ids)} job(s) to finish ...")
    if not job_ids:
        return job_ids
    with tempfile.NamedTemporaryFile(mode="w", prefix="mqwait_", suffix=".ids") as f:
        f.write('\n'.join(job_ids) + '\n')
        f.flush()
        extern.run(f"mqwait -i {f.name}")
    return job_ids


def _parse_walltimes_from_job_ids(job_ids):
    """
    Given a list of PBS job IDs (e.g. ['21290083.aqua', ...]), locate each
    job's .OU file via `qstat -fx <jobid>` and parse the Wall time field.
    Returns a list of wall-clock times in minutes.
    """
    import subprocess
    walltimes = []
    for job_id in job_ids:
        try:
            result = subprocess.run(
                ["qstat", "-fx", job_id],
                capture_output=True, text=True, timeout=30,
            )
            # Extract Output_Path (may be line-wrapped with a tab continuation).
            output = result.stdout.replace('\n\t', '').replace('\t', '')
            ou_dir = None
            for line in output.splitlines():
                if 'Output_Path' in line and ':' in line:
                    ou_dir = line.split(':', 1)[1].strip()
                    break
            if ou_dir is None:
                continue
            import glob
            for ou_file in glob.glob(ou_dir + '/*.OU'):
                try:
                    with open(ou_file) as fh:
                        for line in fh:
                            if 'Wall time' in line and ':' in line:
                                time_part = line.split(': ', 1)[1].strip()
                                parts = time_part.split(':')
                                if len(parts) == 3:
                                    h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
                                    walltimes.append(h * 60 + m + s / 60)
                                    break
                except OSError:
                    pass
        except (subprocess.TimeoutExpired, OSError, ValueError):
            pass
    return walltimes


# ---------------------------------------------------------------------------
# mqsub worker entry points
# ---------------------------------------------------------------------------

def worker_extract_sample(archive_path, marker_domains_tsv, sample_fasta_dir, markers_of_interest):
    """
    Entry point for per-sample extraction mqsub jobs.
    Reads marker_domains from the TSV written by the orchestrator.
    """
    marker_domains = {}
    with open(marker_domains_tsv) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            marker, domains_str = line.split("\t", 1)
            marker_domains[marker] = set(domains_str.split(","))

    seen = extract_sample_fastas(
        archive_path, marker_domains, markers_of_interest, sample_fasta_dir
    )
    logging.info(f"Extracted {len(seen)} marker(s) from {archive_path}")


def worker_cluster_chunk(input_fasta, output_clusters_tsv, output_reps_fasta, max_divergence, threads=16):
    """Entry point for a single-chunk smafa cluster mqsub job.

    Runs smafa on input_fasta, writes clusters.tsv and (if output_reps_fasta is
    given) a reps.fasta containing the header+sequence for each unique centroid.
    """
    if not os.path.exists(input_fasta) or os.path.getsize(input_fasta) == 0:
        logging.warning(f"Empty or missing chunk FASTA: {input_fasta} -- writing empty output.")
        open(output_clusters_tsv, "w").close()
        if output_reps_fasta:
            open(output_reps_fasta, "w").close()
        return
    n_seqs = sum(1 for ln in open(input_fasta) if ln.startswith(">"))
    logging.info(f"Running smafa on {n_seqs} sequence(s) from {input_fasta} ...")
    cluster_tmp = _run_smafa_cluster_file(input_fasta, max_divergence, "chunk", threads=threads)
    try:
        shutil.move(cluster_tmp, output_clusters_tsv)
    except Exception:
        if os.path.exists(cluster_tmp):
            os.unlink(cluster_tmp)
        raise
    centroid_seqs = set()
    with open(output_clusters_tsv) as cf:
        for line in cf:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                centroid_seqs.add(parts[1])
    n_clusters = len(centroid_seqs)
    logging.info(f"Chunk done -- {n_clusters} cluster(s) from {n_seqs} sequence(s).")

    # Write reps.fasta: scan input_fasta for sequences that are centroids.
    if output_reps_fasta:
        with open(output_reps_fasta, "w") as out_fh:
            cur_hdr = None
            seq_parts = []
            with open(input_fasta) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if cur_hdr is not None and seq_parts:
                            seq = "".join(seq_parts)
                            if seq.replace('-', 'N') in centroid_seqs:
                                out_fh.write(cur_hdr + "\n" + seq + "\n")
                        cur_hdr = line
                        seq_parts = []
                    else:
                        seq_parts.append(line)
            if cur_hdr is not None and seq_parts:
                seq = "".join(seq_parts)
                if seq.replace('-', 'N') in centroid_seqs:
                    out_fh.write(cur_hdr + "\n" + seq + "\n")


def worker_cluster_marker(fasta_path, marker_name, output_dir, max_divergence, chunk_size=2000000, cluster_threads=16):
    """Entry point for per-marker smafa cluster mqsub jobs (file-based, no fasta_str)."""
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        logging.warning(f"[{marker_name}] Empty or missing FASTA -- skipping.")
        return
    cluster_tsv, rep_fasta = _cluster_marker_file_chunked(
        fasta_path, max_divergence, marker_name, chunk_size=chunk_size, threads=cluster_threads
    )
    result = (marker_name, cluster_tsv, rep_fasta)
    _write_cluster_results([result], {marker_name: fasta_path}, output_dir)
    n_clusters = len(set(
        ln.split("\t")[0] for ln in cluster_tsv.splitlines() if ln.strip()
    ))
    n_seqs = sum(1 for ln in cluster_tsv.splitlines() if ln.strip())
    logging.info(f"[{marker_name}] Done -- {n_clusters} cluster(s) from {n_seqs} OTU(s).")
    logging.info(f"[{marker_name}] Cluster job complete.")  # summary.tsv written by orchestrator


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def collate_and_cluster(
    archive_paths,
    metapackage_path,
    output_dir,
    max_divergence,
    threads,
    markers_of_interest,
    run_through_mqsub,
    this_script_path,
    cluster_memory=64,
    cluster_threads=16,
    chunk_size=2000000,
):
    os.makedirs(output_dir, exist_ok=True)

    # Build marker->domain map once (used by all paths)
    if metapackage_path:
        logging.info("Loading metapackage for domain filtering ...")
        marker_domains = build_marker_domain_map(metapackage_path)
        # Restrict to markers that target both Bacteria and Archaea
        both_domains = {
            m for m, domains in marker_domains.items()
            if "Bacteria" in domains and "Archaea" in domains
        }
        logging.info(
            f"{len(both_domains)} marker(s) target both Bacteria and Archaea."
        )
        markers_of_interest = (
            markers_of_interest & both_domains if markers_of_interest else both_domains
        )
    else:
        marker_domains = {}
        logging.warning("No metapackage provided -- skipping off-target filtering.")

    if run_through_mqsub:
        # Phase 1: submit one extraction job per archive (skip already-done ones)
        phase_1_sentinel = os.path.join(output_dir, SAMPLE_FASTA_SUBDIR, "done")
        if not os.path.exists(phase_1_sentinel):
            pending = [p for p in archive_paths if not _already_extracted(output_dir, p)]
            skipped = len(archive_paths) - len(pending)
            if skipped:
                logging.info(f"Skipping {skipped} already-extracted archive(s).")
            submit_extraction_jobs(
                pending, marker_domains, markers_of_interest,
                output_dir, this_script_path
            )
            open(phase_1_sentinel, "w").close()
        else:
            logging.info("Skipping Phase 1 extraction -- already done.")

        # Phase 2: cat per-sample FASTAs into per-marker FASTAs (local, fast)
        marker_to_fasta = cat_sample_fastas(archive_paths, output_dir, markers_of_interest, threads)
        if not marker_to_fasta:
            logging.warning("No markers with sequences -- nothing to cluster.")
            return

        # Phase 3: submit one smafa cluster job per (marker, round) pair via mqsub
        cluster_markers_via_mqsub(
            marker_to_fasta, output_dir, max_divergence, this_script_path,
            cluster_memory=cluster_memory,
            cluster_threads=cluster_threads,
            chunk_size=chunk_size,
            target_walltime_hours=args.target_walltime_hours,
        )
    else:
        # Local path: extract all samples in parallel, then cluster in parallel
        # Phase 1: extract all samples in parallel (skip already-done ones)
        phase_1_sentinel = os.path.join(output_dir, SAMPLE_FASTA_SUBDIR, "done")
        if not os.path.exists(phase_1_sentinel):
            pending = [p for p in archive_paths if not _already_extracted(output_dir, p)]
            skipped = len(archive_paths) - len(pending)
            if skipped:
                logging.info(f"Skipping {skipped} already-extracted archive(s).")
            logging.info(
                f"Extracting {len(pending)} archive(s) locally "
                f"across {threads} thread(s) ..."
            )
            worker_args = [
                (
                    archive_path,
                    marker_domains,
                    markers_of_interest,
                    _sample_fasta_dir(output_dir, archive_path),
                )
                for archive_path in pending
            ]
            if threads > 1:
                with Pool(threads) as pool:
                    pool.map(extract_sample_fastas_local_worker, worker_args)
            else:
                list(map(extract_sample_fastas_local_worker, worker_args))
            open(phase_1_sentinel, "w").close()
        else:
            logging.info("Skipping Phase 1 extraction -- already done.")

        # Phase 2: cat per-sample FASTAs into per-marker FASTAs (local, fast)
        marker_to_fasta = cat_sample_fastas(archive_paths, output_dir, markers_of_interest, threads)
        if not marker_to_fasta:
            logging.warning("No markers with sequences -- nothing to cluster.")
            return

        # Phase 3: smafa cluster
        cluster_markers_locally(
            marker_to_fasta, output_dir, max_divergence, threads,
            chunk_size=chunk_size, cluster_threads=cluster_threads
        )

    logging.info("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    _this_script = os.path.abspath(__file__)

    parser = argparse.ArgumentParser(
        description=(
            "Collate OTU sequences from SingleM archive OTU tables per marker, "
            "then cluster with smafa."
        )
    )
    parser.add_argument("--debug", action="store_true", help="Output debug information")
    parser.add_argument("--quiet", action="store_true", help="Only output errors")

    # --- Normal invocation inputs ---
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--input-archive-otu-tables",
        nargs="+",
        metavar="FILE",
        help="One or more archive OTU table files (.json or .json.gz)",
    )
    input_group.add_argument(
        "--input-archive-otu-table-list",
        metavar="FILE",
        help="File listing archive OTU table paths, one per line",
    )

    parser.add_argument(
        "--metapackage",
        metavar="PATH",
        help=(
            "SingleM metapackage used to determine target domains per marker "
            "(required for off-target filtering; omit to skip filtering)"
        ),
    )
    parser.add_argument(
        "--output-directory",
        required=False,
        metavar="DIR",
        help="Directory to write per-marker cluster outputs",
    )
    parser.add_argument(
        "--max-divergence",
        type=int,
        default=2,
        metavar="N",
        help="smafa cluster --max-divergence value (default: 2)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="N",
        help="Number of markers to process in parallel when not using mqsub (default: 1)",
    )
    parser.add_argument(
        "--markers",
        nargs="+",
        metavar="MARKER",
        help="Restrict processing to these marker names (default: all markers)",
    )
    parser.add_argument(
        "--run-through-mqsub",
        action="store_true",
        help=(
            "Submit smafa cluster jobs to the HPC queue via mqsub. "
            "Phase 1 (collation) still runs locally; one job is submitted per marker."
        ),
    )
    parser.add_argument(
        "--cluster-memory",
        type=int,
        default=64,
        metavar="GB",
        help="Memory (GB) to request for each smafa cluster mqsub job (default: 64)",
    )

    parser.add_argument(
        "--cluster-threads",
        type=int,
        default=16,
        metavar="N",
        help="Threads for smafa per chunk (and per mqsub job) (default: 16)",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2000000,
        metavar="N",
        help="Max sequences per smafa chunk for round 0 (default: 2000000). "
             "Subsequent rounds are sized adaptively by --target-walltime-hours.",
    )

    parser.add_argument(
        "--target-walltime-hours",
        type=float,
        default=24.0,
        metavar="H",
        help="Target per-job walltime in hours; chunk size is adapted each round "
             "to hit this target (default: 24.0).",
    )

    # --- Internal: used only when this script is re-invoked by an mqsub worker ---
    parser.add_argument("--_cluster-chunk-fasta", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_chunk-clusters-tsv", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_chunk-reps-fasta", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_cluster-threads", type=int, default=16, help=argparse.SUPPRESS)
    # Kept for backward-compat with any running jobs; not used by new orchestrator
    parser.add_argument("--_cluster-marker-fasta", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_marker-name", metavar="NAME", help=argparse.SUPPRESS)
    parser.add_argument("--_chunk-size", type=int, default=2000000, help=argparse.SUPPRESS)
    # Per-sample extraction worker args
    parser.add_argument("--_extract-sample-archive", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_marker-domains-tsv", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_sample-fasta-dir", metavar="DIR", help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Logging setup
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # --- Worker path: per-sample extraction (invoked by mqsub) ---
    if args._extract_sample_archive:
        if not args._marker_domains_tsv or not args._sample_fasta_dir:
            parser.error(
                "--_extract-sample-archive requires --_marker-domains-tsv and --_sample-fasta-dir"
            )
        worker_extract_sample(
            archive_path=args._extract_sample_archive,
            marker_domains_tsv=args._marker_domains_tsv,
            sample_fasta_dir=args._sample_fasta_dir,
            markers_of_interest=set(args.markers) if args.markers else None,
        )
        sys.exit(0)

    # --- Worker path: single-chunk smafa cluster job (invoked by mqsub) ---
    if args._cluster_chunk_fasta:
        if not args._chunk_clusters_tsv:
            parser.error("--_cluster-chunk-fasta requires --_chunk-clusters-tsv")
        worker_cluster_chunk(
            input_fasta=args._cluster_chunk_fasta,
            output_clusters_tsv=args._chunk_clusters_tsv,
            output_reps_fasta=args._chunk_reps_fasta,
            max_divergence=args.max_divergence,
            threads=args._cluster_threads,
        )
        sys.exit(0)

    # --- Worker path: per-marker smafa cluster (backward-compat) ---
    if args._cluster_marker_fasta:
        if not args._marker_name:
            parser.error("--_cluster-marker-fasta requires --_marker-name")
        worker_cluster_marker(
            fasta_path=args._cluster_marker_fasta,
            marker_name=args._marker_name,
            output_dir=args.output_directory,
            max_divergence=args.max_divergence,
            chunk_size=args._chunk_size,
            cluster_threads=args._cluster_threads,
        )
        sys.exit(0)

    # --- Normal path ---
    if not args.output_directory:
        parser.error("--output-directory is required.")
    if not args.input_archive_otu_tables and not args.input_archive_otu_table_list:
        parser.error(
            "One of --input-archive-otu-tables or --input-archive-otu-table-list is required."
        )

    if args.input_archive_otu_tables:
        archive_paths = args.input_archive_otu_tables
    else:
        with open(args.input_archive_otu_table_list) as fh:
            archive_paths = [line.strip() for line in fh if line.strip()]

    collate_and_cluster(
        archive_paths=archive_paths,
        metapackage_path=args.metapackage,
        output_dir=args.output_directory,
        max_divergence=args.max_divergence,
        threads=args.threads,
        markers_of_interest=set(args.markers) if args.markers else None,
        run_through_mqsub=args.run_through_mqsub,
        this_script_path=_this_script,
        cluster_memory=args.cluster_memory,
        cluster_threads=args.cluster_threads,
        chunk_size=args.chunk_size,
    )
