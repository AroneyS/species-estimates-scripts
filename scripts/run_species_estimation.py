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

    'taxonomy_assignment_method' == 'singlem_query_based'  →  known (lower sort value = first).
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
    """Extract (unknown_rank, hits) from a FASTA header written by otus_to_fasta."""
    try:
        parts = header.lstrip(">").split("|")
        unknown_rank = int(next((p.split("=")[1] for p in parts if p.startswith("unknown=")), "1"))
        hits = int(next((p.split("=")[1] for p in parts if p.startswith("hits=")), "0"))
    except (ValueError, IndexError):
        unknown_rank, hits = 1, 0
    return unknown_rank, hits


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
      seq_to_stats:  seq -> (unknown_rank, hits)
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
    Read a smafa clusters TSV (member\tcentroid) and propagate accumulated
    cluster membership from previous rounds.

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
        n_seqs = len(members)
        n_query = sum(1 for m in members if seq_to_stats.get(m, (1, 0))[0] == 0)
        max_hits = max((seq_to_stats.get(m, (1, 0))[1] for m in members), default=0)
        total_hits = sum(seq_to_stats.get(m, (1, 0))[1] for m in members)
        orig_header = seq_to_header.get(rep_seq)
        if orig_header is None:
            logging.warning(f"[{marker_name}] Rep sequence not found in FASTA — skipping.")
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
    using three streaming passes.

    Pass 1 (cluster file):  collect rep_seqs set  (O(n_clusters) memory).
    Pass 2 (FASTA):         build rep_to_header   (O(n_clusters) — reps only)
                            build seq_to_stats    (O(n_unique_seqs) × 2 ints — no header strings).
    Pass 3 (cluster file):  aggregate per-cluster stats using seq_to_stats.
    """
    # Pass 1: get the set of representative sequences.
    # smafa writes  sequence\tcentroid  so col 2 is the centroid/rep.
    rep_seqs = set()
    with open(cluster_path) as cf:
        for line in cf:
            line = line.rstrip("\n")
            if line:
                parts = line.split("\t")
                if len(parts) >= 2:
                    rep_seqs.add(parts[1])

    # Pass 2: stream the FASTA once.
    #   - seq_to_stats: seq -> (unknown_rank, hits)  for every unique sequence
    #     (values are two ints, not header strings)
    #   - rep_to_header: seq -> header  only for sequences that are cluster reps
    seq_to_stats = {}   # seq -> (unknown_rank, hits)
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
            s["n_seqs"] += 1
            unk, hits = seq_to_stats.get(member_seq, (1, 0))
            s["n_query"] += 1 if unk == 0 else 0
            s["max_hits"] = max(s["max_hits"], hits)
            s["total_hits"] += hits

    # Build enriched rep FASTA.
    rep_fasta_lines = []
    for rep_seq, stats in rep_stats.items():
        orig_header = rep_to_header.get(rep_seq)
        if orig_header is None:
            logging.warning(f"[{marker_name}] Rep sequence not found in FASTA — skipping.")
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

    Sequences are NOT globally sorted here — sorting across samples requires
    the full dataset and is done after cat (phase 2).

    Returns the set of marker names that received ≥1 sequence.
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

    # Serialise marker_domains for worker jobs (pass metapackage path instead —
    # workers rebuild it themselves, which avoids pickling singlem objects).
    # Instead we write a simple marker→domains TSV that workers can read.
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
                f"All {len(existing)} collated FASTA(s) already exist — skipping concatenation."
            )
            return existing

    # Discover which markers need (re-)collation.
    marker_to_sample_fastas = defaultdict(list)
    for archive_path in archive_paths:
        sample_dir = _sample_fasta_dir(output_dir, archive_path)
        if not os.path.isdir(sample_dir):
            logging.warning(f"No output directory for {archive_path} — skipping in cat.")
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
            logging.debug(f"[{marker_name}] Collated FASTA already exists — skipping.")
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
        logging.info(f"Collation complete — {len(pending_args)} marker(s) concatenated.")
    else:
        logging.info("All collated FASTA(s) already exist — skipping concatenation.")

    logging.info(f"{len(marker_to_fasta)} marker(s) ready for clustering.")
    return marker_to_fasta


def _header_sort_key(header):
    """
    Reconstruct a sort key from the FASTA header written by otus_to_fasta.
    Header format: >otu{i}|{sample}|{gene}|unknown={0or1}|hits={n}
    Falls back gracefully if fields are missing.
    """
    # We embed sort fields in the header so we can recover them post-cat.
    # See otus_to_fasta() — the enriched format is set there.
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
        logging.info(f"[{marker_name}] Done — {n_clusters} cluster(s) from {n_seqs} OTU(s).")
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
    Stream a FASTA file and return {seq -> (unknown_rank, hits)} for all sequences.
    Only the two stat integers (not the header string) are stored per sequence.
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


def cluster_markers_via_mqsub(
    marker_to_fasta, output_dir, max_divergence, this_script_path, cluster_memory=64,
    cluster_threads=16, chunk_size=2000000,
):
    """
    Submit smafa cluster jobs via mqsub, one job per (marker, chunk) pair.

    All markers' chunk-N jobs are submitted together as a single mqsub batch,
    waited on, then chunk N+1 is submitted, etc.

    The orchestrator NEVER loads full sequence data into memory.  Instead:
      Phase A — each marker's collated FASTA is split into chunk files on disk
                 (streaming, one sequence at a time).
      Phase B — per round, prev_reps.fasta (small) is catenated with the next
                 chunk file; workers run smafa and write clusters.tsv +
                 reps.fasta; the orchestrator streams clusters.tsv to
                 accumulate per-cluster stat counters (not full member lists).
      Phase C — representatives.fasta is written from the final reps.fasta +
                 accumulated stats.  No full clusters.tsv is written for the
                 large-scale mqsub path (it would be tens of GB per marker).
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
        # ------------------------------------------------------------------
        # Phase A: pre-split collated FASTAs into chunk files (streaming).
        # One marker at a time keeps peak memory at chunk_size sequences.
        # ------------------------------------------------------------------
        n_chunks_per_marker = {}
        for marker_name, fasta_path in pending.items():
            chunk_dir = os.path.join(abs_output_dir, marker_name, "chunks")
            # Resume: count existing chunk files so we don't re-split.
            existing = sorted(
                f for f in (os.listdir(chunk_dir) if os.path.isdir(chunk_dir) else [])
                if f.startswith("chunk_") and f.endswith(".fasta") and "input" not in f
            )
            if existing:
                n_chunks_per_marker[marker_name] = len(existing)
                logging.info(
                    f"[{marker_name}] Resuming — found {len(existing)} existing chunk file(s)."
                )
            else:
                logging.info(
                    f"[{marker_name}] Splitting collated FASTA into chunk(s) of {chunk_size} ..."
                )
                n = _split_fasta_to_chunk_files(fasta_path, chunk_dir, chunk_size)
                n_chunks_per_marker[marker_name] = n
                logging.info(f"[{marker_name}] Split into {n} chunk(s).")

        max_rounds = max(n_chunks_per_marker.values(), default=0)
        logging.info(
            f"Clustering {len(pending)} marker(s) over {max_rounds} round(s) via mqsub "
            f"(streaming orchestrator — no sequences in orchestrator memory)."
        )

        # Per-marker in-memory state (only small stats dicts, no sequence strings).
        # marker_rep_stats:          {marker -> {rep_seq -> {n_seqs,n_query,max_hits,total_hits}}}
        # marker_current_reps_fasta: {marker -> path to most-recent round_*_reps.fasta}
        marker_rep_stats = {m: {} for m in pending}
        marker_current_reps_fasta = {m: None for m in pending}

        # ------------------------------------------------------------------
        # Phase B: round-by-round mqsub submission.
        # ------------------------------------------------------------------
        for chunk_idx in range(max_rounds):
            # Build the list of (marker, paths) for this round.
            round_jobs = []  # (marker_name, initial_chunk, round_input, round_clusters, round_reps)
            for marker_name in list(pending.keys()):
                if chunk_idx >= n_chunks_per_marker[marker_name]:
                    continue
                chunk_dir      = os.path.join(abs_output_dir, marker_name, "chunks")
                initial_chunk  = os.path.join(chunk_dir, f"chunk_{chunk_idx:04d}.fasta")
                round_input    = os.path.join(chunk_dir, f"round_{chunk_idx:04d}_input.fasta")
                round_clusters = os.path.join(chunk_dir, f"round_{chunk_idx:04d}_clusters.tsv")
                round_reps     = os.path.join(chunk_dir, f"round_{chunk_idx:04d}_reps.fasta")

                # Resume: if the clusters file already exists this round is done.
                if os.path.exists(round_clusters) and os.path.getsize(round_clusters) > 0:
                    logging.info(
                        f"[{marker_name}] Round {chunk_idx + 1} already done — will propagate."
                    )
                else:
                    # Build round_input by streaming prev_reps + initial_chunk.
                    prev_reps = marker_current_reps_fasta[marker_name]
                    sources = ([prev_reps] if prev_reps else []) + [initial_chunk]
                    with open(round_input, "w") as out_fh:
                        for src in sources:
                            with open(src) as in_fh:
                                shutil.copyfileobj(in_fh, out_fh)

                round_jobs.append(
                    (marker_name, initial_chunk, round_input, round_clusters, round_reps)
                )

            if not round_jobs:
                continue

            # Submit only the jobs whose cluster output doesn't exist yet.
            jobs_to_submit = [
                (mn, ic, ri, rc, rr) for (mn, ic, ri, rc, rr) in round_jobs
                if not (os.path.exists(rc) and os.path.getsize(rc) > 0)
            ]

            if jobs_to_submit:
                logging.info(
                    f"Round {chunk_idx + 1}/{max_rounds}: submitting "
                    f"{len(jobs_to_submit)} smafa chunk job(s) via mqsub ..."
                )
                with tempfile.NamedTemporaryFile(
                    mode="w", prefix=f"smafa_chunk{chunk_idx:04d}_mqsub_",
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
                        f" --name smafa_chunk_{chunk_idx:04d}"
                        f" --segregated-log-files --hours 48"
                        f" --command-file {cmd_file_path}"
                        f" --chunk-size 1 2>&1"
                    )
                    logging.info(f"Running: {mqsub_cmd}")
                    mqsub_stdout = extern.run(mqsub_cmd)
                    _mqwait(mqsub_stdout)
                finally:
                    os.unlink(cmd_file_path)

                logging.info(f"Round {chunk_idx + 1}/{max_rounds} mqsub jobs complete.")

            # Propagate state: stream clusters.tsv to update per-cluster stats.
            failed = []
            for marker_name, initial_chunk, _round_input, round_clusters, round_reps in round_jobs:
                if not os.path.exists(round_clusters) or os.path.getsize(round_clusters) == 0:
                    failed.append(marker_name)
                    continue

                old_rep_stats = marker_rep_stats[marker_name]

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
                            # Member was a previous round's rep: merge accumulated stats.
                            ms = old_rep_stats[member_seq]
                            s["n_seqs"]     += ms["n_seqs"]
                            s["n_query"]    += ms["n_query"]
                            s["max_hits"]    = max(s["max_hits"], ms["max_hits"])
                            s["total_hits"] += ms["total_hits"]
                        else:
                            # Member is a new sequence from this chunk.
                            unk, hits = chunk_seq_to_stats.get(member_seq, (1, 0))
                            s["n_seqs"]     += 1
                            s["n_query"]    += (1 if unk == 0 else 0)
                            s["max_hits"]    = max(s["max_hits"], hits)
                            s["total_hits"] += hits

                marker_rep_stats[marker_name] = new_rep_stats
                marker_current_reps_fasta[marker_name] = round_reps
                logging.info(
                    f"[{marker_name}] After round {chunk_idx + 1}: "
                    f"{len(new_rep_stats)} cluster(s)."
                )

            if failed:
                raise RuntimeError(
                    f"{len(failed)} smafa chunk job(s) failed in round {chunk_idx + 1}:\n"
                    + "\n".join(f"  {m}" for m in sorted(failed))
                )

        # ------------------------------------------------------------------
        # Phase C: write representatives.fasta per marker from accumulated stats.
        # ------------------------------------------------------------------
        for marker_name in list(pending.keys()):
            fasta_path = pending[marker_name]
            rep_stats  = marker_rep_stats[marker_name]
            reps_fasta = marker_current_reps_fasta[marker_name]

            marker_dir = os.path.join(abs_output_dir, marker_name)
            os.makedirs(marker_dir, exist_ok=True)

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
                            f"[{marker_name}] Rep sequence not found in reps.fasta — skipping: {rep_seq}"
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

            n_otus = 0
            if os.path.exists(fasta_path):
                with open(fasta_path) as fh:
                    n_otus = sum(1 for ln in fh if ln.startswith(">"))
            logging.info(
                f"[{marker_name}] Done — {n_clusters} cluster(s) from {n_otus} OTU(s)."
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
    """Wait for all job IDs parsed from mqsub stdout."""
    r = re.compile(r'^qsub stdout: (\d+\.aqua)$')
    job_ids = []
    for line in mqsub_log.split('\n'):
        m = r.match(line)
        if m:
            job_ids.append(m.group(1))
    logging.info(f"Waiting for {len(job_ids)} job(s) to finish ...")
    if not job_ids:
        return
    with tempfile.NamedTemporaryFile(mode="w", prefix="mqwait_", suffix=".ids") as f:
        f.write('\n'.join(job_ids) + '\n')
        f.flush()
        extern.run(f"mqwait -i {f.name}")


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
        logging.warning(f"Empty or missing chunk FASTA: {input_fasta} — writing empty output.")
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
    logging.info(f"Chunk done — {n_clusters} cluster(s) from {n_seqs} sequence(s).")

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
        logging.warning(f"[{marker_name}] Empty or missing FASTA — skipping.")
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
    logging.info(f"[{marker_name}] Done — {n_clusters} cluster(s) from {n_seqs} OTU(s).")
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

    # Build marker→domain map once (used by all paths)
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
        logging.warning("No metapackage provided — skipping off-target filtering.")

    if run_through_mqsub:
        # Phase 1: submit one extraction job per archive (skip already-done ones)
        pending = [p for p in archive_paths if not _already_extracted(output_dir, p)]
        skipped = len(archive_paths) - len(pending)
        if skipped:
            logging.info(f"Skipping {skipped} already-extracted archive(s).")
        submit_extraction_jobs(
            pending, marker_domains, markers_of_interest,
            output_dir, this_script_path
        )
        # Phase 2: cat per-sample FASTAs into per-marker FASTAs (local, fast)
        marker_to_fasta = cat_sample_fastas(archive_paths, output_dir, markers_of_interest, threads)
        if not marker_to_fasta:
            logging.warning("No markers with sequences — nothing to cluster.")
            return
        # Phase 3: submit one smafa cluster job per marker
        cluster_markers_via_mqsub(
            marker_to_fasta, output_dir, max_divergence, this_script_path,
            cluster_memory=cluster_memory,
            cluster_threads=cluster_threads,
            chunk_size=chunk_size,
        )
    else:
        # Local path: extract all samples in parallel, then cluster in parallel
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

        marker_to_fasta = cat_sample_fastas(archive_paths, output_dir, markers_of_interest, threads)
        if not marker_to_fasta:
            logging.warning("No markers with sequences — nothing to cluster.")
            return
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
        help="Number of sequences per smafa clustering chunk (default: 2000000)",
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
