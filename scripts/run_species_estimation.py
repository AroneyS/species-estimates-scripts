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


def _run_smafa_cluster_file(fasta_path, max_divergence, marker_name):
    """
    Run smafa cluster directly on fasta_path, writing output to a sibling temp file.
    Returns the path to the cluster TSV file (caller is responsible for deletion).
    """
    cluster_tmp = fasta_path + ".smafa_clusters.tsv"
    cmd = [
        "smafa", "cluster",
        "--max-divergence", str(max_divergence),
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
    rep_seqs = set()
    with open(cluster_path) as cf:
        for line in cf:
            line = line.rstrip("\n")
            if line:
                parts = line.split("\t")
                if parts:
                    rep_seqs.add(parts[0].lstrip(">"))

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
            rep_seq = parts[0].lstrip(">")
            member_seq = parts[1].lstrip(">")
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
    FASTAs into one file per marker.  Also performs the global sort
    (GlobDB-known first, then num_hits desc) by re-reading sequence headers.

    Returns a dict: marker_name -> collated_fasta_path.
    """
    fasta_dir = os.path.join(output_dir, FASTA_SUBDIR)
    os.makedirs(fasta_dir, exist_ok=True)

    # Discover which markers exist across all samples
    sample_root = os.path.join(output_dir, SAMPLE_FASTA_SUBDIR)
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

    logging.info(
        f"Concatenating sample FASTAs for {len(marker_to_sample_fastas)} marker(s) "
        f"across {threads} thread(s) ..."
    )
    worker_args = [
        (marker_name, sample_fastas, os.path.join(fasta_dir, f"{marker_name}.fasta"))
        for marker_name, sample_fastas in marker_to_sample_fastas.items()
    ]
    marker_to_fasta = {}
    if threads > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for marker_name, collated_path in executor.map(_cat_one_marker, worker_args):
                marker_to_fasta[marker_name] = collated_path
    else:
        for marker_name, collated_path in map(_cat_one_marker, worker_args):
            marker_to_fasta[marker_name] = collated_path

    logging.info(f"Collation complete — {len(marker_to_fasta)} marker(s) ready for clustering.")
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

def cluster_markers_locally(marker_to_fasta, output_dir, max_divergence, threads):
    """Run smafa cluster for each marker in parallel (file-based, no fasta_str)."""

    def _cluster_one(args_tuple):
        (marker_name, fasta_path, max_divergence) = args_tuple
        num_seqs = sum(1 for ln in open(fasta_path) if ln.startswith(">"))
        logging.info(f"[{marker_name}] Clustering {num_seqs} OTU(s) ...")
        cluster_tmp = _run_smafa_cluster_file(fasta_path, max_divergence, marker_name)
        try:
            cluster_tsv, rep_fasta = _process_cluster_file(fasta_path, cluster_tmp, marker_name)
        finally:
            if os.path.exists(cluster_tmp):
                os.unlink(cluster_tmp)
        n_clusters = len(set(
            ln.split("\t")[0] for ln in cluster_tsv.splitlines() if ln.strip()
        ))
        logging.info(f"[{marker_name}] Done — {n_clusters} cluster(s) from {num_seqs} OTU(s).")
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
# Phase 3 (mqsub): submit one smafa cluster job per marker
# ---------------------------------------------------------------------------

def _already_clustered(output_dir, marker_name):
    """Return True if this marker's clusters.tsv already exists and is non-empty."""
    clusters_path = os.path.join(output_dir, marker_name, "clusters.tsv")
    return os.path.exists(clusters_path) and os.path.getsize(clusters_path) > 0


def cluster_markers_via_mqsub(
    marker_to_fasta, output_dir, max_divergence, this_script_path, cluster_memory=64
):
    """Submit one mqsub job per marker for smafa cluster."""
    abs_output_dir = os.path.abspath(output_dir)

    # Skip markers already clustered successfully
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
        logging.info(
            f"Submitting {len(pending)} smafa cluster job(s) via mqsub ..."
        )
        with tempfile.NamedTemporaryFile(
            mode="w", prefix="smafa_cluster_mqsub_", suffix=".cmds", delete=False
        ) as cmd_file:
            cmd_file_path = cmd_file.name
            for marker_name, fasta_path in pending.items():
                cmd = (
                    f"python3 {this_script_path}"
                    f" --_cluster-marker-fasta {os.path.abspath(fasta_path)}"
                    f" --_marker-name {marker_name}"
                    f" --output-directory {abs_output_dir}"
                    f" --max-divergence {max_divergence}"
                )
                cmd_file.write(cmd + "\n")

        try:
            mqsub_cmd = (
                f"mqsub -m {cluster_memory} --name smafa_cluster"
                f" --segregated-log-files"
                f" --hours 12"
                f" --command-file {cmd_file_path}"
                f" --chunk-size 1 2>&1"
            )
            logging.info(f"Running: {mqsub_cmd}")
            mqsub_stdout = extern.run(mqsub_cmd)
            logging.info(f"Submitted {len(pending)} mqsub job(s).")
            _mqwait(mqsub_stdout)
        finally:
            os.unlink(cmd_file_path)

        logging.info("All smafa cluster jobs finished.")

        # Check for failures: any pending marker missing its clusters.tsv
        failed = [
            m for m in pending
            if not _already_clustered(abs_output_dir, m)
        ]
        if failed:
            raise RuntimeError(
                f"{len(failed)} smafa cluster job(s) failed (no output produced):\n"
                + "\n".join(f"  {m}" for m in sorted(failed))
            )

    # Aggregate summary across all marker dirs (each worker only wrote its own
    # cluster files; writing summary here ensures all markers are included).
    summary_rows = []
    for marker_name, fasta_path in marker_to_fasta.items():
        marker_dir = os.path.join(abs_output_dir, marker_name)
        clusters_path = os.path.join(marker_dir, "clusters.tsv")
        num_clusters = 0
        if os.path.exists(clusters_path):
            with open(clusters_path) as fh:
                num_clusters = sum(1 for ln in fh if ln.strip())
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


def worker_cluster_marker(fasta_path, marker_name, output_dir, max_divergence):
    """Entry point for per-marker smafa cluster mqsub jobs (file-based, no fasta_str)."""
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        logging.warning(f"[{marker_name}] Empty or missing FASTA — skipping.")
        return
    num_seqs = sum(1 for ln in open(fasta_path) if ln.startswith(">"))
    logging.info(f"[{marker_name}] Clustering {num_seqs} OTU(s) ...")
    cluster_tmp = _run_smafa_cluster_file(fasta_path, max_divergence, marker_name)
    try:
        cluster_tsv, rep_fasta = _process_cluster_file(fasta_path, cluster_tmp, marker_name)
    finally:
        if os.path.exists(cluster_tmp):
            os.unlink(cluster_tmp)
    result = (marker_name, cluster_tsv, rep_fasta)
    _write_cluster_results([result], {marker_name: fasta_path}, output_dir)
    n_clusters = len(set(
        ln.split("\t")[0] for ln in cluster_tsv.splitlines() if ln.strip()
    ))
    logging.info(f"[{marker_name}] Done — {n_clusters} cluster(s) from {num_seqs} OTU(s).")
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
            marker_to_fasta, output_dir, max_divergence, threads
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
        required=True,
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

    # --- Internal: used only when this script is re-invoked by an mqsub worker ---
    parser.add_argument("--_cluster-marker-fasta", metavar="FILE", help=argparse.SUPPRESS)
    parser.add_argument("--_marker-name", metavar="NAME", help=argparse.SUPPRESS)
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

    # --- Worker path: per-marker smafa cluster (invoked by mqsub) ---
    if args._cluster_marker_fasta:
        if not args._marker_name:
            parser.error("--_cluster-marker-fasta requires --_marker-name")
        worker_cluster_marker(
            fasta_path=args._cluster_marker_fasta,
            marker_name=args._marker_name,
            output_dir=args.output_directory,
            max_divergence=args.max_divergence,
        )
        sys.exit(0)

    # --- Normal path ---
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
    )
