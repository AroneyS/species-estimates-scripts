#!/usr/bin/env python3

###############################################################################
#
#    Summarise a species_estimation output folder (e.g.
#    results/species_estimation/20260704) into a tidy TSV of:
#        marker  stratum  stratum_value  sample_count  species_count
#    where species_count is the number of species (clusters) that appear in
#    exactly sample_count samples, for each marker and each metadata stratum:
#      - stratum="all", stratum_value="all": the unstratified count, from
#        {marker}/{marker}.tsv's n_samples column.
#      - stratum="year"/"host", stratum_value=<year>/<host or ecological>:
#        recomputed per stratum value by joining clusters.tsv and the
#        pre-dedup collated FASTA against sample metadata (year,
#        host_or_not) exported from the sandpiper duckdb.
#      - stratum="domain"/"phylum", stratum_value=<taxon>: assigned per
#        cluster representative by running representatives.fasta through
#        `singlem renew` against a SingleM metapackage (default: GlobDB_r232).
#
#    Strata don't sum to the "all" total: a species can appear in samples from
#    several years (or both host/ecological), so it's counted once per stratum
#    value it appears in, not once overall.
#
#    HEAVY STEPS RUN VIA MQSUB, ONE JOB PER MARKER (see run_species_estimation.py
#    for the same submit/wait pattern): both the year/host recomputation
#    (disk-based sort/join over tens-of-GB per-marker files) and the
#    domain/phylum query are too heavy to run on the login node, so this
#    script re-invokes itself as an mqsub worker via hidden --_stratify-marker
#    / --_taxonomy-marker flags, exactly as run_species_estimation.py does
#    for its own per-marker phases.
#
###############################################################################

import argparse
import csv
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter

MARKER_DIR_RE = re.compile(r"^S\d+\.")

STRATIFIED_COUNTS_FILENAME = "stratified_sample_counts.tsv"
TAXONOMY_COUNTS_FILENAME = "taxonomy_stratified_counts.tsv"

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)

DEFAULT_DUCKDB = "/scratch/microbiome/sandpiper2-globdb-renew/sandpiper_db_v2/sandpiper_39.duckdb"
DEFAULT_SAMPLE_METADATA = os.path.join(REPO_ROOT, "results", "sample_metadata",
                                        "sandpiper_39_sample_metadata.tsv")
DEFAULT_METAPACKAGE_DIR = "/work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg"
DEFAULT_SINGLEM_BIN = os.path.join(REPO_ROOT, ".pixi", "envs", "default", "bin", "singlem")
# singlem renew shells out to a bare `smafa` binary; this repo's own pixi env
# bundles one (.pixi/envs/default/bin/smafa) right alongside singlem itself,
# it just isn't on PATH unless this directory is added.
DEFAULT_SMAFA_BIN_DIR = os.path.join(REPO_ROOT, ".pixi", "envs", "default", "bin")

NA = "NA"


# ---------------------------------------------------------------------------
# Shared mqsub submit/wait helpers (same pattern as run_species_estimation.py)
# ---------------------------------------------------------------------------

def _submit_mqsub_batch(cmds, name, memory_gb, hours, threads=1, chunk_size=1):
    """
    Submit one mqsub job per command in cmds (via --command-file), returning
    the list of PBS job IDs immediately (does NOT wait) so multiple batches
    can be submitted before waiting on any of them.
    """
    if not cmds:
        return []
    with tempfile.NamedTemporaryFile(
        mode="w", prefix=f"{name}_mqsub_", suffix=".cmds", delete=False
    ) as cmd_file:
        cmd_file_path = cmd_file.name
        for cmd in cmds:
            cmd_file.write(cmd + "\n")
    try:
        mqsub_cmd = (
            f"mqsub -m {memory_gb} -t {threads} --name {name}"
            f" --segregated-log-files --hours {hours}"
            f" --command-file {cmd_file_path} --chunk-size {chunk_size}"
        )
        logging.info(f"Running: {mqsub_cmd}")
        # mqsub writes its log output -- including the "qsub stdout: X.aqua"
        # job-ID lines _parse_job_ids looks for -- to STDERR, not stdout, so
        # stderr must be merged into stdout (not just captured separately)
        # or job_ids silently comes back empty and _wait_for_jobs never
        # actually waits for anything.
        result = subprocess.run(
            mqsub_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"mqsub submission failed:\n{result.stdout}")
        return _parse_job_ids(result.stdout)
    finally:
        os.unlink(cmd_file_path)


def _parse_job_ids(mqsub_stdout):
    r = re.compile(r'^qsub stdout: (\d+\.aqua)$')
    job_ids = []
    for line in mqsub_stdout.split("\n"):
        m = r.match(line)
        if m:
            job_ids.append(m.group(1))
    return job_ids


def _wait_for_jobs(job_ids):
    if not job_ids:
        return
    logging.info(f"Waiting for {len(job_ids)} job(s) to finish ...")
    with tempfile.NamedTemporaryFile(mode="w", prefix="mqwait_", suffix=".ids") as f:
        f.write("\n".join(job_ids) + "\n")
        f.flush()
        subprocess.run(f"mqwait -i {f.name}", shell=True, check=True)


# ---------------------------------------------------------------------------
# Sample metadata export (from the sandpiper duckdb)
# ---------------------------------------------------------------------------

def export_sample_metadata(duckdb_path, output_path, host_column="host_or_not_mature"):
    """Export (acc, year, host_or_not) from the sandpiper duckdb to output_path."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    query = (
        "COPY (SELECT m.acc AS acc, a.collection_year AS year, "
        f"a.{host_column} AS host_or_not "
        "FROM ncbi_metadata m "
        "LEFT JOIN parsed_sample_attributes a ON a.run_id = m.id) "
        f"TO '{output_path}' (FORMAT CSV, DELIMITER '\\t', HEADER)"
    )
    cmd = ["duckdb", "-readonly", duckdb_path, "-c", query]
    logging.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"duckdb export failed:\n{result.stderr}")
    logging.info(f"Sample metadata written to {output_path}")


def load_sample_metadata(tsv_path):
    """Return {sample_acc: (year_str, host_str)}, missing values mapped to 'NA'."""
    metadata = {}
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        assert header == ["acc", "year", "host_or_not"], f"Unexpected header: {header}"
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            acc = parts[0]
            year = parts[1] if len(parts) > 1 and parts[1] else NA
            host = parts[2] if len(parts) > 2 and parts[2] else NA
            metadata[acc] = (year, host)
    return metadata


# ---------------------------------------------------------------------------
# Worker: year/host stratified sample counts (disk-based sort/join, O(1) mem)
#
#   clusters.tsv maps a deduplicated member sequence -> its cluster
#   representative sequence (both raw sequence text, not IDs).
#   collated_fastas/{marker}.fasta still has the pre-dedup, per-sample OTU
#   sequences (one entry per sample*OTU occurrence), with the sample name
#   embedded in the header (see otus_to_fasta()/_parse_sample_from_header in
#   run_species_estimation.py). Both files can be tens of GB at full scale, so
#   this follows the same O(1)-memory, disk-based sort/join pattern as
#   run_species_estimation.py's _compute_rep_sample_counts -- no per-cluster
#   membership dict is ever held in memory.
#
#   Pipeline (all via GNU sort/join subprocesses, LC_ALL=C for determinism):
#     1. sort clusters.tsv by member sequence (col 1).
#     2. stream collated_fastas/{marker}.fasta -> (sequence, sample, year,
#        host) tuples (year/host looked up from the small in-memory sample
#        metadata dict -- ~913k rows, trivial), then sort by sequence.
#     3. join (1) and (2) on sequence -> (rep_seq, sample, year, host) rows,
#        emitting two tagged tuples per row: (rep, "year", <year>, sample)
#        and (rep, "host", <host>, sample).
#     4. sort -u the tagged-tuple stream (dedupes a sample contributing
#        multiple divergent-but-clustered sequences, matching the "union not
#        sum" semantics of the existing n_samples field).
#     5. stream-count: consecutive (rep, stratum, value) runs -> distinct
#        sample count -> tally into a species-count histogram.
# ---------------------------------------------------------------------------

def _parse_sample_from_header(header):
    """
    Extract the sample name from a per-sample OTU header written by
    run_species_estimation.py's otus_to_fasta():
        >otu{i}|{sample}|{gene}|unknown={0or1}|hits={n}
    """
    parts = header.lstrip(">").split("|")
    return parts[1] if len(parts) > 1 else None


def _sorted_copy(src_path, dest_path, key_field=1, tmp_dir=None, threads=1):
    """LC_ALL=C sort src_path by the given tab-separated field, writing dest_path."""
    env = dict(os.environ, LC_ALL="C")
    cmd = ["sort", f"-k{key_field},{key_field}", "-t", "\t", "--buffer-size=1G",
           f"--parallel={threads}"]
    if tmp_dir:
        cmd += ["-T", tmp_dir]
    cmd += [src_path, "-o", dest_path]
    logging.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, env=env)


def _write_otu_tuples(collated_fasta_path, sample_metadata, out_path):
    """
    Stream collated_fastas/{marker}.fasta and write one line per OTU entry:
        sequence \t sample \t year \t host
    Sequence has '-' replaced with 'N' to match smafa's dash-to-N transform,
    the same way _write_rep_sample_pair does in run_species_estimation.py.
    """
    n_written = 0
    n_missing_sample = 0
    with open(collated_fasta_path) as in_fh, open(out_path, "w") as out_fh:
        cur_hdr = None
        seq_parts = []

        def _flush(hdr, parts):
            nonlocal n_written, n_missing_sample
            seq = "".join(parts)
            if hdr is None or not seq:
                return
            sample = _parse_sample_from_header(hdr)
            if sample is None:
                return
            year, host = sample_metadata.get(sample, (NA, NA))
            if sample not in sample_metadata:
                n_missing_sample += 1
            out_fh.write(f"{seq.replace('-', 'N')}\t{sample}\t{year}\t{host}\n")
            n_written += 1

        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush(cur_hdr, seq_parts)
                cur_hdr = line
                seq_parts = []
            else:
                seq_parts.append(line)
        _flush(cur_hdr, seq_parts)

    if n_missing_sample:
        logging.warning(
            f"{n_missing_sample}/{n_written} OTU entries had a sample accession "
            f"not found in the sample metadata table (treated as year=NA, host=NA)."
        )
    return n_written


def _join_and_tag(clusters_sorted_path, otu_sorted_path, tagged_path):
    """
    join clusters_sorted_path (member_seq \t rep_seq) with otu_sorted_path
    (seq \t sample \t year \t host) on sequence, then write two tagged rows
    per match to tagged_path:
        rep_seq \t year \t <year> \t sample
        rep_seq \t host \t <host> \t sample
    """
    env = dict(os.environ, LC_ALL="C")
    join_cmd = [
        "join", "-t", "\t", "-1", "1", "-2", "1",
        clusters_sorted_path, otu_sorted_path,
    ]
    logging.info(f"Running: {' '.join(join_cmd)}")
    n_rows = 0
    with subprocess.Popen(join_cmd, stdout=subprocess.PIPE, text=True, env=env) as proc, \
            open(tagged_path, "w") as out_fh:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            _seq, rep_seq, sample, year, host = parts[0], parts[1], parts[2], parts[3], parts[4]
            out_fh.write(f"{rep_seq}\tyear\t{year}\t{sample}\n")
            out_fh.write(f"{rep_seq}\thost\t{host}\t{sample}\n")
            n_rows += 1
        ret = proc.wait()
    if ret != 0:
        raise RuntimeError(f"join failed with exit code {ret}")
    return n_rows


def _dedup_sort(tagged_path, dedup_path, tmp_dir=None, threads=1):
    """sort -u tagged_path (rep, stratum, value, sample) -- collapses a sample
    contributing >1 divergent-but-clustered sequence to the same cluster."""
    env = dict(os.environ, LC_ALL="C")
    cmd = [
        "sort", "-u",
        "-k1,1", "-k2,2", "-k3,3", "-k4,4", "-t", "\t", "--buffer-size=1G",
        f"--parallel={threads}",
    ]
    if tmp_dir:
        cmd += ["-T", tmp_dir]
    cmd += [tagged_path, "-o", dedup_path]
    logging.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, env=env)


def _tally_histogram(dedup_sorted_path):
    """
    Stream the sorted-unique (rep, stratum, value, sample) file and, for each
    consecutive (rep, stratum, value) run, count distinct samples, then tally
    into counts[(stratum, value, sample_count)] += 1.
    """
    counts = Counter()
    cur_key = None
    cur_count = 0

    def _flush():
        if cur_key is not None and cur_count > 0:
            stratum, value = cur_key[1], cur_key[2]
            counts[(stratum, value, cur_count)] += 1

    with open(dedup_sorted_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            rep, stratum, value = parts[0], parts[1], parts[2]
            key = (rep, stratum, value)
            if key == cur_key:
                cur_count += 1
            else:
                _flush()
                cur_key = key
                cur_count = 1
        _flush()

    return counts


def write_histogram(counts, output_path):
    with open(output_path, "w") as fh:
        fh.write("stratum\tstratum_value\tsample_count\tspecies_count\n")
        for (stratum, value, sample_count), species_count in sorted(counts.items()):
            if sample_count > 0 and species_count > 0:
                fh.write(f"{stratum}\t{value}\t{sample_count}\t{species_count}\n")


def worker_stratify_marker(output_dir, marker_name, sample_metadata_path, tmp_dir=None, threads=1):
    """mqsub worker: recompute year/host stratified sample counts for one marker."""
    marker_dir = os.path.join(output_dir, marker_name)
    collated_fasta_path = os.path.join(output_dir, "collated_fastas", f"{marker_name}.fasta")
    clusters_path = os.path.join(marker_dir, "clusters.tsv")
    output_path = os.path.join(marker_dir, STRATIFIED_COUNTS_FILENAME)
    tmp_dir = tmp_dir or marker_dir

    logging.info(f"Loading sample metadata from {sample_metadata_path} ...")
    sample_metadata = load_sample_metadata(sample_metadata_path)
    logging.info(f"Loaded metadata for {len(sample_metadata)} sample(s).")

    clusters_sorted = os.path.join(tmp_dir, f".{marker_name}.clusters_by_seq.tsv")
    otu_tuples = os.path.join(tmp_dir, f".{marker_name}.otu_tuples.tsv")
    otu_sorted = os.path.join(tmp_dir, f".{marker_name}.otu_tuples_sorted.tsv")
    tagged = os.path.join(tmp_dir, f".{marker_name}.tagged.tsv")
    dedup_sorted = os.path.join(tmp_dir, f".{marker_name}.tagged_dedup.tsv")

    try:
        logging.info(f"[{marker_name}] Sorting clusters.tsv by member sequence ...")
        _sorted_copy(clusters_path, clusters_sorted, key_field=1, tmp_dir=tmp_dir, threads=threads)

        logging.info(f"[{marker_name}] Building (sequence, sample, year, host) tuples ...")
        n_otus = _write_otu_tuples(collated_fasta_path, sample_metadata, otu_tuples)
        logging.info(f"[{marker_name}] Wrote {n_otus} OTU tuple(s). Sorting by sequence ...")
        _sorted_copy(otu_tuples, otu_sorted, key_field=1, tmp_dir=tmp_dir, threads=threads)

        logging.info(f"[{marker_name}] Joining clusters <-> OTU tuples on sequence ...")
        n_joined = _join_and_tag(clusters_sorted, otu_sorted, tagged)
        if n_joined < n_otus:
            logging.warning(
                f"[{marker_name}] {n_otus - n_joined}/{n_otus} OTU sequence(s) had no "
                f"matching entry in clusters.tsv (e.g. sequences smafa drops -- see "
                f"_process_cluster_file's docstring in run_species_estimation.py) -- "
                f"excluded from stratified counts."
            )
        logging.info(f"[{marker_name}] {n_joined} joined row(s). Deduplicating ...")
        _dedup_sort(tagged, dedup_sorted, tmp_dir=tmp_dir, threads=threads)

        logging.info(f"[{marker_name}] Tallying stratified sample-count histogram ...")
        counts = _tally_histogram(dedup_sorted)
    finally:
        for f in (clusters_sorted, otu_tuples, otu_sorted, tagged, dedup_sorted):
            if os.path.exists(f):
                os.unlink(f)

    write_histogram(counts, output_path)
    logging.info(f"[{marker_name}] Stratified sample counts written to {output_path}")


# ---------------------------------------------------------------------------
# Worker: domain/phylum stratified counts via `singlem renew`
#
#   Domain/phylum, unlike year/host, is a property of the cluster
#   REPRESENTATIVE itself (one value per cluster) -- no per-sample join is
#   needed. This runs representatives.fasta through `singlem renew` against a
#   SingleM metapackage to get a taxonomy string per representative, extracts
#   the d__/p__ tokens, and joins directly against {marker}.tsv's n_samples
#   column by representative id (see _representatives_id_from_header).
# ---------------------------------------------------------------------------

def _representatives_id_from_header(header):
    """
    Extract the same 'id' that _representatives_fasta_to_tsv (run_species_estimation.py)
    writes for this header: the '|'-joined leading tokens WITHOUT an '=' sign.

    This must match exactly, since it's later used to join singlem renew's
    output back to {marker}.tsv's 'id' column. Two header formats exist:
      mqsub path:  >seq123|unknown=0|hits=5|...                  -> id = "seq123"
      local path:  >otu0|sample1|gene1|unknown=0|hits=5|...      -> id = "otu0|sample1|gene1"
    Taking just the first '|'-token (as if only the mqsub format existed)
    would silently produce the wrong id -- and hence an all-NA taxonomy join
    -- for output produced via the local (non-mqsub) clustering path.
    """
    tokens = header.lstrip(">").split("|")
    return "|".join(tok for tok in tokens if "=" not in tok)


ARCHIVE_OTU_TABLE_FIELDS = [
    "gene", "sample", "sequence", "num_hits", "coverage", "taxonomy",
    "read_names", "nucleotides_aligned", "taxonomy_by_known?",
    "read_unaligned_sequences", "equal_best_hit_taxonomies",
    "taxonomy_assignment_method",
]


def _iter_representative_batches(rep_fasta_path, batch_size):
    """
    Stream representatives.fasta and yield lists of (rep_id, sequence) tuples,
    each up to batch_size long (the last batch may be shorter).

    Chunking representatives.fasta this way -- rather than converting the
    whole file to one archive OTU table -- is what keeps memory bounded: see
    _write_archive_otu_table_batch and the OOM note on worker_taxonomy_marker.
    """
    batch = []
    cur_hdr = None
    seq_parts = []

    def _flush(hdr, parts):
        seq = "".join(parts)
        if hdr is None or not seq:
            return None
        rep_id = _representatives_id_from_header(hdr)
        return (rep_id, seq)

    with open(rep_fasta_path) as in_fh:
        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                entry = _flush(cur_hdr, seq_parts)
                if entry is not None:
                    batch.append(entry)
                    if len(batch) >= batch_size:
                        yield batch
                        batch = []
                cur_hdr = line
                seq_parts = []
            else:
                seq_parts.append(line)
        entry = _flush(cur_hdr, seq_parts)
        if entry is not None:
            batch.append(entry)
    if batch:
        yield batch


def _write_archive_otu_table_batch(batch, marker_name, archive_path):
    """
    Write one batch of (rep_id, sequence) tuples as a minimal SingleM archive
    OTU table (JSON, version 4), suitable as
    `singlem renew --input-archive-otu-table` input.

    `renew` (singlem/renew.py) only needs read_names / read_unaligned_sequences
    / nucleotides_aligned present and length-matched -- it never checks that
    the "unaligned" sequence is actually a longer raw read, so the
    representative's own (already-windowed) sequence is reused for both.
    "sample" is set to the representative's id (matching {marker}.tsv's 'id'
    column, see _representatives_id_from_header) so the renewed output can be
    joined back the same way. alignment_hmm_sha256s/singlem_package_sha256s
    are set to the placeholder "na" -- ArchiveOtuTable.read() never validates
    them (same precedent as singlem/condense.py).
    """
    with open(archive_path, "w") as out_fh:
        out_fh.write('{"version": 4, "alignment_hmm_sha256s": "na", '
                     '"singlem_package_sha256s": "na", "fields": ')
        json.dump(ARCHIVE_OTU_TABLE_FIELDS, out_fh)
        out_fh.write(', "otus": [')
        for i, (rep_id, seq) in enumerate(batch):
            row = [marker_name, rep_id, seq, 1, 1.0, "", [rep_id], [len(seq)],
                   False, [seq], [], ""]
            if i:
                out_fh.write(",")
            json.dump(row, out_fh)
        out_fh.write("]}")


def _run_singlem_renew(archive_path, metapackage_dir, results_path, singlem_bin,
                        smafa_bin_dir, threads):
    """
    Run `singlem renew` against metapackage_dir, writing a renewed archive
    OTU table to results_path.

    The default --assignment-method (smafa_naive_then_diamond) tries the fast
    smafa nearest-neighbour lookup first (same mechanism as singlem query),
    then falls back to DIAMOND blastx for anything unmatched -- unlike
    query's strict max-divergence cutoff (which leaves most representatives
    with no hit at all), this gives essentially every representative a
    taxonomy call. Validated live: 10/10 on a small test, vs. partial
    coverage from query on the same input.

    Crucially, renew operates directly on the already-extracted/aligned OTU
    sequence -- it does NOT re-run graftM's HMM-search/ORF-calling steps (see
    singlem/renew.py), so it doesn't hit the process-substitution/`smafa`
    issues that `singlem pipe` does on this representative-sequence input.

    smafa_bin_dir is prepended to PATH because singlem shells out to a bare
    `smafa` binary; this repo's own pixi env already bundles one
    (.pixi/envs/default/bin/smafa) alongside singlem itself.
    """
    env = dict(os.environ)
    env["PATH"] = f"{smafa_bin_dir}:{env.get('PATH', '')}"
    cmd = [
        singlem_bin, "renew",
        "--input-archive-otu-table", archive_path,
        "--metapackage", metapackage_dir,
        "--archive-otu-table", results_path,
        "--threads", str(threads),
    ]
    logging.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"singlem renew failed:\n{result.stderr}")


def _parse_taxonomy_token(taxonomy, prefix):
    """Extract the token starting with prefix (e.g. 'd__', 'p__') from a
    ';'-separated taxonomy string (renew's output is prefixed with a "Root;"
    token, which simply won't match either prefix), or NA if absent."""
    for token in taxonomy.split(";"):
        token = token.strip()
        if token.startswith(prefix):
            return token
    return NA


def _load_rep_taxonomy(results_path):
    """
    Return {rep_id: (domain, phylum)} from a renewed archive OTU table.

    Loads the whole JSON document at once (not streamed) -- fine as long as
    results_path holds one CHUNK's worth of representatives (see
    worker_taxonomy_marker), not an entire ~20M-representative marker: a
    real full-marker attempt OOM-killed singlem renew's own
    ArchiveOtuTable.read() (also a whole-document json.load, outside this
    script's control) at a 16GB job memory limit, before renew even started
    assigning taxonomy.
    """
    with open(results_path) as fh:
        data = json.load(fh)
    fields = data["fields"]
    sample_idx = fields.index("sample")
    taxonomy_idx = fields.index("taxonomy")

    rep_taxonomy = {}
    for row in data["otus"]:
        rep_id = row[sample_idx]
        taxonomy = row[taxonomy_idx] or ""
        domain = _parse_taxonomy_token(taxonomy, "d__")
        phylum = _parse_taxonomy_token(taxonomy, "p__")
        rep_taxonomy[rep_id] = (domain, phylum)
    return rep_taxonomy


def _load_marker_n_samples(marker_tsv_path):
    """Return {id: n_samples} for every cluster representative in {marker}.tsv."""
    n_samples_by_id = {}
    with open(marker_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is not None and "n_samples" not in reader.fieldnames:
            raise ValueError(f"{marker_tsv_path} has no 'n_samples' column -- cannot summarise.")
        for row in reader:
            n_samples_by_id[row["id"]] = int(row["n_samples"])
    return n_samples_by_id


DEFAULT_TAXONOMY_CHUNK_SIZE = 2_000_000


def worker_taxonomy_marker(output_dir, marker_name, metapackage_dir, singlem_bin,
                            smafa_bin_dir, threads=1, tmp_dir=None,
                            chunk_size=DEFAULT_TAXONOMY_CHUNK_SIZE):
    """
    mqsub worker: assign domain/phylum per cluster representative via
    singlem renew, then tally a stratified sample-count histogram.

    Processes representatives.fasta in chunks of chunk_size, each as its own
    renew invocation, rather than one archive OTU table for the whole marker.
    This is not an optimisation -- a whole-marker attempt on a ~20M-
    representative marker OOM-killed inside singlem renew's own
    ArchiveOtuTable.read() (a whole-document json.load, before renew even
    started assigning taxonomy) at a 16GB job memory limit. Chunking keeps
    per-invocation memory bounded regardless of marker size; the
    {marker}.tsv id -> n_samples lookup is still loaded once up front (one
    dict of scalars, not the much larger archive-OTU-table structure).
    """
    marker_dir = os.path.join(output_dir, marker_name)
    rep_fasta_path = os.path.join(marker_dir, "representatives.fasta")
    marker_tsv_path = os.path.join(marker_dir, f"{marker_name}.tsv")
    output_path = os.path.join(marker_dir, TAXONOMY_COUNTS_FILENAME)
    tmp_dir = tmp_dir or marker_dir

    logging.info(f"[{marker_name}] Loading n_samples per representative ...")
    n_samples_by_id = _load_marker_n_samples(marker_tsv_path)
    logging.info(f"[{marker_name}] Loaded {len(n_samples_by_id)} representative(s).")

    archive_path = os.path.join(tmp_dir, f".{marker_name}.rep_archive_otu_table.json")
    results_path = os.path.join(tmp_dir, f".{marker_name}.renew_results.json")

    counts = Counter()
    n_total = 0
    n_with_taxonomy = 0
    try:
        for batch_idx, batch in enumerate(_iter_representative_batches(rep_fasta_path, chunk_size)):
            logging.info(f"[{marker_name}] Chunk {batch_idx + 1}: {len(batch)} representative(s) ...")
            _write_archive_otu_table_batch(batch, marker_name, archive_path)
            _run_singlem_renew(archive_path, metapackage_dir, results_path,
                                singlem_bin, smafa_bin_dir, threads)
            rep_taxonomy = _load_rep_taxonomy(results_path)

            for rep_id, _seq in batch:
                n_total += 1
                domain, phylum = rep_taxonomy.get(rep_id, (NA, NA))
                if rep_id in rep_taxonomy:
                    n_with_taxonomy += 1
                n_samples = n_samples_by_id.get(rep_id)
                if n_samples is None:
                    continue
                counts[("domain", domain, n_samples)] += 1
                counts[("phylum", phylum, n_samples)] += 1
    finally:
        for f in (archive_path, results_path):
            if os.path.exists(f):
                os.unlink(f)

    logging.info(f"[{marker_name}] {n_with_taxonomy}/{n_total} representative(s) got a taxonomy call.")
    write_histogram(counts, output_path)
    logging.info(f"[{marker_name}] Taxonomy stratified counts written to {output_path}")


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def find_marker_dirs(input_dir, restrict=None):
    """Return a sorted list of (marker_name, marker_dir) for every marker subdirectory."""
    markers = []
    for name in sorted(os.listdir(input_dir)):
        marker_dir = os.path.join(input_dir, name)
        if not os.path.isdir(marker_dir) or not MARKER_DIR_RE.match(name):
            continue
        if restrict and name not in restrict:
            continue
        markers.append((name, marker_dir))
    return markers


def count_all_sample_counts(marker_tsv_path):
    """
    Stream a marker's {marker}.tsv and return a Counter mapping
    n_samples -> number of species (rows) with that n_samples value.
    """
    counts = Counter()
    with open(marker_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is not None and "n_samples" not in reader.fieldnames:
            raise ValueError(f"{marker_tsv_path} has no 'n_samples' column -- cannot summarise.")
        for row in reader:
            counts[int(row["n_samples"])] += 1
    return counts


def read_histogram(tsv_path):
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            yield (
                row["stratum"],
                row["stratum_value"],
                int(row["sample_count"]),
                int(row["species_count"]),
            )


def summarise(input_dir, markers):
    """Yield (marker, stratum, stratum_value, sample_count, species_count) rows for all markers."""
    for marker_name, marker_dir in markers:
        marker_tsv_path = os.path.join(marker_dir, f"{marker_name}.tsv")
        if not os.path.isfile(marker_tsv_path):
            logging.warning(f"No {marker_name}.tsv found in {marker_dir} -- skipping.")
            continue

        logging.info(f"Aggregating {marker_name} ...")
        counts = count_all_sample_counts(marker_tsv_path)
        for sample_count in sorted(counts):
            species_count = counts[sample_count]
            if sample_count > 0 and species_count > 0:
                yield marker_name, "all", "all", sample_count, species_count

        for filename in (STRATIFIED_COUNTS_FILENAME, TAXONOMY_COUNTS_FILENAME):
            path = os.path.join(marker_dir, filename)
            if os.path.isfile(path):
                for stratum, value, sample_count, species_count in read_histogram(path):
                    if sample_count > 0 and species_count > 0:
                        yield marker_name, stratum, value, sample_count, species_count
            else:
                logging.info(f"[{marker_name}] No {filename} found -- those stratum rows will be absent.")


def run_orchestrator(args):
    markers = find_marker_dirs(args.input_dir, restrict=set(args.markers) if args.markers else None)
    if not markers:
        raise SystemExit(f"No marker subdirectories found in {args.input_dir}")

    stratify_job_ids = []
    taxonomy_job_ids = []

    if not args.skip_stratify_metadata:
        if not os.path.isfile(args.sample_metadata):
            logging.info(f"{args.sample_metadata} not found -- exporting from {args.duckdb} ...")
            export_sample_metadata(args.duckdb, args.sample_metadata, args.host_column)

        pending = [
            (name, d) for name, d in markers
            if not os.path.isfile(os.path.join(d, STRATIFIED_COUNTS_FILENAME))
        ]
        if pending:
            logging.info(f"Submitting {len(pending)} year/host stratification job(s) via mqsub ...")
            cmds = [
                f"python3 {THIS_SCRIPT_PATH} {args.input_dir}"
                f" --_stratify-marker {marker_name}"
                f" --sample-metadata {args.sample_metadata}"
                f" --threads {args.sort_threads}"
                for marker_name, _ in pending
            ]
            stratify_job_ids = _submit_mqsub_batch(
                cmds, "stratify_metadata", args.stratify_memory, args.stratify_hours,
                threads=args.sort_threads,
            )
        else:
            logging.info("All markers already have year/host stratified counts.")

    if not args.skip_taxonomy:
        pending = [
            (name, d) for name, d in markers
            if not os.path.isfile(os.path.join(d, TAXONOMY_COUNTS_FILENAME))
        ]
        if pending:
            logging.info(f"Submitting {len(pending)} domain/phylum taxonomy job(s) via mqsub ...")
            cmds = [
                f"python3 {THIS_SCRIPT_PATH} {args.input_dir}"
                f" --_taxonomy-marker {marker_name}"
                f" --metapackage {args.metapackage}"
                f" --singlem-bin {args.singlem_bin}"
                f" --smafa-bin-dir {args.smafa_bin_dir}"
                f" --threads {args.taxonomy_threads}"
                f" --taxonomy-chunk-size {args.taxonomy_chunk_size}"
                for marker_name, _ in pending
            ]
            taxonomy_job_ids = _submit_mqsub_batch(
                cmds, "taxonomy_renew", args.taxonomy_memory, args.taxonomy_hours,
                threads=args.taxonomy_threads,
            )
        else:
            logging.info("All markers already have domain/phylum stratified counts.")

    _wait_for_jobs(stratify_job_ids + taxonomy_job_ids)

    output_path = args.output or os.path.join(args.input_dir, "sample_count_summary.tsv")
    with open(output_path, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["marker", "stratum", "stratum_value", "sample_count", "species_count"])
        for row in summarise(args.input_dir, markers):
            writer.writerow(row)
    logging.info(f"Summary written to {output_path}")


THIS_SCRIPT_PATH = os.path.abspath(__file__)


def main():
    parser = argparse.ArgumentParser(
        description="Summarise a species_estimation output folder into a tidy TSV of "
                     "marker, stratum, stratum_value, sample_count, species_count -- "
                     "submitting mqsub sub-jobs per marker for the heavy year/host and "
                     "domain/phylum stratification steps."
    )
    parser.add_argument(
        "input_dir",
        help="Path to a species_estimation output folder "
             "(e.g. results/species_estimation/20260704)",
    )
    parser.add_argument("-o", "--output", default=None,
                         help="Output TSV path (default: <input_dir>/sample_count_summary.tsv)")
    parser.add_argument("--markers", nargs="+", metavar="MARKER",
                         help="Restrict processing to these marker names (default: all markers)")

    # --- Year/host stratification ---
    parser.add_argument("--skip-stratify-metadata", action="store_true",
                         help="Skip year/host stratification (only 'all' rows are written).")
    parser.add_argument("--duckdb", default=DEFAULT_DUCKDB,
                         help="Sandpiper duckdb path, used only if --sample-metadata doesn't "
                              "already exist.")
    parser.add_argument("--sample-metadata", default=DEFAULT_SAMPLE_METADATA,
                         help="Sample metadata TSV (acc, year, host_or_not); exported from "
                              "--duckdb if missing.")
    parser.add_argument("--host-column", default="host_or_not_mature",
                         help="parsed_sample_attributes column for host_or_not "
                              "(default: host_or_not_mature).")
    parser.add_argument("--stratify-memory", type=int, default=16, metavar="GB",
                         help="Memory (GB) per year/host stratification mqsub job (default: 16).")
    parser.add_argument("--stratify-hours", type=int, default=48, metavar="H",
                         help="Walltime (hours) per year/host stratification mqsub job (default: 48).")
    parser.add_argument("--sort-threads", type=int, default=4, metavar="N",
                         help="Threads for GNU sort's --parallel in the stratification job (default: 4).")

    # --- Domain/phylum taxonomy ---
    parser.add_argument("--skip-taxonomy", action="store_true",
                         help="Skip domain/phylum assignment via singlem renew.")
    parser.add_argument("--metapackage", default=DEFAULT_METAPACKAGE_DIR,
                         help="SingleM metapackage directory (default: GlobDB_r232). Used with "
                              "singlem renew's default --assignment-method "
                              "(smafa_naive_then_diamond): a fast nearest-neighbour lookup first, "
                              "falling back to DIAMOND blastx for anything unmatched -- unlike "
                              "singlem query's strict divergence cutoff, this assigns essentially "
                              "every representative a taxonomy (validated 10/10 on a live test).")
    parser.add_argument("--singlem-bin", default=DEFAULT_SINGLEM_BIN,
                         help="Path to the singlem executable.")
    parser.add_argument("--smafa-bin-dir", default=DEFAULT_SMAFA_BIN_DIR,
                         help="Directory containing a smafa binary, prepended to PATH for "
                              "singlem renew (which shells out to `smafa`).")
    parser.add_argument("--taxonomy-threads", type=int, default=1, metavar="N",
                         help="Threads for singlem renew (default: 1; unlike singlem query, "
                              "renew's --threads is a real int and safe to raise).")
    parser.add_argument("--taxonomy-chunk-size", type=int, default=DEFAULT_TAXONOMY_CHUNK_SIZE,
                         metavar="N",
                         help="Representatives per singlem renew invocation (default: "
                              f"{DEFAULT_TAXONOMY_CHUNK_SIZE}). REQUIRED, not just a tuning knob: "
                              "a whole-marker attempt (~20M representatives for the largest "
                              "markers) OOM-killed inside singlem renew's own "
                              "ArchiveOtuTable.read() (a whole-document json.load this script "
                              "doesn't control) at a 16GB job memory limit, before renew even "
                              "started assigning taxonomy. Sized to bound the largest marker "
                              "(~20M reps) to ~10 renew invocations: a live test saw one "
                              "200-representative chunk take ~112 minutes under cluster "
                              "contention (vs. ~3.5 min for an identical chunk right after), so "
                              "with --taxonomy-hours capped at this cluster's 48h maximum, too "
                              "many small chunks risks the total exceeding the walltime even "
                              "though no single chunk is slow on its own. Lower this only if "
                              "chunks OOM at your --taxonomy-memory; each reduction multiplies "
                              "chunk count and therefore worst-case total time.")
    parser.add_argument("--taxonomy-memory", type=int, default=32, metavar="GB",
                         help="Memory (GB) per taxonomy mqsub job (default: 32). A 500K-"
                              "representative chunk used under 1GB in testing, so 32GB gives "
                              "headroom at the larger default chunk size above -- but the "
                              "500K-scale test is the only data point; not confirmed at 2M.")
    parser.add_argument("--taxonomy-hours", type=int, default=48, metavar="H",
                         help="Walltime (hours) per taxonomy mqsub job (default: 48 -- this "
                              "cluster's maximum, not a choice). Per-renew-invocation time has "
                              "been highly inconsistent in live testing under cluster contention "
                              "(69s to ~112min for comparable/identical-size inputs), so the "
                              "worst case for ~10 chunks (see --taxonomy-chunk-size) could still "
                              "approach this ceiling; there is no larger value available if it "
                              "doesn't fit.")

    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--quiet", action="store_true")

    # --- Internal: used only when this script is re-invoked as an mqsub worker ---
    parser.add_argument("--_stratify-marker", metavar="NAME", help=argparse.SUPPRESS)
    parser.add_argument("--_taxonomy-marker", metavar="NAME", help=argparse.SUPPRESS)
    parser.add_argument("--threads", type=int, default=1, help=argparse.SUPPRESS)

    args = parser.parse_args()

    level = logging.DEBUG if args.debug else (logging.ERROR if args.quiet else logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if args._stratify_marker:
        worker_stratify_marker(
            args.input_dir, args._stratify_marker, args.sample_metadata, threads=args.threads
        )
        sys.exit(0)

    if args._taxonomy_marker:
        worker_taxonomy_marker(
            args.input_dir, args._taxonomy_marker, args.metapackage, args.singlem_bin,
            args.smafa_bin_dir, threads=args.threads, chunk_size=args.taxonomy_chunk_size,
        )
        sys.exit(0)

    run_orchestrator(args)


if __name__ == "__main__":
    main()
