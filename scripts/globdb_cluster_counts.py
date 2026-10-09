#!/usr/bin/env python3
"""Cluster one representative marker window per GlobDB species with smafa.

The output contains one row per marker with the number of input species and
the number of clusters at the requested divergence cutoff.
"""

import argparse
import collections
import logging
import os
import random
import sqlite3
import subprocess
import tempfile


DEFAULT_DB = ("/mnt/hpccs01/work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg/"
              "new_metapackage.sdb/otus.sqlite3")


def species_of(taxonomy):
    return taxonomy.rsplit(";", 1)[-1].strip()


def run_smafa(fasta, max_divergence, threads):
    cmd = [
        "cargo", "run", "--manifest-path", "/home/aroneys/src/smafa/Cargo.toml", "--",
        "cluster", "--max-divergence", str(max_divergence),
        "--threads", str(threads), "--input", fasta,
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    representatives = set()
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2 and fields[1]:
            representatives.add(fields[1])
    return len(representatives)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--output", default="results/globdb_marker_distances/GlobDB_r232_clusters.tsv")
    parser.add_argument("--marker-prefix", default="S3.")
    parser.add_argument("--max-marker-number", type=int, default=13)
    parser.add_argument("--max-divergence", type=int, default=2)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--subsample", type=int, default=None,
                        help="Randomly retain at most this many species per marker")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed used with --subsample (default: 42)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    markers = [(i, n) for i, n in conn.execute("SELECT id, marker FROM markers")
               if n.startswith(args.marker_prefix) and int(n.split(".")[1]) <= args.max_marker_number]
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    tmp = args.output + ".tmp"
    with open(tmp, "w") as out:
        out.write("marker\tspecies\tclusters\tmax_divergence\n")
        for marker_id, marker in sorted(markers, key=lambda x: x[1]):
            rows = conn.execute(
                "SELECT o.sequence, t.taxonomy FROM otus o JOIN taxonomy t "
                "ON o.taxonomy_id = t.id WHERE o.marker_id = ?", (marker_id,))
            sequences = collections.defaultdict(set)
            lengths = collections.Counter()
            for sequence, taxonomy in rows:
                sequence = sequence.replace("-", "N")
                sequences[species_of(taxonomy)].add(sequence)
                lengths[len(sequence)] += 1
            window = lengths.most_common(1)[0][0]
            representatives = {
                sp: sorted(seq for seq in seqs if len(seq) == window)[0]
                for sp, seqs in sequences.items()
                if any(len(seq) == window for seq in seqs)
            }
            if args.subsample is not None and len(representatives) > args.subsample:
                sampler = random.Random(args.seed)
                sampled_species = sampler.sample(sorted(representatives), args.subsample)
                representatives = {sp: representatives[sp] for sp in sampled_species}
            with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as fasta:
                for i, sequence in enumerate(representatives.values()):
                    fasta.write(f">species_{i}\n{sequence}\n")
                fasta_path = fasta.name
            try:
                logging.info("%s: clustering %d species", marker, len(representatives))
                clusters = run_smafa(fasta_path, args.max_divergence, args.threads)
            finally:
                os.unlink(fasta_path)
            out.write(f"{marker}\t{len(representatives)}\t{clusters}\t{args.max_divergence}\n")
            out.flush()
    os.replace(tmp, args.output)
    logging.info("Written to %s", args.output)


if __name__ == "__main__":
    main()
