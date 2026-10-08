#!/usr/bin/env python3

###############################################################################
#
#    Within- and between-species distances of SingleM marker windows in a
#    SingleM reference database (default: GlobDB_r232), to assess the
#    clustering divergence cutoff.
#
#    For each marker:
#      - intra_pair: Hamming distances between distinct window sequences of
#        genomes in the same species (up to MAX_INTRA_SEQS per species).
#      - intra_diameter: per species, the largest such distance (0 when all
#        of a species' genomes share one sequence, or it has one genome).
#      - inter_sampled_all: distances from one randomly selected window in
#        sampled species to one representative window from every species.
#
#    Output is a long TSV: marker, metric, distance, count. Distances are
#    capped at MAX_DISTANCE ("MAX_DISTANCE or more").
#
#    Run through the queue -- reads a multi-GB SQLite database:
#      mqsub -t 1 -m 32 --hours 4 -- pixi run -e default python \
#        scripts/globdb_marker_distances.py
#
###############################################################################

import argparse
import collections
import logging
import os
import random
import sqlite3

import numpy as np

DEFAULT_DB = ("/mnt/hpccs01/work/microbiome/db/singlem/GlobDB_r232.metapackage_v4.smpkg/"
              "new_metapackage.sdb/otus.sqlite3")
DEFAULT_OUTPUT = "results/globdb_marker_distances/GlobDB_r232.tsv"
MAX_DISTANCE = 15
MAX_INTRA_SEQS = 30
QUERY_BATCH = 64


def species_of(taxonomy):
    return taxonomy.rsplit(";", 1)[-1].strip()


def encode(seqs):
    return np.frombuffer("".join(seqs).encode(), dtype=np.uint8).reshape(len(seqs), -1)


def marker_distances(conn, marker_id, n_inter_species, rng):
    rows = conn.execute(
        "SELECT o.sequence, t.taxonomy FROM otus o JOIN taxonomy t ON o.taxonomy_id = t.id "
        "WHERE o.marker_id = ?", (marker_id,))
    species_seqs = collections.defaultdict(set)
    n_rows = 0
    for seq, taxonomy in rows:
        species_seqs[species_of(taxonomy)].add(seq.replace("-", "N"))
        n_rows += 1

    lengths = collections.Counter(len(s) for v in species_seqs.values() for s in v)
    window = lengths.most_common(1)[0][0]
    for sp in species_seqs:
        species_seqs[sp] = {s for s in species_seqs[sp] if len(s) == window}
    species = sorted(sp for sp, v in species_seqs.items() if v)

    counts = collections.Counter()
    counts[("n_genome_windows", 0)] = n_rows
    counts[("n_species", 0)] = len(species)

    # Within species
    for sp in species:
        seqs = list(species_seqs[sp])
        if len(seqs) > MAX_INTRA_SEQS:
            seqs = rng.sample(seqs, MAX_INTRA_SEQS)
        if len(seqs) == 1:
            counts[("intra_diameter", 0)] += 1
            continue
        a = encode(seqs)
        d = (a[:, None, :] != a[None, :, :]).sum(-1)
        upper = d[np.triu_indices(len(seqs), 1)]
        for x in np.minimum(upper, MAX_DISTANCE):
            counts[("intra_pair", int(x))] += 1
        counts[("intra_diameter", int(min(upper.max(), MAX_DISTANCE)))] += 1

    # Compare sampled species against all species. Row batches bound the
    # temporary distance matrix, while the target array contains one
    # representative window per species.
    sampled_species = rng.sample(range(len(species)), min(n_inter_species, len(species)))
    sampled_seqs = [rng.choice(sorted(species_seqs[species[i]])) for i in sampled_species]
    target_seqs = [rng.choice(sorted(species_seqs[sp])) for sp in species]
    sampled_arr = encode(sampled_seqs)
    target_arr = encode(target_seqs)
    for start in range(0, len(sampled_seqs), QUERY_BATCH):
        stop = min(start + QUERY_BATCH, len(sampled_seqs))
        d = (sampled_arr[start:stop, None, :] != target_arr[None, :, :]).sum(-1)
        for row, query_i in zip(d, range(start, stop)):
            query_species_index = sampled_species[query_i]
            nearest = None
            for target_i, distance in enumerate(row):
                if target_i == query_species_index:
                    continue
                counts[("inter_sampled_all", int(min(distance, MAX_DISTANCE)))] += 1
                nearest = distance if nearest is None else min(nearest, distance)
            counts[("inter_sampled_nearest", int(min(nearest, MAX_DISTANCE)))] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--marker-prefix", default="S3.", help="Markers to include (default: S3.)")
    parser.add_argument("--max-marker-number", type=int, default=13,
                        help="Highest S3.N marker to include (default: 13, the markers analysed)")
    parser.add_argument("--inter-species", type=int, default=3000,
                        help="Species sampled per marker for sampled-to-all distance histogram")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    markers = [
        (marker_id, name) for marker_id, name in conn.execute("SELECT id, marker FROM markers")
        if name.startswith(args.marker_prefix)
        and int(name.split(".")[1]) <= args.max_marker_number
    ]
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    tmp = args.output + ".tmp"
    with open(tmp, "w") as out:
        out.write("marker\tmetric\tdistance\tcount\n")
        for marker_id, name in sorted(markers, key=lambda m: m[1]):
            logging.info(f"{name} ...")
            counts = marker_distances(conn, marker_id, args.inter_species, random.Random(args.seed))
            for (metric, distance), count in sorted(counts.items()):
                out.write(f"{name}\t{metric}\t{distance}\t{count}\n")
            out.flush()
            logging.info(f"{name}: {counts[('n_species', 0)]} species")
    os.replace(tmp, args.output)
    logging.info(f"Written to {args.output}")


if __name__ == "__main__":
    main()
