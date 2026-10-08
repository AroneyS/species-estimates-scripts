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
#      - inter_nearest: for a random sample of species, the distance from one
#        of its sequences to the nearest sequence of any other species
#        (0 when the sequence is shared with another species).
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
QUERY_BATCH = 8


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

    # Nearest other species. owner[i] is the species index of sequence i, or
    # -1 when the sequence is shared by several species.
    seq_owner = {}
    for i, sp in enumerate(species):
        for s in species_seqs[sp]:
            seq_owner[s] = i if seq_owner.get(s, i) == i else -1
    all_seqs = list(seq_owner)
    owner = np.array([seq_owner[s] for s in all_seqs])
    all_arr = encode(all_seqs)
    seq_index = {s: i for i, s in enumerate(all_seqs)}

    query_species = rng.sample(range(len(species)), min(n_inter_species, len(species)))
    for start in range(0, len(query_species), QUERY_BATCH):
        batch = query_species[start:start + QUERY_BATCH]
        query_seqs = [rng.choice(sorted(species_seqs[species[i]])) for i in batch]
        d = (encode(query_seqs)[:, None, :] != all_arr[None, :, :]).sum(-1)
        for row, sp_i, q in zip(d, batch, query_seqs):
            if owner[seq_index[q]] == -1:
                nearest = 0
            else:
                nearest = row[owner != sp_i].min()
            counts[("inter_nearest", int(min(nearest, MAX_DISTANCE)))] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--marker-prefix", default="S3.", help="Markers to include (default: S3.)")
    parser.add_argument("--max-marker-number", type=int, default=13,
                        help="Highest S3.N marker to include (default: 13, the markers analysed)")
    parser.add_argument("--inter-species", type=int, default=3000,
                        help="Species sampled per marker for nearest other-species distance")
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
