#!/usr/bin/env python3

###############################################################################
#
#    Summarise a species_estimation output folder (e.g.
#    results/species_estimation/20260704) into a TSV of:
#        marker  sample_count  species_count
#    where species_count is the number of species (clusters, i.e. rows of
#    {marker}/{marker}.tsv) whose n_samples equals sample_count, for each
#    marker directory found in the output folder.
#
###############################################################################

import argparse
import csv
import logging
import os
import re
import sys
from collections import Counter

MARKER_DIR_RE = re.compile(r"^S\d+\.")


def find_marker_tsvs(input_dir):
    """
    Return a sorted list of (marker_name, tsv_path) for every marker
    subdirectory of input_dir that has a {marker}.tsv file.
    """
    markers = []
    for name in sorted(os.listdir(input_dir)):
        marker_dir = os.path.join(input_dir, name)
        if not os.path.isdir(marker_dir) or not MARKER_DIR_RE.match(name):
            continue
        tsv_path = os.path.join(marker_dir, f"{name}.tsv")
        if os.path.isfile(tsv_path):
            markers.append((name, tsv_path))
        else:
            logging.warning(f"No {name}.tsv found in {marker_dir} -- skipping.")
    return markers


def count_sample_counts(tsv_path):
    """
    Stream a marker's {marker}.tsv and return a Counter mapping
    n_samples -> number of species (rows) with that n_samples value.
    """
    counts = Counter()
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is not None and "n_samples" not in reader.fieldnames:
            raise ValueError(f"{tsv_path} has no 'n_samples' column -- cannot summarise.")
        for row in reader:
            counts[int(row["n_samples"])] += 1
    return counts


def summarise(input_dir):
    """Yield (marker, sample_count, species_count) rows for all markers."""
    for marker_name, tsv_path in find_marker_tsvs(input_dir):
        logging.info(f"Processing {marker_name} ...")
        counts = count_sample_counts(tsv_path)
        for sample_count in sorted(counts):
            species_count = counts[sample_count]
            if sample_count > 0 and species_count > 0:
                yield marker_name, sample_count, species_count


def main():
    parser = argparse.ArgumentParser(
        description="Summarise a species_estimation output folder into a TSV of "
                     "marker, sample_count, species_count."
    )
    parser.add_argument(
        "input_dir",
        help="Path to a species_estimation output folder "
             "(e.g. results/species_estimation/20260704)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output TSV path (default: <input_dir>/sample_count_summary.tsv)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    output_path = args.output or os.path.join(args.input_dir, "sample_count_summary.tsv")

    with open(output_path, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["marker", "sample_count", "species_count"])
        for row in summarise(args.input_dir):
            writer.writerow(row)

    logging.info(f"Summary written to {output_path}")


if __name__ == "__main__":
    main()
