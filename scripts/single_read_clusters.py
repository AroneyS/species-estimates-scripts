#!/usr/bin/env python3

"""
Author: Samuel Aroney
Count clusters found with exactly one read hit in a sample scope (all samples or human
gut), per marker. A cluster with one total hit has one occurrence in one sample, so these
are the single-read f1 clusters that an error correction of the f1 count applies to.

Inputs: counts table (sequence, n_occurrences, sum_hits; gaps already converted to N, see
counts_table.py) and smafa clusters.tsv (member TAB representative, no header). Members
are joined to the counts, summed per representative, and clusters with no hits in the
scope are absent.
"""

import sys
import argparse
import logging
import polars as pl


def pipeline(counts, clusters):
    """counts: sequence, n_occurrences, sum_hits; clusters: member, rep. DataFrame or LazyFrame."""
    per_cluster = (
        clusters
        .join(counts, left_on="member", right_on="sequence", how="inner")
        .group_by("rep")
        .agg(total_hits=pl.col("sum_hits").sum(), total_occurrences=pl.col("n_occurrences").sum())
    )
    summary = per_cluster.select(
        n_clusters=pl.len(),
        n_single_read_clusters=(pl.col("total_hits") == 1).sum(),
        n_single_occurrence_clusters=(pl.col("total_occurrences") == 1).sum(),
    )
    return summary.collect() if isinstance(summary, pl.LazyFrame) else summary


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--counts-table", required=True, help="TSV: sequence, n_occurrences, sum_hits")
    parser.add_argument("--clusters", required=True, help="smafa clusters.tsv (member TAB rep)")
    parser.add_argument("--marker", required=True, help="Marker name, written to output")
    parser.add_argument("--scope", required=True, help="Sample scope (e.g. all, human_gut), written to output")
    parser.add_argument("--output", required=True, help="Output TSV")

    args = parser.parse_args(arguments)

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")
    logging.info(args)

    counts = pl.scan_csv(args.counts_table, separator="\t").select("sequence", "n_occurrences", "sum_hits")
    clusters = pl.scan_csv(args.clusters, separator="\t", has_header=False, new_columns=["member", "rep"]).with_columns(
        pl.col("member").str.replace_all("-", "N"))
    result = pipeline(counts, clusters).with_columns(scope=pl.lit(args.scope), marker=pl.lit(args.marker))
    result.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
