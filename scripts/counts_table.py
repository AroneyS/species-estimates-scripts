#!/usr/bin/env python3

"""
Author: Samuel Aroney
Aggregate per-occurrence (or per-unique) rows into a counts table per sequence:
sequence, n_occurrences, sum_hits. Input TSV has a header with columns
sequence, n_occurrences, sum_hits (awk-extracted from the dedup FASTA or from human gut
collated rows; gaps already converted to N to match smafa's dash-to-N transform).
Sequences differing only by gap vs N are merged.
"""

import sys
import argparse
import logging
import polars as pl


def pipeline(df):
    """Works on DataFrame or LazyFrame."""
    return (
        df.group_by("sequence")
        .agg(n_occurrences=pl.col("n_occurrences").sum(), sum_hits=pl.col("sum_hits").sum())
    )


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--input", required=True, help="Input TSV")
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

    lf = pl.scan_csv(args.input, separator="\t").select("sequence", "n_occurrences", "sum_hits")
    pipeline(lf).sink_csv(args.output, separator="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
