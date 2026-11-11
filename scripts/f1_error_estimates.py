#!/usr/bin/env python3

"""
Author: Samuel Aroney
Estimate the sequencing-error fraction among rare-sequence mismatches from the collated
signature table, as a two-component mixture of errors and real variation.

Real variation is taken from the reference class (f100+, pooled across single-read status
and parent abundance) within each scope and marker. For each class (scope, marker, f_class,
single_read, parent_class), with observed fraction p_obs of a statistic, reference p_ref and
error expectation p_err:

    error_fraction = (p_ref - p_obs) / (p_ref - p_err)

Statistics (error expectation under uniform random substitution):
  codon_position: fraction at codon position 3; expectation 1/3
  synonymous:     fraction synonymous; expectation from the class's own parent codons
  stop_codon:     fraction stop-gain; expectation from the class's own parent codons
By default the reference is from the same scope. With --reference-scope (e.g. all), the
reference for every scope is that scope's f100+ pairs for the same marker; this avoids a
thin or error-contaminated reference in a deep-sampled scope such as human gut.
Standard errors are binomial on p_obs only (reference uncertainty and non-independence
of codons within sequences and studies are ignored).
"""

import sys
import argparse
import logging
import polars as pl

KEYS = ["scope", "marker"]
CLASS_KEYS = KEYS + ["f_class", "single_read", "parent_class"]
REF_CLASS = "f100+"

# statistic: (count column, expectation expression)
STATISTICS = {
    "codon_position": ("pos3", pl.lit(1 / 3)),
    "synonymous": ("syn", pl.col("exp_syn") / pl.col("n_single_nt")),
    "stop_codon": ("stop", pl.col("exp_stop") / pl.col("n_single_nt")),
}


def pipeline(signature: pl.DataFrame, reference_scope: str | None = None) -> pl.DataFrame:
    """signature: collated signature table (scope, marker, strata, counts)."""
    sums = [c for c in signature.columns if c not in CLASS_KEYS + ["sample_fraction"]]
    pooled = signature.group_by(
        KEYS + ["f_class", "single_read"]
    ).agg(pl.col(sums).sum()).with_columns(parent_class=pl.lit("all"))
    per_parent = signature.select(CLASS_KEYS + sums)
    classes = pl.concat([pooled.select(per_parent.columns), per_parent])

    ref_keys = KEYS if reference_scope is None else ["marker"]
    ref_rows = signature.filter(pl.col("f_class") == REF_CLASS)
    if reference_scope is not None:
        ref_rows = ref_rows.filter(pl.col("scope") == reference_scope)
    ref = (
        ref_rows.group_by(ref_keys).agg(pl.col(sums).sum())
        .rename({c: f"ref_{c}" for c in sums})
    )
    data = classes.filter(pl.col("f_class") != REF_CLASS).join(ref, on=ref_keys, how="left").filter(
        pl.col("ref_n_single_nt").is_not_null()
    )

    out = []
    for name, (count, expected) in STATISTICS.items():
        out.append(
            data.filter(pl.col("n_single_nt") > 0).with_columns(
                p_obs=pl.col(count) / pl.col("n_single_nt"),
                p_ref=pl.col(f"ref_{count}") / pl.col("ref_n_single_nt"),
                p_err=expected,
            ).with_columns(
                error_fraction=(pl.col("p_ref") - pl.col("p_obs")) / (pl.col("p_ref") - pl.col("p_err")),
                error_fraction_se=(
                    (pl.col("p_obs") * (1 - pl.col("p_obs")) / pl.col("n_single_nt")).sqrt()
                    / (pl.col("p_ref") - pl.col("p_err")).abs()
                ),
                statistic=pl.lit(name),
                reference=pl.lit(reference_scope or "same_scope"),
            ).select(
                CLASS_KEYS + ["statistic", "reference", "n_pairs", "n_single_nt", "p_obs", "p_ref", "p_err",
                              "error_fraction", "error_fraction_se"]
            )
        )
    return pl.concat(out).sort(CLASS_KEYS + ["statistic"])


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--signature", required=True, help="Collated signature TSV")
    parser.add_argument("--reference-scope", default=None,
                        help="Take the f100+ reference from this scope for all scopes (default: same scope)")
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

    pipeline(pl.read_csv(args.signature, separator="\t"), args.reference_scope).write_csv(args.output, separator="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
