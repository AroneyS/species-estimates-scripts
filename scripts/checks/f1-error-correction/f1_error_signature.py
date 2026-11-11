#!/usr/bin/env python3

"""
Author: Samuel Aroney
Codon-level mismatch signature of rare SingleM window sequences vs their cluster
representative, stratified by sample-count class (f1, f2, f3, ...), single-read
status and representative abundance.

Sequencing errors predict ~uniform mismatches over codon positions, non-synonymous
excess and stop gains; real variation predicts third-position/synonymous excess.

Inputs: dedup table (seq_id, unknown, hits, n_occurrences, sum_hits, sequence) and
smafa clusters.tsv (representative TAB member, no header). Windows are 60 bp, in frame.
Codons containing N or gaps in either sequence are excluded from mismatch counts.
Standard genetic code (table 1); ignores table 4/25 TGA reassignments.
"""

import sys
import argparse
import logging
import itertools
import polars as pl

WINDOW = 60
N_CODONS = WINDOW // 3
REF_CLASS = "f100+"
STRATA = ["f_class", "single_read", "parent_class"]

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {"".join(c): a for c, a in zip(itertools.product(BASES, repeat=3), AMINO)}
TRANSITIONS = ["AG", "GA", "CT", "TC"]


def _expected_table() -> pl.DataFrame:
    """Per parent codon: fraction of the 9 single-nt substitutions that are syn/nonsyn/stop-gain."""
    rows = []
    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        syn = non = stop = 0
        for i in range(3):
            for b in "ACGT":
                if b == codon[i]:
                    continue
                new = CODON_TABLE[codon[:i] + b + codon[i + 1:]]
                if new == aa:
                    syn += 1
                elif new == "*":
                    stop += 1
                else:
                    non += 1
        rows.append([codon, syn / 9, non / 9, stop / 9])
    return pl.DataFrame(rows, schema={"rep_codon": str, "exp_syn": float, "exp_nonsyn": float, "exp_stop": float}, orient="row")


def f_class_expr(col: str = "n_occurrences") -> pl.Expr:
    c = pl.col(col)
    return (
        pl.when(c == 1).then(pl.lit("f1"))
        .when(c == 2).then(pl.lit("f2"))
        .when(c == 3).then(pl.lit("f3"))
        .when(c < 10).then(pl.lit("f4-9"))
        .when(c < 100).then(pl.lit("f10-99"))
        .otherwise(pl.lit(REF_CLASS))
    )


def parent_class_expr(col: str = "rep_n_occurrences") -> pl.Expr:
    c = pl.col(col)
    return (
        pl.when(c < 10).then(pl.lit("p<10"))
        .when(c < 100).then(pl.lit("p10-99"))
        .when(c < 1000).then(pl.lit("p100-999"))
        .otherwise(pl.lit("p1000+"))
    )


def _base_at_pos(codon: str) -> pl.Expr:
    return (
        pl.when(pl.col("pos") == 1).then(pl.col(codon).str.slice(0, 1))
        .when(pl.col("pos") == 2).then(pl.col(codon).str.slice(1, 1))
        .otherwise(pl.col(codon).str.slice(2, 1))
    )


def pipeline(pairs: pl.DataFrame) -> pl.DataFrame:
    """
    pairs: rep, member (60 bp), n_occurrences, sum_hits, rep_n_occurrences.
    Returns one row per stratum (f_class, single_read, parent_class) with counts.
    """
    pairs = pairs.filter(
        (pl.col("rep").str.len_chars() == WINDOW) & (pl.col("member").str.len_chars() == WINDOW)
    ).with_columns(
        f_class=f_class_expr(),
        single_read=pl.col("sum_hits") == pl.col("n_occurrences"),
        parent_class=parent_class_expr(),
        pair_id=pl.int_range(pl.len()),
        ambiguous=pl.col("rep").str.contains("[^ACGT]") | pl.col("member").str.contains("[^ACGT]"),
    )

    pair_counts = pairs.group_by(STRATA).agg(
        n_pairs=pl.len(),
        n_pairs_ambiguous=pl.col("ambiguous").sum(),
    )

    codons = pl.concat([
        pairs.select(
            "pair_id", *STRATA,
            codon_idx=pl.lit(i),
            rep_codon=pl.col("rep").str.slice(3 * i, 3),
            mem_codon=pl.col("member").str.slice(3 * i, 3),
        )
        for i in range(N_CODONS)
    ]).filter(
        (pl.col("rep_codon") != pl.col("mem_codon"))
        & pl.col("rep_codon").str.contains("^[ACGT]{3}$")
        & pl.col("mem_codon").str.contains("^[ACGT]{3}$")
    ).with_columns(
        n_diff=sum(
            (pl.col("rep_codon").str.slice(p, 1) != pl.col("mem_codon").str.slice(p, 1)).cast(pl.Int32)
            for p in range(3)
        ),
    )

    multi = codons.filter(pl.col("n_diff") > 1).group_by(STRATA).agg(n_codons_multi_nt=pl.len())

    single = codons.filter(pl.col("n_diff") == 1).with_columns(
        pos=pl.when(pl.col("rep_codon").str.slice(0, 1) != pl.col("mem_codon").str.slice(0, 1)).then(1)
        .when(pl.col("rep_codon").str.slice(1, 1) != pl.col("mem_codon").str.slice(1, 1)).then(2)
        .otherwise(3),
        rep_aa=pl.col("rep_codon").replace_strict(CODON_TABLE),
        mem_aa=pl.col("mem_codon").replace_strict(CODON_TABLE),
    ).filter(pl.col("rep_aa") != "*").with_columns(
        effect=pl.when(pl.col("rep_aa") == pl.col("mem_aa")).then(pl.lit("syn"))
        .when(pl.col("mem_aa") == "*").then(pl.lit("stop"))
        .otherwise(pl.lit("nonsyn")),
        change=pl.concat_str(_base_at_pos("rep_codon"), _base_at_pos("mem_codon")),
    ).with_columns(
        is_transition=pl.col("change").is_in(TRANSITIONS),
    ).join(_expected_table(), on="rep_codon", how="left")

    single_counts = single.group_by(STRATA).agg(
        n_single_nt=pl.len(),
        pos1=(pl.col("pos") == 1).sum(),
        pos2=(pl.col("pos") == 2).sum(),
        pos3=(pl.col("pos") == 3).sum(),
        syn=(pl.col("effect") == "syn").sum(),
        nonsyn=(pl.col("effect") == "nonsyn").sum(),
        stop=(pl.col("effect") == "stop").sum(),
        exp_syn=pl.col("exp_syn").sum(),
        exp_nonsyn=pl.col("exp_nonsyn").sum(),
        exp_stop=pl.col("exp_stop").sum(),
        transitions=pl.col("is_transition").sum(),
    )

    out = pair_counts.join(single_counts, on=STRATA, how="left").join(multi, on=STRATA, how="left")
    fill = [c for c in out.columns if c not in STRATA]
    return out.with_columns(pl.col(fill).fill_null(0)).sort(STRATA)


def load_pairs(dedup_table: str, clusters: str, sample_fraction: float, seed: int) -> pl.DataFrame:
    """Sample non-representative members from clusters and attach counts and parent abundance."""
    dedup = pl.scan_csv(dedup_table, separator="\t").select("sequence", "n_occurrences", "sum_hits")
    clu = pl.scan_csv(clusters, separator="\t", has_header=False, new_columns=["rep", "member"])

    modulus = 1_000_000
    keep = (pl.col("sequence").hash(seed) % modulus) < int(sample_fraction * modulus)
    members = (
        dedup.filter((f_class_expr() == REF_CLASS) | keep)
        .rename({"sequence": "member"})
    )
    pairs = (
        clu.filter(pl.col("rep") != pl.col("member"))
        .join(members, on="member", how="inner")
        .collect(engine="streaming")
    )
    reps = pairs.select("rep").unique().lazy()
    rep_counts = (
        dedup.join(reps, left_on="sequence", right_on="rep", how="semi")
        .select(rep="sequence", rep_n_occurrences="n_occurrences")
        .collect(engine="streaming")
    )
    pairs = pairs.join(rep_counts, on="rep", how="inner").filter(
        pl.col("rep_n_occurrences") > pl.col("n_occurrences")
    )
    logging.info("pairs with a more abundant representative: %d", pairs.height)
    return pairs


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--dedup-table", required=True, help="TSV from deduplicated.fasta (with header)")
    parser.add_argument("--clusters", required=True, help="smafa clusters.tsv (rep TAB member)")
    parser.add_argument("--marker", required=True, help="Marker name, written to output")
    parser.add_argument("--sample-fraction", type=float, default=0.02,
                        help="Hash-sampled fraction of members with n_occurrences < 100 (f100+ always kept)")
    parser.add_argument("--seed", type=int, default=0)
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

    pairs = load_pairs(args.dedup_table, args.clusters, args.sample_fraction, args.seed)
    result = pipeline(pairs).with_columns(
        marker=pl.lit(args.marker), sample_fraction=pl.lit(args.sample_fraction))
    result.write_csv(args.output, separator="\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
