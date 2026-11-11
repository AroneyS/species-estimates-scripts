#!/usr/bin/env python3

import unittest
import polars as pl
from polars.testing import assert_frame_equal
import os
import tempfile
from f1_error_signature import pipeline, combine, build_pairs, summarise_pairs

INPUT_COLUMNS = {
    "rep": str, "member": str, "n_occurrences": int, "sum_hits": int, "rep_n_occurrences": int,
}

REP = "ATG" + "GCT" * 19  # Met, then 19 x Ala


def mutate(seq, codon_idx, new_codon):
    return seq[:3 * codon_idx] + new_codon + seq[3 * codon_idx + 3:]


class Tests(unittest.TestCase):
    def test_single_nt_effects(self):
        # GCT->GCC syn at pos 3; GCT->ACT nonsyn at pos 1; GCT->GAT nonsyn pos 2 (Ala->Asp)
        input = pl.DataFrame([
            [REP, mutate(REP, 1, "GCC"), 1, 1, 50],
            [REP, mutate(REP, 2, "ACT"), 1, 1, 50],
            [REP, mutate(REP, 3, "GAT"), 1, 1, 50],
        ], schema=INPUT_COLUMNS, orient="row")
        out = pipeline(input)
        self.assertEqual(out.height, 1)
        row = out.row(0, named=True)
        self.assertEqual((row["f_class"], row["single_read"], row["parent_class"]), ("f1", True, "p10-99"))
        self.assertEqual(row["n_pairs"], 3)
        self.assertEqual((row["pos1"], row["pos2"], row["pos3"]), (1, 1, 1))
        self.assertEqual((row["syn"], row["nonsyn"], row["stop"]), (1, 2, 0))
        self.assertEqual(row["transitions"], 2)  # T->C and G->A; C->A is a transversion

    def test_stop_gain_and_expected(self):
        # TAC->TAA not possible from Ala; use TGG (Trp) -> TGA stop
        rep = "ATG" + "TGG" * 19
        input = pl.DataFrame([[rep, mutate(rep, 4, "TGA"), 2, 3, 5]], schema=INPUT_COLUMNS, orient="row")
        row = pipeline(input).row(0, named=True)
        self.assertEqual((row["f_class"], row["single_read"]), ("f2", False))
        self.assertEqual(row["stop"], 1)
        self.assertAlmostEqual(row["exp_stop"], 2 / 9)  # TGG -> TAG, TGA

    def test_ambiguous_codons_excluded(self):
        member = mutate(REP, 1, "GCN")
        input = pl.DataFrame([[REP, member, 1, 1, 50]], schema=INPUT_COLUMNS, orient="row")
        row = pipeline(input).row(0, named=True)
        self.assertEqual(row["n_pairs_ambiguous"], 1)
        self.assertEqual(row["n_single_nt"], 0)

    def test_multi_nt_codon(self):
        input = pl.DataFrame([[REP, mutate(REP, 1, "GGG"), 1, 1, 50]], schema=INPUT_COLUMNS, orient="row")
        row = pipeline(input).row(0, named=True)
        self.assertEqual((row["n_codons_multi_nt"], row["n_single_nt"]), (1, 0))

    def test_combine_sums_strata(self):
        input = pl.DataFrame([[REP, mutate(REP, 1, "GCC"), 1, 1, 50]], schema=INPUT_COLUMNS, orient="row")
        out = pipeline(input)
        combined = combine(pl.concat([out, out]))
        self.assertEqual(combined.row(0, named=True)["n_pairs"], 2)
        self.assertEqual(combined.row(0, named=True)["pos3"], 2)

    def test_build_pairs_normalises_gaps_and_filters_abundance(self):
        gapped = mutate(REP, 1, "G-T")
        normalised = gapped.replace("-", "N")
        member2 = mutate(REP, 2, "GCC")
        with tempfile.TemporaryDirectory() as d:
            counts = os.path.join(d, "counts.tsv")
            clusters = os.path.join(d, "clusters.tsv")
            out = os.path.join(d, "pairs.parquet")
            pl.DataFrame({
                "sequence": [REP, normalised, member2],
                "n_occurrences": [50, 1, 60],  # member2 more abundant than its rep: dropped
                "sum_hits": [100, 1, 100],
            }).write_csv(counts, separator="\t")
            with open(clusters, "w") as fh:
                fh.write(f"{REP}\t{REP}\n{REP}\t{normalised}\n{REP}\t{member2}\n")
            self.assertEqual(build_pairs(counts, clusters, out), 1)
            row = pl.read_parquet(out).row(0, named=True)
            self.assertEqual(row["member"], normalised)
            self.assertEqual(row["rep_n_occurrences"], 50)
            summary = summarise_pairs(out, batch_rows=1)
            self.assertEqual(summary["n_pairs_ambiguous"].sum(), 1)

    def test_empty(self):
        out = pipeline(pl.DataFrame(schema=INPUT_COLUMNS))
        self.assertEqual(out.height, 0)


if __name__ == "__main__":
    unittest.main()
