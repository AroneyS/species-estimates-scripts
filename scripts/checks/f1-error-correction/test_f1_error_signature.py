#!/usr/bin/env python3

import unittest
import polars as pl
from polars.testing import assert_frame_equal
from f1_error_signature import pipeline

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

    def test_empty(self):
        out = pipeline(pl.DataFrame(schema=INPUT_COLUMNS))
        self.assertEqual(out.height, 0)


if __name__ == "__main__":
    unittest.main()
