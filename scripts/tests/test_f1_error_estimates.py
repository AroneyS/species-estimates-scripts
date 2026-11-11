import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
#!/usr/bin/env python3

import unittest
import polars as pl
from f1_error_estimates import pipeline

COLUMNS = {
    "scope": str, "marker": str, "f_class": str, "single_read": bool, "parent_class": str,
    "n_pairs": int, "n_pairs_ambiguous": int, "n_single_nt": int, "pos1": int, "pos2": int, "pos3": int,
    "syn": int, "nonsyn": int, "stop": int, "exp_syn": float, "exp_nonsyn": float, "exp_stop": float,
    "transitions": int, "n_codons_multi_nt": int, "sample_fraction": float,
}


def row(f_class, single_read, n, pos3, syn, stop, exp_syn, exp_stop, parent="p<10"):
    return ["all", "m1", f_class, single_read, parent, n, 0, n, 0, 0, pos3, syn, n - syn - stop, stop,
            exp_syn, 0.0, exp_stop, 0, 0, 1.0]


class Tests(unittest.TestCase):
    def test_half_error_mixture(self):
        # reference: pos3 0.9, syn 0.9, stop 0; error expectation 1/3, syn 0.3, stop 0.03
        # class: 50% error -> pos3 0.6167, syn 0.6, stop 0.015
        sig = pl.DataFrame([
            row("f100+", False, 1000, 900, 900, 0, 300.0, 30.0),
            row("f1", True, 1000, 617, 600, 15, 300.0, 30.0),
        ], schema=COLUMNS, orient="row")
        out = pipeline(sig).filter(pl.col("parent_class") == "all")
        est = dict(zip(out["statistic"], out["error_fraction"]))
        self.assertAlmostEqual(est["codon_position"], 0.5, places=2)
        self.assertAlmostEqual(est["synonymous"], 0.5, places=2)
        self.assertAlmostEqual(est["stop_codon"], 0.5, places=2)
        self.assertEqual(set(out["f_class"]), {"f1"})

    def test_reference_pooled_and_parent_classes(self):
        sig = pl.DataFrame([
            row("f100+", False, 500, 450, 450, 0, 150.0, 15.0, "p100-999"),
            row("f100+", False, 500, 450, 450, 0, 150.0, 15.0, "p1000+"),
            row("f1", True, 1000, 617, 600, 15, 300.0, 30.0),
        ], schema=COLUMNS, orient="row")
        out = pipeline(sig)
        self.assertEqual(set(out["parent_class"]), {"all", "p<10"})

    def test_all_scope_reference_for_human_gut(self):
        hg = lambda *a, **k: ["human_gut"] + row(*a, **k)[1:]
        sig = pl.DataFrame([
            row("f100+", False, 1000, 900, 900, 0, 300.0, 30.0),
            hg("f100+", False, 1000, 700, 700, 0, 300.0, 30.0),
            hg("f1", True, 1000, 517, 500, 15, 300.0, 30.0),
        ], schema=COLUMNS, orient="row")
        same = pipeline(sig).filter((pl.col("scope") == "human_gut") & (pl.col("statistic") == "codon_position")
                                    & (pl.col("parent_class") == "all"))
        ref_all = pipeline(sig, reference_scope="all").filter(
            (pl.col("scope") == "human_gut") & (pl.col("statistic") == "codon_position")
            & (pl.col("parent_class") == "all"))
        self.assertAlmostEqual(same["error_fraction"][0], 0.5, places=2)   # (0.7-0.517)/(0.7-1/3)
        self.assertAlmostEqual(ref_all["error_fraction"][0], 0.675, places=2)  # (0.9-0.517)/(0.9-1/3)
        self.assertEqual(ref_all["reference"][0], "all")

    def test_skips_classes_without_single_nt_changes(self):
        sig = pl.DataFrame([
            row("f100+", False, 1000, 900, 900, 0, 300.0, 30.0),
            row("f1", True, 0, 0, 0, 0, 0.0, 0.0),
        ], schema=COLUMNS, orient="row")
        self.assertEqual(pipeline(sig).height, 0)


if __name__ == "__main__":
    unittest.main()
