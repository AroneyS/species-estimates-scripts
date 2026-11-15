import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
#!/usr/bin/env python3

import unittest
import polars as pl
from polars.testing import assert_frame_equal
from single_read_clusters import pipeline

COUNTS = {"sequence": str, "n_occurrences": int, "sum_hits": int}
CLUSTERS = {"member": str, "rep": str}
OUTPUT = {"n_clusters": int, "n_single_read_clusters": int, "n_single_occurrence_clusters": int}


class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False, check_row_order=False)

    def test_single_read_clusters(self):
        counts = pl.DataFrame([
            ["A", 1, 1],   # cluster A: one hit -> single-read
            ["B", 1, 3],   # cluster B: one occurrence, three hits
            ["C", 1, 1], ["D", 1, 1],  # cluster C: two members, one hit each -> two hits
        ], schema=COUNTS, orient="row")
        clusters = pl.DataFrame([
            ["A", "A"], ["B", "B"], ["C", "C"], ["D", "C"], ["E", "E"],  # E has no hits in scope
        ], schema=CLUSTERS, orient="row")
        expected = pl.DataFrame([[3, 1, 2]], schema=OUTPUT, orient="row")
        self.assertDataFrameEqual(expected, pipeline(counts, clusters))

    def test_lazy_input(self):
        counts = pl.DataFrame([["A", 1, 1]], schema=COUNTS, orient="row").lazy()
        clusters = pl.DataFrame([["A", "A"]], schema=CLUSTERS, orient="row").lazy()
        expected = pl.DataFrame([[1, 1, 1]], schema=OUTPUT, orient="row")
        self.assertDataFrameEqual(expected, pipeline(counts, clusters))

    def test_empty(self):
        counts = pl.DataFrame(schema=COUNTS)
        clusters = pl.DataFrame([["A", "A"]], schema=CLUSTERS, orient="row")
        expected = pl.DataFrame([[0, 0, 0]], schema=OUTPUT, orient="row")
        self.assertDataFrameEqual(expected, pipeline(counts, clusters))


if __name__ == "__main__":
    unittest.main()
