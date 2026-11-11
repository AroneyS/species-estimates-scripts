import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
#!/usr/bin/env python3

import unittest
import polars as pl
from polars.testing import assert_frame_equal
from counts_table import pipeline

INPUT_COLUMNS = {"sequence": str, "n_occurrences": int, "sum_hits": int}


class Tests(unittest.TestCase):
    def assertDataFrameEqual(self, a, b):
        assert_frame_equal(a, b, check_dtype=False, check_row_order=False)

    def test_sums_per_sequence(self):
        input = pl.DataFrame([["AAA", 1, 3], ["AAA", 1, 2], ["CCC", 4, 4]], schema=INPUT_COLUMNS, orient="row")
        expected = pl.DataFrame([["AAA", 2, 5], ["CCC", 4, 4]], schema=INPUT_COLUMNS, orient="row")
        self.assertDataFrameEqual(expected, pipeline(input))

    def test_empty(self):
        self.assertEqual(pipeline(pl.DataFrame(schema=INPUT_COLUMNS)).height, 0)


if __name__ == "__main__":
    unittest.main()
