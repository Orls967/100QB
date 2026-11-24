#!/usr/bin/env python3
import math
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from percolation import (
    binom_weights, canonical_eval, check_figures, eval_exact_poly,
    exact_enumeration, make_all_figures, python_newman_ziff,
)


class TestPercolationC(unittest.TestCase):
    def test_exact_self_dual_half(self) -> None:
        self.assertEqual(exact_enumeration(2, 3)["p_half"], Fraction(1, 2))
        self.assertNotEqual(exact_enumeration(3, 2)["p_half"], Fraction(1, 2))

    def test_binomial_weights_sum_to_one(self) -> None:
        for p in (0.0, 0.2, 0.5, 0.85, 1.0):
            self.assertAlmostEqual(float(binom_weights(127, p).sum()), 1.0, places=11)

    def test_nz_python_vs_exact(self) -> None:
        ex = exact_enumeration(2, 3)
        hist = python_newman_ziff(2, 3, 10000, 20261124)
        ps = [round(0.1 * k, 1) for k in range(1, 10)]
        vals, _ = canonical_eval(hist, ps)
        for k, pv in zip(range(1, 10), vals):
            pe = float(eval_exact_poly(ex["counts"], Fraction(k, 10)))
            self.assertLessEqual(abs(pv - pe), 3.5 * math.sqrt(pe * (1.0 - pe) / 10000.0))

    def test_figures_valid(self) -> None:
        sdoc = {"rows": 8, "cols": 9, "open_bonds": [[0, 1], [1, 2], [2, 11]]}
        curves = [{"L": l, "p_grid": [0.4, 0.5, 0.6], "p_cross": [0.2, 0.5, 0.8], "dp_half": l ** 0.75, "mean_largest": l ** 1.89} for l in (8, 16)]
        with tempfile.TemporaryDirectory() as td:
            figs = make_all_figures(sdoc, curves, {"nu_hat": 1.333}, {"d_hat": 1.896}, Path(td), prefix="t")
            self.assertTrue(check_figures(figs))


if __name__ == "__main__":
    unittest.main()
