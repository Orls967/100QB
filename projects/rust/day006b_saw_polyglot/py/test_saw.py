#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path

from saw_analysis import (
    A001411_18, brute_force_saw, check_figures, fit_ratio, fit_scaling, make_figures, validate_saw,
)


class TestSawB(unittest.TestCase):
    def test_brute_force(self) -> None:
        self.assertEqual(brute_force_saw(8)[1:], A001411_18[:8])

    def test_validator(self) -> None:
        self.assertTrue(validate_saw([[0, 0], [1, 0], [1, 1]], 2))
        self.assertFalse(validate_saw([[0, 0], [1, 0], [0, 0]], 2))
        self.assertFalse(validate_saw([[0, 0], [2, 0]], 1))

    def test_synthetic_regressions(self) -> None:
        c = [1.0]
        for n in range(1, 20):
            c.append(c[-1] * 2.638 * (1.0 + 0.34375 / n))
        rf = fit_ratio(c, 11, 18)
        self.assertAlmostEqual(rf["mu_hat"], 2.638, places=8)
        self.assertAlmostEqual(rf["gamma_hat"], 1.34375, places=8)
        ns = [50, 100, 200, 400]
        sf = fit_scaling(ns, [0.77 * (n ** 1.5) for n in ns])
        self.assertAlmostEqual(sf["nu_hat"], 0.75, places=8)

    def test_figures_xml_valid(self) -> None:
        cnts = [1] + A001411_18
        data = {
            "nmax": 18, "go_counts": cnts,
            "rust_single": {"heatmap_n": 4, "heatmap": [[1] * 9 for _ in range(9)]},
            "pivots": [{"n": n, "r2_mean": 0.77 * (n ** 1.5)} for n in (50, 100, 200, 400)],
            "samples_n300": {"samples": [[[i, 0] for i in range(301)] for _ in range(100)]},
        }
        rf = fit_ratio(cnts, 11, 18)
        sf = fit_scaling([50, 100, 200, 400], [p["r2_mean"] for p in data["pivots"]])
        with tempfile.TemporaryDirectory() as td:
            figs = make_figures(data, rf, sf, Path(td), prefix="test")
            self.assertTrue(check_figures(figs))
            self.assertEqual(len(figs), 4)


if __name__ == "__main__":
    unittest.main()
