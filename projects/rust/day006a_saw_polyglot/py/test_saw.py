#!/usr/bin/env python3
"""Unit tests for Python SAW reference, validator, regressions, and visual artifact checks."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from saw_analysis import (
    MEMORY_A001411_N1_20,
    brute_force_saw_counts,
    estimate_mu_gamma_from_counts,
    estimate_nu_from_pivots,
    generate_all_figures,
    generate_interactive_html,
    validate_saw_walk,
    verify_visual_artifacts,
)


class TestSawAnalysis(unittest.TestCase):
    def test_brute_force_small_n(self) -> None:
        counts = brute_force_saw_counts(9)
        self.assertEqual(counts[0], 1)
        self.assertEqual(counts[1:10], MEMORY_A001411_N1_20[:9])

    def test_validator_detects_valid_collision_and_jump(self) -> None:
        valid = [[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1]]
        ok, msg = validate_saw_walk(valid, expected_n=4)
        self.assertTrue(ok, msg)

        collided = [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]
        ok_c, _ = validate_saw_walk(collided, expected_n=4)
        self.assertFalse(ok_c)

        jumped = [[0, 0], [2, 0], [2, 1]]
        ok_j, _ = validate_saw_walk(jumped, expected_n=2)
        self.assertFalse(ok_j)

    def test_loglog_scaling_regression_synthetic(self) -> None:
        true_amp = 0.78
        true_nu = 0.75
        ns = [50, 100, 200, 400, 800]
        r2 = [true_amp * (n ** (2.0 * true_nu)) for n in ns]
        fit = estimate_nu_from_pivots(ns, r2)
        self.assertAlmostEqual(float(fit["nu_hat"]), true_nu, places=10)
        self.assertAlmostEqual(float(fit["amplitude_hat"]), true_amp, places=10)

    def test_ratio_extrapolation_synthetic_1_over_n_series(self) -> None:
        # Construct exact synthetic floats with r_n = mu * (1 + (gamma - 1)/n)
        mu_true = 2.638
        gamma_true = 1.34375
        c = [1.0]
        for n in range(1, 25):
            r_n = mu_true * (1.0 + (gamma_true - 1.0) / n)
            c.append(c[-1] * r_n)
        fit = estimate_mu_gamma_from_counts(c, 16, 24)
        self.assertAlmostEqual(float(fit["mu_hat"]), mu_true, places=9)
        self.assertAlmostEqual(float(fit["gamma_hat"]), gamma_true, places=9)

    def test_figure_and_html_generation_and_validation(self) -> None:
        pivots = []
        for n in (50, 100, 200, 400, 800):
            samples = [[[i, 0] for i in range(n + 1)] for _ in range(120 if n == 400 else 4)]
            pivots.append({"n": n, "seed": 1, "r2_mean": 0.77 * (n ** 1.5), "r2_se": 0.5, "samples": samples})
        counts = [1] + MEMORY_A001411_N1_20[:16]
        hm_side = 2 * 8 + 1
        hm = [[1 for _ in range(hm_side)] for _ in range(hm_side)]
        bundle = {
            "nmax": 16,
            "prefix_k": 6,
            "num_prefix_tasks": 195,
            "base_seed": 20261124,
            "rust_single": {"nmax": 16, "counts": counts, "heatmap_n": 8, "heatmap": hm},
            "rust_single_sec": 0.01,
            "go_sharded_counts": counts,
            "worker_timings": [
                {"workers": 1, "elapsed_sec": 0.04, "speedup": 1.0, "counts_match": True},
                {"workers": 2, "elapsed_sec": 0.022, "speedup": 1.82, "counts_match": True},
            ],
            "worker_invariant": True,
            "pivots": pivots,
        }
        r_fit = estimate_mu_gamma_from_counts(counts, 8, 16)
        s_fit = estimate_nu_from_pivots([50, 100, 200, 400, 800], [p["r2_mean"] for p in pivots])
        with tempfile.TemporaryDirectory() as tmp:
            fig_dir = Path(tmp)
            figs = generate_all_figures(bundle, r_fit, s_fit, fig_dir)
            html_path = fig_dir / "day006a_saw_viewer.html"
            generate_interactive_html(bundle, html_path)
            ok, issues = verify_visual_artifacts(figs, html_path)
            self.assertTrue(ok, f"artifact issues: {issues}")
            self.assertEqual(len(figs), 5)


if __name__ == "__main__":
    unittest.main()
