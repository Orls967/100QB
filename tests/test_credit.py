"""Unit test untuk modul qb100.credit."""

from __future__ import annotations

import numpy as np

from qb100 import credit


def test_exact_pmf_mean_and_bivariate_variance():
    for case, expected_mean in (("A", 2.0), ("B", 6.0)):
        exp = credit.make_exposures(case)
        pmf = credit.exact_loss_pmf(exp)
        assert abs(float(np.sum(pmf)) - 1.0) < 1e-12

        metrics = credit.exact_risk_metrics(pmf)
        assert abs(metrics["mean"] - expected_mean) < 1e-12

        var_formula = credit.bivariate_loss_variance(exp)
        rel_err = abs(metrics["var"] - var_formula) / var_formula
        assert rel_err < 1e-6


def test_homogeneous_pmf_matches_quad():
    exp_a = credit.make_exposures("A")
    pmf_a = credit.exact_loss_pmf(exp_a)
    k_vals = (0, 1, 2, 5, 10, 20)
    quad_vals = credit.exact_pmf_quad_homogeneous(k_vals)
    assert np.allclose(pmf_a[list(k_vals)], quad_vals, atol=1e-10)


def test_lhp_approaches_exact_for_large_n():
    n_large = 5000
    exp_large = credit.make_exposures("A", n_obligors=n_large)
    pmf = credit.exact_loss_pmf(exp_large)
    cdf_exact = np.cumsum(pmf)
    x_grid = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.30])
    idx = np.round(x_grid * n_large).astype(int)
    tail_exact = 1.0 - cdf_exact[idx - 1]
    tail_lhp = credit.lhp_tail_prob(x_grid)
    assert np.max(np.abs(tail_exact - tail_lhp)) < 5e-3


def test_three_estimators_unbiased_reproducible_and_is_reduces_se():
    exp_b = credit.make_exposures("B")
    pmf_b = credit.exact_loss_pmf(exp_b)
    thr = 90  # 30% dari total eksposur 300
    exact_p = float(np.sum(pmf_b[thr:]))
    spline = credit.build_conditional_tail_interpolator(exp_b, [thr])[thr]

    res_naive = credit.estimate_naive_mc(exp_b, thr, n_samples=20_000, chunk_size=10_000, seed=42, exact_val=exact_p)
    res_cond = credit.estimate_conditional_mc(spline, n_samples=20_000, seed=42, exact_val=exact_p)
    theta, _, _ = credit.find_optimal_is_theta(spline, n_pilot=5_000, seed=999)
    res_is1 = credit.estimate_importance_sampling_mc(spline, theta, n_samples=20_000, seed=42, exact_val=exact_p)
    res_is2 = credit.estimate_importance_sampling_mc(spline, theta, n_samples=20_000, seed=42, exact_val=exact_p)

    assert res_naive["z_score"] < 3.0
    assert res_cond["z_score"] < 3.0
    assert res_is1["z_score"] < 3.0
    assert res_is1["se"] < res_naive["se"]
    assert res_is1["estimate"] == res_is2["estimate"]
