"""Unit test untuk modul qb100.covariance."""

from __future__ import annotations

import numpy as np

from qb100 import covariance as cov


def test_factor_model_is_symmetric_positive_definite():
    sigma_true, mu = cov.generate_factor_model(n_assets=50, n_factors=3, seed=42)
    assert sigma_true.shape == (50, 50)
    assert mu.shape == (50,)
    assert np.allclose(sigma_true, sigma_true.T, atol=1e-14)
    eigs = np.linalg.eigvalsh(sigma_true)
    assert np.all(eigs > 0.0)


def test_gmv_closed_form_matches_numerical():
    sigma_true, _ = cov.generate_factor_model(n_assets=20, n_factors=3, seed=42)
    w_cf, var_cf = cov.gmv_closed_form(sigma_true)
    w_num, var_num = cov.gmv_numerical(sigma_true)
    assert np.allclose(w_cf, w_num, atol=1e-8)
    assert abs(var_cf - var_num) < 1e-10
    assert abs(np.sum(w_cf) - 1.0) < 1e-12


def test_ledoit_wolf_properties_and_large_t():
    sigma_true, mu = cov.generate_factor_model(n_assets=10, n_factors=3, seed=42)
    rets = cov.simulate_returns(mu, sigma_true, n_obs=100, n_trials=1, seed=42)[0]
    s_t = cov.sample_covariance(rets, ddof=0)
    sigma_lw, intensity = cov.ledoit_wolf(rets)

    assert 0.0 <= intensity <= 1.0
    assert np.allclose(sigma_lw, sigma_lw.T, atol=1e-12)
    assert np.all(np.linalg.eigvalsh(sigma_lw) > 0.0)
    assert abs(np.trace(sigma_lw) - np.trace(s_t)) < 1e-10

    rets_large = cov.simulate_returns(mu, sigma_true, n_obs=20_000, n_trials=1, seed=99)[0]
    _, intensity_large = cov.ledoit_wolf(rets_large)
    assert 0.0 <= intensity_large < 0.05


def test_wishart_in_sample_ratio_and_reproducibility():
    sigma_true, mu = cov.generate_factor_model(n_assets=20, n_factors=3, seed=42)
    exp1 = cov.run_wishart_lw_experiment(sigma_true, mu, n_obs=100, n_trials=400, batch_size=200, seed=42)
    exp2 = cov.run_wishart_lw_experiment(sigma_true, mu, n_obs=100, n_trials=400, batch_size=200, seed=42)

    assert np.allclose(exp1["in_sample_ratio"], exp2["in_sample_ratio"], atol=1e-15)
    theo_mean = (100 - 20) / (100 - 1)
    sample_mean = float(np.mean(exp1["in_sample_ratio"]))
    se = float(np.std(exp1["in_sample_ratio"], ddof=1) / np.sqrt(400))
    assert abs(sample_mean - theo_mean) < 3.0 * se


def test_marchenko_pastur_eigenvalues_within_bounds():
    n_assets, n_obs = 100, 400
    q = n_assets / n_obs
    lam_min, lam_max = cov.marchenko_pastur_bounds(q)
    eigs = cov.sample_correlation_eigenvalues(n_assets=n_assets, n_obs=n_obs, seed=42)
    tw_tol = 3.0 * (n_assets ** (-2.0 / 3.0))
    assert eigs.min() >= lam_min - tw_tol
    assert eigs.max() <= lam_max + tw_tol
