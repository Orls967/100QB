"""Tes unit untuk Bagian 4: Greeks (delta, gamma, vega-bucket) dan kalibrasi balik LM."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from qb100.day003a_local_vol.bs import bs_call_greeks
from qb100.day003a_local_vol.calibration import (
    _ssvi_iv_and_jacobian,
    calibrate_ssvi_lm,
    generate_noisy_ssvi_surface,
)
from qb100.day003a_local_vol.greeks import pde_delta_gamma, vega_bucket_analysis
from qb100.day003a_local_vol.ssvi import SSVIParams


def test_pde_delta_gamma_flat_vs_bs():
    S0, K, T, r, q, sigma = 100.0, 100.0, 1.0, 0.02, 0.01, 0.2
    vol_fn = lambda s, t: np.full_like(s, sigma, dtype=np.float64)
    res = pde_delta_gamma(S0, K, T, r, q, vol_fn, n_space=601, n_time=300, fd_bump=0.05)
    bs = bs_call_greeks(S0, K, T, r, q, sigma)

    assert abs(res["delta_grid"] - bs["delta"]) < 1e-5
    assert abs(res["delta_fd"] - bs["delta"]) < 1e-5
    assert abs(res["gamma_grid"] - bs["gamma"]) < 1e-4
    assert abs(res["gamma_fd"] - bs["gamma"]) < 1e-4


def test_vega_bucket_linear_consistency():
    res = vega_bucket_analysis(n_space=301, n_time=200)
    assert res["linear_consistent"], f"rel_diff={res['rel_diff']}"


def test_calibration_jacobian_vs_fd():
    y_vals = np.log(np.array([0.8, 0.9, 1.0, 1.1, 1.2]))
    T_grid = np.array([0.25, 0.5, 1.0, 1.5, 2.0])
    beta = np.array([-0.6, 0.9, 0.45, 0.015, 0.028, 0.05, 0.07, 0.09])

    _, jac_ana = _ssvi_iv_and_jacobian(beta, y_vals, T_grid)
    jac_fd = np.zeros_like(jac_ana)
    h = 1e-6
    for j in range(len(beta)):
        bp = beta.copy()
        bm = beta.copy()
        bp[j] += h
        bm[j] -= h
        iv_p, _ = _ssvi_iv_and_jacobian(bp, y_vals, T_grid)
        iv_m, _ = _ssvi_iv_and_jacobian(bm, y_vals, T_grid)
        jac_fd[:, j] = (iv_p - iv_m) / (2.0 * h)

    assert np.max(np.abs(jac_ana - jac_fd)) < 1e-6


def test_calibration_noiseless_and_noisy_crosscheck():
    params = SSVIParams()
    clean = generate_noisy_ssvi_surface(params=params, noise_std=0.0, seed=20261121)
    fit_clean = calibrate_ssvi_lm(clean["iv_true"], clean["y_vals"], clean["T_grid"])
    assert fit_clean["rms_vol_error"] < 1e-10
    assert abs(fit_clean["rho"] - params.rho) < 1e-8
    assert abs(fit_clean["eta"] - params.eta) < 1e-8
    assert abs(fit_clean["gamma"] - params.gamma) < 1e-8

    noisy = generate_noisy_ssvi_surface(params=params, noise_std=5e-4, seed=20261121)
    fit_noisy = calibrate_ssvi_lm(noisy["iv_noisy"], noisy["y_vals"], noisy["T_grid"])
    assert 2e-4 < fit_noisy["rms_vol_error"] < 1e-3

    # Pembanding silang dengan scipy.optimize.least_squares
    obs_vec = noisy["iv_noisy"].T.reshape(-1)
    x0 = np.concatenate([[-0.3, 0.6, 0.3], np.full(len(noisy["T_grid"]), 0.04)])
    scipy_fit = least_squares(
        lambda b: _ssvi_iv_and_jacobian(b, noisy["y_vals"], noisy["T_grid"])[0] - obs_vec,
        x0,
        bounds=([-0.995, 1e-5, 0.01] + [1e-6] * 5, [0.995, 10.0, 0.99] + [5.0] * 5),
    )
    assert np.max(np.abs(fit_noisy["beta"] - scipy_fit.x)) < 1e-6
