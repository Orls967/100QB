"""Tes unit untuk Bagian 1 (SSVI) dan Bagian 2 (Volatilitas Lokal Dupire)."""

from __future__ import annotations

import numpy as np
import pytest

from qb100.day003a_local_vol.bs import bs_call_from_total_variance, bs_call_greeks, implied_vol_call
from qb100.day003a_local_vol.dupire import (
    dupire_local_var_from_prices_fd4,
    dupire_local_var_ssvi,
)
from qb100.day003a_local_vol.ssvi import (
    SSVIParams,
    check_butterfly_arbitrage,
    check_calendar_arbitrage,
    check_gatheral_jacquier_sufficient,
    forward_price,
    log_moneyness,
    ssvi_butterfly_g,
    ssvi_d2w_dy2,
    ssvi_dw_dT,
    ssvi_dw_dy,
    ssvi_theta_prime,
    ssvi_total_variance,
)


def test_bs_and_implied_vol_roundtrip():
    g = bs_call_greeks(100.0, 105.0, 1.0, 0.02, 0.01, 0.25)
    iv = implied_vol_call(g["price"], 100.0, 105.0, 1.0, 0.02, 0.01)
    assert abs(iv - 0.25) < 1e-11


def test_ssvi_analytic_derivatives_vs_fd():
    params = SSVIParams()
    y_grid = np.linspace(-1.5, 1.5, 31)
    T_grid = np.linspace(0.05, 2.0, 20)
    Y, TT = np.meshgrid(y_grid, T_grid, indexing="ij")

    hy = 1e-5
    hT = 1e-6

    dw_dy_ana = ssvi_dw_dy(Y, TT, params)
    dw_dy_fd = (
        ssvi_total_variance(Y + hy, TT, params) - ssvi_total_variance(Y - hy, TT, params)
    ) / (2.0 * hy)
    rel_wy = np.max(np.abs(dw_dy_ana - dw_dy_fd) / np.maximum(np.abs(dw_dy_fd), 1e-12))
    assert rel_wy < 1e-6

    hy2 = 1.5e-3
    d2w_dy2_ana = ssvi_d2w_dy2(Y, TT, params)
    d2w_dy2_fd = (
        -ssvi_total_variance(Y + 2.0 * hy2, TT, params)
        + 16.0 * ssvi_total_variance(Y + hy2, TT, params)
        - 30.0 * ssvi_total_variance(Y, TT, params)
        + 16.0 * ssvi_total_variance(Y - hy2, TT, params)
        - ssvi_total_variance(Y - 2.0 * hy2, TT, params)
    ) / (12.0 * hy2 * hy2)
    rel_wyy = np.max(np.abs(d2w_dy2_ana - d2w_dy2_fd) / np.maximum(np.abs(d2w_dy2_fd), 1e-8))
    assert rel_wyy < 1e-6

    dw_dT_ana = ssvi_dw_dT(Y, TT, params)
    dw_dT_fd = (
        ssvi_total_variance(Y, TT + hT, params) - ssvi_total_variance(Y, TT - hT, params)
    ) / (2.0 * hT)
    rel_wT = np.max(np.abs(dw_dT_ana - dw_dT_fd) / np.maximum(np.abs(dw_dT_fd), 1e-12))
    assert rel_wT < 1e-6


def test_ssvi_arbitrage_free_and_counterfactual():
    params = SSVIParams()
    y_grid = np.linspace(-1.5, 1.5, 61)
    T_grid = np.linspace(0.05, 2.0, 40)

    bf = check_butterfly_arbitrage(y_grid, T_grid, params)
    cal = check_calendar_arbitrage(y_grid, T_grid, params)
    gj = check_gatheral_jacquier_sufficient(T_grid, params)

    assert bf["arbitrage_free"] and bf["min_g"] > 0.0
    assert cal["arbitrage_free"] and cal["min_dw_dT"] > 0.0
    assert gj["theorem_4_2_satisfied"]
    assert gj["powerlaw_eta_le_4"]

    # Kontrafaktual: parameter melanggar syarat cukup (eta=4.0, rho=-0.9)
    bad_params = SSVIParams(eta=4.0, rho=-0.9)
    bad_bf = check_butterfly_arbitrage(np.linspace(-2.5, 2.5, 101), T_grid, bad_params)
    assert not bad_bf["arbitrage_free"]
    assert bad_bf["min_g"] < 0.0


def test_call_monotonicity_and_convexity():
    params = SSVIParams()
    K_grid = np.linspace(40.0, 180.0, 71)
    for T in (0.1, 0.5, 1.0, 2.0):
        y = log_moneyness(K_grid, T, S0=params.S0, r=params.r, q=params.q)
        w = ssvi_total_variance(y, T, params)
        C = bs_call_from_total_variance(params.S0, K_grid, T, params.r, params.q, w)
        dC = np.diff(C)
        d2C = np.diff(C, n=2)
        assert np.all(dC <= 1e-12)
        assert np.all(d2C >= -1e-12)


def test_dupire_eta_zero_anchor():
    params_flat = SSVIParams(eta=0.0)
    y_grid = np.linspace(-1.5, 1.5, 21)
    T_grid = np.linspace(0.05, 2.0, 20)
    Y, TT = np.meshgrid(y_grid, T_grid, indexing="ij")
    loc_var = dupire_local_var_ssvi(Y, TT, params_flat)
    expected = ssvi_theta_prime(TT, params_flat)
    assert np.max(np.abs(loc_var - expected)) < 1e-12


def test_dupire_consistency_fd4():
    params = SSVIParams()
    for T in (0.25, 0.5, 1.0, 1.5):
        F = forward_price(T, S0=params.S0, r=params.r, q=params.q)
        for m in (0.8, 0.9, 1.0, 1.1, 1.2):
            K = m * F
            y = log_moneyness(K, T, S0=params.S0, r=params.r, q=params.q)
            var_ana = dupire_local_var_ssvi(y, T, params)
            var_fd4 = dupire_local_var_from_prices_fd4(K, T, params)
            rel_err = abs(var_fd4 - var_ana) / var_ana
            assert rel_err < 1e-4, f"Mismatch at m={m}, T={T}: ana={var_ana}, fd4={var_fd4}"
