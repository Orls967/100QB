"""Unit test untuk fungsi karakteristik, Fourier pricing (Gil-Pelaez & COS), dan implied vol."""

from __future__ import annotations

import numpy as np

from qb100 import heston


def test_martingale_and_large_t_u_stability():
    p = heston.DEFAULT_PARAMS
    for tau in (0.25, 1.0, 5.0):
        phi_neg_i = heston.heston_cf(-1j, tau, 100.0, 0.03, p)
        assert abs(phi_neg_i - 100.0 * np.exp(0.03 * tau)) < 1e-12

    u_large = np.linspace(10.0, 500.0, 200)
    phi_vals = heston.heston_cf(u_large, 5.0, 100.0, 0.03, p)
    assert np.all(np.isfinite(phi_vals))
    assert np.all(np.abs(phi_vals) <= 1.0 + 1e-12)


def test_gil_pelaez_vs_cos_and_parity_and_arbitrage():
    p = heston.DEFAULT_PARAMS
    ks, ts = heston.market_grid(100.0, 0.03)
    calls_gp = heston.heston_gil_pelaez(100.0, ks, 0.03, ts, p, option_type="call")
    calls_cos = heston.heston_cos(100.0, ks, 0.03, ts, p, n_cos=256, option_type="call")
    assert np.max(np.abs(calls_gp - calls_cos)) < 1e-8

    puts_cos = heston.heston_cos(100.0, ks, 0.03, ts, p, n_cos=256, option_type="put")
    parity_diff = calls_cos - puts_cos - (100.0 - ks * np.exp(-0.03 * ts))
    assert np.max(np.abs(parity_diff)) < 1e-8

    lower_bound = np.maximum(100.0 - ks * np.exp(-0.03 * ts), 0.0)
    assert np.all(calls_cos >= lower_bound - 1e-10)
    assert np.all(calls_cos <= 100.0 + 1e-10)
    # Konveksitas terhadap strike per maturitas
    second_diff = np.diff(calls_cos, n=2, axis=1)
    assert np.all(second_diff >= -1e-10)


def test_xi_to_zero_converges_to_black_scholes():
    bs_ref = heston.bs_price(100.0, 100.0, 0.03, 1.0, 0.20, option_type="call")
    errs = []
    for xi_val in (0.1, 0.05, 0.02, 0.01):
        p = heston.HestonParams(v0=0.04, kappa=2.0, theta=0.04, xi=xi_val, rho=-0.7)
        c_val = heston.heston_cos(100.0, 100.0, 0.03, 1.0, p, n_cos=256)
        errs.append(abs(c_val - bs_ref))
    assert all(errs[i] > errs[i + 1] for i in range(len(errs) - 1))
    assert errs[-1] < 2e-3


def test_implied_vol_round_trip():
    ks, ts = heston.market_grid(100.0, 0.03)
    calls = heston.heston_cos(100.0, ks, 0.03, ts, heston.DEFAULT_PARAMS, n_cos=256)
    ivs = heston.implied_vol(calls, 100.0, ks, 0.03, ts)
    reconstructed = heston.bs_price(100.0, ks, 0.03, ts, ivs)
    ivs_rt = heston.implied_vol(reconstructed, 100.0, ks, 0.03, ts)
    assert np.max(np.abs(ivs - ivs_rt)) < 1e-10
