"""Unit test untuk kalibrasi Heston dan Monte Carlo QE Andersen."""

from __future__ import annotations

import numpy as np

from qb100 import heston


def test_noise_free_calibration_recovers_params():
    ks, ts = heston.market_grid(100.0, 0.03)
    calls = heston.heston_cos(100.0, ks, 0.03, ts, heston.DEFAULT_PARAMS, n_cos=256)
    market_iv = heston.implied_vol(calls, 100.0, ks, 0.03, ts)

    res = heston.calibrate_heston(market_iv, 100.0, ks, 0.03, ts, n_cos=256)
    assert res["rmse"] < 1e-7
    p_rec = res["params"]
    for est, true_val in zip(p_rec.as_tuple(), heston.DEFAULT_PARAMS.as_tuple()):
        assert abs(est - true_val) / abs(true_val) < 1e-3


def test_qe_monte_carlo_martingale_and_price_and_seed():
    ref = float(heston.heston_cos(100.0, 100.0, 0.03, 1.0, heston.DEFAULT_PARAMS, n_cos=256))
    sim1 = heston.simulate_qe_paths(100.0, 0.03, 1.0, heston.DEFAULT_PARAMS, n_steps=32, n_paths=20_000, seed=42)
    sim2 = heston.simulate_qe_paths(100.0, 0.03, 1.0, heston.DEFAULT_PARAMS, n_steps=32, n_paths=20_000, seed=42)

    assert np.array_equal(sim1["s_t"], sim2["s_t"])
    assert sim1["min_v"] >= 0.0
    assert 0.0 < sim1["frac_high_psi"] < 1.0

    stats = heston.price_from_paths(sim1["s_t"], 100.0, 0.03, 1.0, 100.0, use_cv=True)
    assert stats["mart_z"] < 3.0
    assert abs(stats["cv_price"] - ref) < 3.0 * stats["cv_se"]
    assert stats["vrf"] > 1.5
