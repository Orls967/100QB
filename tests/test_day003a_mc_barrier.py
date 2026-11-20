"""Tes unit untuk Monte Carlo Euler dan Opsi Barrier (Reiner-Rubinstein, PDE, Brownian Bridge)."""

from __future__ import annotations

import numpy as np

from qb100.day003a_local_vol.barrier import (
    mc_down_and_out_brownian_bridge,
    reiner_rubinstein_down_call,
    solve_backward_down_and_out_call,
)
from qb100.day003a_local_vol.bs import bs_call_from_total_variance, bs_call_price
from qb100.day003a_local_vol.dupire import dupire_local_vol_from_spot
from qb100.day003a_local_vol.mc import mc_euler_european_calls, mc_euler_weak_convergence
from qb100.day003a_local_vol.ssvi import SSVIParams, forward_price, log_moneyness, ssvi_total_variance


def test_reiner_rubinstein_parity_and_limit():
    S0, K, H, T, r, q, sigma = 100.0, 100.0, 85.0, 1.0, 0.02, 0.01, 0.2
    rr = reiner_rubinstein_down_call(S0, K, H, T, r, q, sigma)
    vanilla = bs_call_price(S0, K, T, r, q, sigma)
    assert abs((rr["down_out"] + rr["down_in"]) - vanilla) < 1e-12

    # Uji juga kasus K < H
    rr_kh = reiner_rubinstein_down_call(S0, 75.0, 85.0, T, r, q, sigma)
    vanilla_kh = bs_call_price(S0, 75.0, T, r, q, sigma)
    assert abs((rr_kh["down_out"] + rr_kh["down_in"]) - vanilla_kh) < 1e-12

    # Limit H -> 0 menuju vanila
    rr_tiny = reiner_rubinstein_down_call(S0, K, 1e-6, T, r, q, sigma)
    assert abs(rr_tiny["down_out"] - vanilla) < 1e-10


def test_barrier_pde_vs_reiner_rubinstein_flat():
    S0, K, H, T, r, q, sigma = 100.0, 100.0, 85.0, 1.0, 0.02, 0.01, 0.2
    rr = reiner_rubinstein_down_call(S0, K, H, T, r, q, sigma)["down_out"]
    vol_fn = lambda s, t: np.full_like(s, sigma, dtype=np.float64)
    pde = solve_backward_down_and_out_call(
        S0, K, H, T, r, q, vol_fn, n_space=1601, n_time=800
    )
    assert abs(pde - rr) < 1e-5, f"pde={pde}, rr={rr}, diff={abs(pde - rr)}"


def test_barrier_ssvi_pde_vs_mc_bridge():
    params = SSVIParams()
    vol_fn = lambda s, t: dupire_local_vol_from_spot(s, t, params)
    H, K, T = 85.0, 100.0, 1.0
    pde = solve_backward_down_and_out_call(
        params.S0, K, H, T, params.r, params.q, vol_fn, n_space=501, n_time=300, time_power=2.0
    )
    mc = mc_down_and_out_brownian_bridge(
        params.S0, K, H, T, params.r, params.q, vol_fn, n_paths=100_000, steps_per_year=150, seed=20261121
    )
    z = abs(mc["price"] - pde) / mc["se"]
    assert z < 3.0, f"pde={pde}, mc={mc['price']} +/- {mc['se']}, z={z}"


def test_mc_euler_ssvi_within_3se():
    params = SSVIParams()
    vol_fn = lambda s, t: dupire_local_vol_from_spot(s, t, params)
    T = 1.0
    F = forward_price(T, S0=params.S0, r=params.r, q=params.q)
    strikes = [m * F for m in (0.8, 0.9, 1.0, 1.1, 1.2)]
    exact = [
        bs_call_from_total_variance(
            params.S0,
            k,
            T,
            params.r,
            params.q,
            ssvi_total_variance(log_moneyness(k, T, params.S0, params.r, params.q), T, params),
        )
        for k in strikes
    ]
    mc = mc_euler_european_calls(
        params.S0, strikes, T, params.r, params.q, vol_fn, n_paths=80_000, steps_per_year=150, seed=20261121
    )
    z_scores = np.abs(mc["prices"] - np.asarray(exact)) / mc["se"]
    assert np.all(z_scores < 3.0), f"z_scores={z_scores}"
