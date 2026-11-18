"""Unit tests untuk modul qb100.options."""

import math

import numpy as np
import pytest

from qb100 import options as opt


def test_bs_put_call_parity():
    S0, K, r, sigma, T = 100.0, 105.0, 0.05, 0.2, 1.0
    c = opt.bs_call(S0, K, r, sigma, T)
    p = opt.bs_put(S0, K, r, sigma, T)
    assert c - p == pytest.approx(S0 - K * math.exp(-r * T), abs=1e-10)


def test_geometric_asian_matches_large_mc():
    S0, K, r, sigma, T, n_steps = 100.0, 100.0, 0.05, 0.2, 1.0, 252
    exact = opt.geometric_asian_call(S0, K, r, sigma, T, n_steps)
    paths = opt.simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths=50_000, seed=42)
    geom_avg = np.exp(np.log(paths).mean(axis=1))
    payoffs = np.maximum(geom_avg - K, 0.0) * math.exp(-r * T)
    mc_mean = payoffs.mean()
    mc_se = payoffs.std(ddof=1) / math.sqrt(len(payoffs))
    assert abs(mc_mean - exact) < 3.0 * mc_se


def test_control_variate_beats_naive_se():
    S0, K, r, sigma, T, n_steps, n_paths = 100.0, 100.0, 0.05, 0.2, 1.0, 252, 20_000
    naive = opt.price_asian_naive(S0, K, r, sigma, T, n_steps, n_paths, seed=42)
    cv = opt.price_asian_control_variate(S0, K, r, sigma, T, n_steps, n_paths, seed=42)
    assert cv["se"] < naive["se"] / 10.0


def test_options_mc_reproducible_same_seed():
    a = opt.price_asian_control_variate(100.0, 100.0, 0.05, 0.2, 1.0, 50, 5_000, seed=7)
    b = opt.price_asian_control_variate(100.0, 100.0, 0.05, 0.2, 1.0, 50, 5_000, seed=7)
    assert a["price"] == b["price"]
    assert a["se"] == b["se"]


def test_delta_pathwise_matches_bump():
    S0, K, r, sigma, T, n_steps, n_paths = 100.0, 100.0, 0.05, 0.2, 1.0, 252, 50_000
    paths = opt.simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=42)
    pw = opt.delta_asian_pathwise(S0, K, r, sigma, T, n_steps, n_paths, paths=paths)
    bmp = opt.delta_asian_bump(S0, K, r, sigma, T, n_steps, n_paths, h_rel=0.01, paths=paths)
    assert abs(pw["delta"] - bmp["delta"]) < 3.0 * math.hypot(pw["se"], bmp["se"])
