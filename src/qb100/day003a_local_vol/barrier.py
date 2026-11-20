"""Opsi barrier Down-and-Out / Down-and-In: Reiner-Rubinstein, PDE backward, dan MC jembatan Brownian."""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.special import ndtr

from .bs import bs_call_price
from .grids import cubic_interp_1d, nonuniform_fd_weights, sinh_grid
from .pde_backward import _build_spatial_operator, _solve_tridiag

LocalVolFn = Callable[[np.ndarray, float], np.ndarray]


def reiner_rubinstein_down_call(
    S0: float,
    K: float,
    H: float,
    T: float,
    r: float,
    q: float,
    sigma: float,
) -> dict[str, float]:
    """Rumus tertutup Reiner-Rubinstein (1991) untuk Down-and-Out & Down-and-In Call (S0 > H)."""
    if H <= 0.0:
        vanilla = float(bs_call_price(S0, K, T, r, q, sigma))
        return {"down_out": vanilla, "down_in": 0.0, "vanilla": vanilla}
    if S0 <= H:
        vanilla = float(bs_call_price(S0, K, T, r, q, sigma))
        return {"down_out": 0.0, "down_in": vanilla, "vanilla": vanilla}

    sqrt_T = np.sqrt(T)
    vol_sqrt_T = sigma * sqrt_T
    mu = (r - q - 0.5 * sigma * sigma) / (sigma * sigma)
    disc_q = np.exp(-q * T)
    disc_r = np.exp(-r * T)

    x1 = np.log(S0 / K) / vol_sqrt_T + (1.0 + mu) * vol_sqrt_T
    x2 = np.log(S0 / H) / vol_sqrt_T + (1.0 + mu) * vol_sqrt_T
    y1 = np.log((H * H) / (S0 * K)) / vol_sqrt_T + (1.0 + mu) * vol_sqrt_T
    y2 = np.log(H / S0) / vol_sqrt_T + (1.0 + mu) * vol_sqrt_T

    h_ratio = H / S0
    pow1 = h_ratio ** (2.0 * (mu + 1.0))
    pow2 = h_ratio ** (2.0 * mu)

    A = S0 * disc_q * float(ndtr(x1)) - K * disc_r * float(ndtr(x1 - vol_sqrt_T))
    B = S0 * disc_q * float(ndtr(x2)) - K * disc_r * float(ndtr(x2 - vol_sqrt_T))
    C = S0 * disc_q * pow1 * float(ndtr(y1)) - K * disc_r * pow2 * float(ndtr(y1 - vol_sqrt_T))
    D = S0 * disc_q * pow1 * float(ndtr(y2)) - K * disc_r * pow2 * float(ndtr(y2 - vol_sqrt_T))

    if K >= H:
        down_in = C
        down_out = A - C
    else:
        down_in = A - B + D
        down_out = B - D

    return {
        "down_out": float(down_out),
        "down_in": float(down_in),
        "vanilla": float(A),
    }


def solve_backward_down_and_out_call(
    S0: float,
    K: float,
    H: float,
    T: float,
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    n_space: int = 401,
    n_time: int = 200,
    x_width_stds: float = 6.0,
    concentration: float = 0.12,
    rannacher_steps: int = 4,
    time_power: float = 1.0,
) -> float:
    """Hitung harga Down-and-Out Call dengan PDE backward (batas absorbing di ``x = ln H``)."""
    if S0 <= H:
        return 0.0

    ref_vol = 0.25
    x0 = float(np.log(S0))
    xh = float(np.log(H))
    xk = float(np.log(K))
    x_max = max(x0, xk) + max(x_width_stds * ref_vol * np.sqrt(T), 1.5)

    x_center = xk if xh < xk < x_max else 0.5 * (xh + x_max)
    x = sinh_grid(xh, x_max, x_center, n_space, concentration=concentration)
    S_grid = np.exp(x)
    S_int = S_grid[1:-1]

    a1, b1, c1, a2, b2, c2 = nonuniform_fd_weights(x)
    w1 = (a1, b1, c1)
    w2 = (a2, b2, c2)

    V = np.maximum(S_grid - K, 0.0)
    V[0] = 0.0

    def bc_low(tau: float) -> float:
        return 0.0

    def bc_high(tau: float) -> float:
        return float(S_grid[-1] * np.exp(-q * tau) - K * np.exp(-r * tau))

    s_grid = np.linspace(1.0, 0.0, n_time + 1)
    tau_grid = T * (1.0 - np.power(s_grid, time_power))
    n_ran_full = min(rannacher_steps // 2, n_time)

    tau = 0.0
    for step_idx in range(n_ran_full):
        dt = float(tau_grid[step_idx + 1] - tau_grid[step_idx])
        half_dt = 0.5 * dt
        for _ in range(2):
            tau_next = tau + half_dt
            t_next = max(T - tau_next, 1e-6)
            sig_next = local_vol_fn(S_int, t_next)
            lo, di, up = _build_spatial_operator(x[1:-1], sig_next, r, q, w1, w2)

            v_low_next = bc_low(tau_next)
            v_high_next = bc_high(tau_next)

            lhs_lo = -half_dt * lo
            lhs_di = 1.0 - half_dt * di
            lhs_up = -half_dt * up

            rhs = V[1:-1].copy()
            rhs[0] -= lhs_lo[0] * v_low_next
            rhs[-1] -= lhs_up[-1] * v_high_next

            V[1:-1] = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
            V[0] = v_low_next
            V[-1] = v_high_next
            tau = tau_next

    t_curr = max(T - tau, 1e-6)
    sig_curr = local_vol_fn(S_int, t_curr)
    lo_curr, di_curr, up_curr = _build_spatial_operator(x[1:-1], sig_curr, r, q, w1, w2)

    for step_idx in range(n_ran_full, n_time):
        tau_next = float(tau_grid[step_idx + 1])
        dt = tau_next - tau
        t_next = max(T - tau_next, 1e-6)
        sig_next = local_vol_fn(S_int, t_next)
        lo_next, di_next, up_next = _build_spatial_operator(x[1:-1], sig_next, r, q, w1, w2)

        v_low_next = bc_low(tau_next)
        v_high_next = bc_high(tau_next)

        half_dt = 0.5 * dt
        rhs = V[1:-1] + half_dt * (lo_curr * V[:-2] + di_curr * V[1:-1] + up_curr * V[2:])
        lhs_lo = -half_dt * lo_next
        lhs_di = 1.0 - half_dt * di_next
        lhs_up = -half_dt * up_next

        rhs[0] -= lhs_lo[0] * v_low_next
        rhs[-1] -= lhs_up[-1] * v_high_next

        V[1:-1] = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
        V[0] = v_low_next
        V[-1] = v_high_next

        tau = tau_next
        lo_curr, di_curr, up_curr = lo_next, di_next, up_next

    return cubic_interp_1d(x, V, x0)


def mc_down_and_out_brownian_bridge(
    S0: float,
    K: float,
    H: float,
    T: float,
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    n_paths: int = 400_000,
    steps_per_year: int = 200,
    seed: int = 20261121,
    batch_pairs: int = 25_000,
) -> dict[str, float]:
    """Monte Carlo Down-and-Out Call dengan koreksi jembatan Brownian (sigma lokal beku per langkah)."""
    if S0 <= H:
        return {"price": 0.0, "se": 0.0}
    n_pairs = n_paths // 2
    n_steps = max(int(round(T * steps_per_year)), 1)
    dt = T / n_steps
    sqrt_dt = np.sqrt(dt)
    disc = np.exp(-r * T)
    xh = float(np.log(H))

    rng = np.random.Generator(np.random.PCG64(seed))
    pair_payoffs = np.empty(n_pairs, dtype=np.float64)

    done = 0
    while done < n_pairs:
        m = min(batch_pairs, n_pairs - done)
        x_p = np.full(m, np.log(S0), dtype=np.float64)
        x_m = np.full(m, np.log(S0), dtype=np.float64)
        surv_p = np.ones(m, dtype=np.float64)
        surv_m = np.ones(m, dtype=np.float64)

        for step in range(n_steps):
            t = max(step * dt, 1e-6)
            z = rng.standard_normal(m)
            sig_p = local_vol_fn(np.exp(x_p), t)
            sig_m = local_vol_fn(np.exp(x_m), t)

            x_p_next = x_p + (r - q - 0.5 * sig_p * sig_p) * dt + sig_p * sqrt_dt * z
            x_m_next = x_m + (r - q - 0.5 * sig_m * sig_m) * dt - sig_m * sqrt_dt * z

            alive_p = (x_p > xh) & (x_p_next > xh)
            alive_m = (x_m > xh) & (x_m_next > xh)

            p_bridge_p = np.where(
                alive_p,
                1.0 - np.exp(-2.0 * (x_p - xh) * (x_p_next - xh) / (sig_p * sig_p * dt)),
                0.0,
            )
            p_bridge_m = np.where(
                alive_m,
                1.0 - np.exp(-2.0 * (x_m - xh) * (x_m_next - xh) / (sig_m * sig_m * dt)),
                0.0,
            )
            surv_p *= p_bridge_p
            surv_m *= p_bridge_m
            x_p = x_p_next
            x_m = x_m_next

        pay_p = surv_p * np.maximum(np.exp(x_p) - K, 0.0)
        pay_m = surv_m * np.maximum(np.exp(x_m) - K, 0.0)
        pair_payoffs[done : done + m] = 0.5 * disc * (pay_p + pay_m)
        done += m

    return {
        "price": float(np.mean(pair_payoffs)),
        "se": float(np.std(pair_payoffs, ddof=1) / np.sqrt(n_pairs)),
    }
