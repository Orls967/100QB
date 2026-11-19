"""Simulasi Monte Carlo Heston: Skema QE (Andersen 2008) & Euler Full-Truncation."""

from __future__ import annotations

import time

import numpy as np

from .cf import HestonParams
from .pricing import heston_cos


def simulate_qe_paths(
    s0: float,
    r: float,
    t: float,
    params: HestonParams,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator = 42,
    psi_c: float = 1.5,
    chunk_size: int = 50_000,
) -> dict[str, np.ndarray | float]:
    """Simulasi jalur terminal ``S_T`` dengan skema Quadratic-Exponential (Andersen 2008)."""
    v0, kappa, theta, xi, rho = params.as_tuple()
    dt = t / float(n_steps)
    exp_kdt = np.exp(-kappa * dt)
    c_s2_v = (xi * xi * exp_kdt * (1.0 - exp_kdt)) / kappa
    c_s2_0 = (theta * xi * xi * ((1.0 - exp_kdt) ** 2)) / (2.0 * kappa)

    gamma1 = 0.5
    gamma2 = 0.5
    k0 = -rho * kappa * theta * dt / xi
    k1 = gamma1 * dt * (kappa * rho / xi - 0.5) - rho / xi
    k2 = gamma2 * dt * (kappa * rho / xi - 0.5) + rho / xi
    k3 = gamma1 * dt * (1.0 - rho * rho)
    k4 = gamma2 * dt * (1.0 - rho * rho)

    s_term = np.empty(n_paths, dtype=float)
    min_v = float(v0)
    high_psi_count = 0
    total_step_count = n_paths * n_steps

    rng_is_gen = isinstance(seed, np.random.Generator)
    n_chunks = 1 if rng_is_gen else (n_paths + chunk_size - 1) // chunk_size
    child_seeds = [seed] if rng_is_gen else np.random.SeedSequence(seed).spawn(n_chunks)

    offset = 0
    for c_idx, c_seed in enumerate(child_seeds):
        cur = n_paths if rng_is_gen else min(chunk_size, n_paths - offset)
        rng = c_seed if rng_is_gen else np.random.default_rng(c_seed)
        ln_s = np.full(cur, np.log(s0), dtype=float)
        v = np.full(cur, v0, dtype=float)

        for _ in range(n_steps):
            z_v = rng.standard_normal(cur)
            z_s = rng.standard_normal(cur)
            u_v = rng.random(cur)

            m = theta + (v - theta) * exp_kdt
            s2 = v * c_s2_v + c_s2_0
            psi = s2 / (m * m)

            mask_low = psi <= psi_c
            high_psi_count += int(np.sum(~mask_low))
            v_next = np.empty(cur, dtype=float)

            if np.any(mask_low):
                inv_psi = 2.0 / psi[mask_low]
                b2 = inv_psi - 1.0 + np.sqrt(inv_psi) * np.sqrt(np.maximum(inv_psi - 1.0, 0.0))
                a = m[mask_low] / (1.0 + b2)
                v_next[mask_low] = a * ((np.sqrt(b2) + z_v[mask_low]) ** 2)

            if np.any(~mask_low):
                psi_hi = psi[~mask_low]
                m_hi = m[~mask_low]
                u_hi = u_v[~mask_low]
                p_qe = (psi_hi - 1.0) / (psi_hi + 1.0)
                beta_qe = (1.0 - p_qe) / m_hi
                v_hi = np.where(u_hi <= p_qe, 0.0, np.log((1.0 - p_qe) / np.maximum(1.0 - u_hi, 1e-15)) / beta_qe)
                v_next[~mask_low] = v_hi

            min_v = min(min_v, float(np.min(v_next)))
            var_s = np.maximum(k3 * v + k4 * v_next, 0.0)
            ln_s += r * dt + k0 + k1 * v + k2 * v_next + np.sqrt(var_s) * z_s
            v = v_next

        s_term[offset : offset + cur] = np.exp(ln_s)
        offset += cur

    return {"s_t": s_term, "min_v": min_v, "frac_high_psi": high_psi_count / float(total_step_count)}


def simulate_euler_ft_paths(
    s0: float,
    r: float,
    t: float,
    params: HestonParams,
    n_steps: int,
    n_paths: int,
    seed: int = 42,
    chunk_size: int = 50_000,
) -> np.ndarray:
    """Simulasi jalur terminal ``S_T`` dengan skema Euler Full-Truncation."""
    v0, kappa, theta, xi, rho = params.as_tuple()
    dt = t / float(n_steps)
    sq_dt = np.sqrt(dt)
    sq_1m_rho2 = np.sqrt(1.0 - rho * rho)

    s_term = np.empty(n_paths, dtype=float)
    n_chunks = (n_paths + chunk_size - 1) // chunk_size
    child_seeds = np.random.SeedSequence(seed).spawn(n_chunks)

    offset = 0
    for c_seed in child_seeds:
        cur = min(chunk_size, n_paths - offset)
        rng = np.random.default_rng(c_seed)
        ln_s = np.full(cur, np.log(s0), dtype=float)
        v = np.full(cur, v0, dtype=float)
        for _ in range(n_steps):
            z1 = rng.standard_normal(cur)
            z_perp = rng.standard_normal(cur)
            z2 = rho * z1 + sq_1m_rho2 * z_perp
            v_pos = np.maximum(v, 0.0)
            sq_vdt = np.sqrt(v_pos) * sq_dt
            ln_s += (r - 0.5 * v_pos) * dt + sq_vdt * z1
            v += kappa * (theta - v_pos) * dt + xi * sq_vdt * z2
        s_term[offset : offset + cur] = np.exp(ln_s)
        offset += cur
    return s_term


def price_from_paths(
    s_t: np.ndarray,
    k: float,
    r: float,
    t: float,
    s0: float,
    use_cv: bool = False,
) -> dict[str, float]:
    """Hitung harga call Eropa, SE, dan uji martingale dari sampel ``S_T``."""
    disc = np.exp(-r * t)
    payoffs = disc * np.maximum(s_t - k, 0.0)
    n = len(payoffs)
    mean_p = float(np.mean(payoffs))
    var_p = float(np.var(payoffs, ddof=1))
    se_p = float(np.sqrt(var_p / n))

    disc_s = disc * s_t
    mart_mean = float(np.mean(disc_s))
    mart_se = float(np.std(disc_s, ddof=1) / np.sqrt(n))

    out = {
        "price": mean_p,
        "se": se_p,
        "var": var_p,
        "mart_mean": mart_mean,
        "mart_se": mart_se,
        "mart_z": abs(mart_mean - s0) / mart_se,
    }
    if use_cv:
        cov_yx = float(np.cov(payoffs, disc_s, ddof=1)[0, 1])
        var_x = float(np.var(disc_s, ddof=1))
        beta = cov_yx / var_x if var_x > 0.0 else 0.0
        cv_samples = payoffs - beta * (disc_s - s0)
        var_cv = float(np.var(cv_samples, ddof=1))
        out["cv_price"] = float(np.mean(cv_samples))
        out["cv_se"] = float(np.sqrt(var_cv / n))
        out["vrf"] = var_p / var_cv if var_cv > 0.0 else np.nan
        out["beta"] = beta
    return out


def fourier_greeks(
    s0: float,
    k: float,
    r: float,
    t: float,
    params: HestonParams,
    h_s_rel: float = 1e-4,
    h_v: float = 1e-4,
) -> tuple[float, float]:
    """Delta dan ``dPrice/dv0`` via beda hingga sentral pada harga COS."""
    hs = h_s_rel * s0
    c_sp = float(heston_cos(s0 + hs, k, r, t, params, n_cos=256))
    c_sm = float(heston_cos(s0 - hs, k, r, t, params, n_cos=256))
    delta = (c_sp - c_sm) / (2.0 * hs)

    p_vp = HestonParams(params.v0 + h_v, params.kappa, params.theta, params.xi, params.rho)
    p_vm = HestonParams(params.v0 - h_v, params.kappa, params.theta, params.xi, params.rho)
    c_vp = float(heston_cos(s0, k, r, t, p_vp, n_cos=256))
    c_vm = float(heston_cos(s0, k, r, t, p_vm, n_cos=256))
    vega_v0 = (c_vp - c_vm) / (2.0 * h_v)
    return delta, vega_v0


def qe_greeks_crn(
    s0: float,
    k: float,
    r: float,
    t: float,
    params: HestonParams,
    n_steps: int = 64,
    n_paths: int = 200_000,
    seed: int = 42,
    h_s_rel: float = 0.01,
    h_v: float = 1e-3,
) -> dict[str, float]:
    """Delta dan ``dPrice/dv0`` QE Monte Carlo dengan Common Random Numbers (CRN)."""
    t0 = time.perf_counter()
    hs = h_s_rel * s0
    disc = np.exp(-r * t)

    # Delta: karena S0 hanya menggeser ln(S_T) secara aditif, satu simulasi QE sudah eksak CRN!
    s_t_base = simulate_qe_paths(s0, r, t, params, n_steps, n_paths, seed=seed)["s_t"]
    s_t_sp = s_t_base * ((s0 + hs) / s0)
    s_t_sm = s_t_base * ((s0 - hs) / s0)
    delta_samples = disc * (np.maximum(s_t_sp - k, 0.0) - np.maximum(s_t_sm - k, 0.0)) / (2.0 * hs)

    # dPrice/dv0: CRN dengan seed sama pada v0 + h_v dan v0 - h_v
    p_vp = HestonParams(params.v0 + h_v, params.kappa, params.theta, params.xi, params.rho)
    p_vm = HestonParams(params.v0 - h_v, params.kappa, params.theta, params.xi, params.rho)
    s_t_vp = simulate_qe_paths(s0, r, t, p_vp, n_steps, n_paths, seed=seed)["s_t"]
    s_t_vm = simulate_qe_paths(s0, r, t, p_vm, n_steps, n_paths, seed=seed)["s_t"]
    vega_samples = disc * (np.maximum(s_t_vp - k, 0.0) - np.maximum(s_t_vm - k, 0.0)) / (2.0 * h_v)

    return {
        "delta": float(np.mean(delta_samples)),
        "delta_se": float(np.std(delta_samples, ddof=1) / np.sqrt(n_paths)),
        "vega_v0": float(np.mean(vega_samples)),
        "vega_v0_se": float(np.std(vega_samples, ddof=1) / np.sqrt(n_paths)),
        "time_sec": time.perf_counter() - t0,
    }
