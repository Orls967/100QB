"""Pricing opsi Eropa & Asia dengan Black-Scholes dan simulasi Monte Carlo tervektorisasi."""

from __future__ import annotations

import math

import numpy as np

from .simulation import get_rng


def _norm_cdf(x: float) -> float:
    """CDF distribusi normal standar."""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def bs_call(S0: float, K: float, r: float, sigma: float, T: float) -> float:
    """Harga opsi call Eropa dengan rumus tertutup Black-Scholes."""
    if T <= 0 or sigma <= 0:
        return max(S0 - K * math.exp(-r * T), 0.0)
    vol_sqrt_t = sigma * math.sqrt(T)
    d1 = (math.log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / vol_sqrt_t
    d2 = d1 - vol_sqrt_t
    return S0 * _norm_cdf(d1) - K * math.exp(-r * T) * _norm_cdf(d2)


def bs_put(S0: float, K: float, r: float, sigma: float, T: float) -> float:
    """Harga opsi put Eropa dengan rumus tertutup Black-Scholes."""
    if T <= 0 or sigma <= 0:
        return max(K * math.exp(-r * T) - S0, 0.0)
    vol_sqrt_t = sigma * math.sqrt(T)
    d1 = (math.log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / vol_sqrt_t
    d2 = d1 - vol_sqrt_t
    return K * math.exp(-r * T) * _norm_cdf(-d2) - S0 * _norm_cdf(-d1)


def geometric_asian_call(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
) -> float:
    """Harga tertutup opsi call Asia geometrik dengan monitoring diskret (Kemna-Vorst)."""
    if n_steps < 1:
        raise ValueError("n_steps harus >= 1")
    mu_ln_g = math.log(S0) + (r - 0.5 * sigma * sigma) * T * (n_steps + 1) / (2.0 * n_steps)
    var_ln_g = (sigma * sigma * T) * (n_steps + 1) * (2 * n_steps + 1) / (6.0 * n_steps * n_steps)
    sigma_ln_g = math.sqrt(var_ln_g)
    expected_g = math.exp(mu_ln_g + 0.5 * var_ln_g)
    d1 = (mu_ln_g - math.log(K) + var_ln_g) / sigma_ln_g
    d2 = d1 - sigma_ln_g
    return math.exp(-r * T) * (expected_g * _norm_cdf(d1) - K * _norm_cdf(d2))


def simulate_gbm_paths(
    S0: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator | None = 42,
    antithetic: bool = False,
) -> np.ndarray:
    """Simulasi jalur GBM tervektorisasi; mengembalikan array ``(n_paths, n_steps)`` titik monitoring."""
    if n_steps < 1 or n_paths < 1:
        raise ValueError("n_steps dan n_paths harus >= 1")
    rng = get_rng(seed)
    dt = T / n_steps
    drift = (r - 0.5 * sigma * sigma) * dt
    vol = sigma * math.sqrt(dt)

    if antithetic:
        if n_paths % 2 != 0:
            raise ValueError("n_paths harus genap untuk antithetic variates")
        half = n_paths // 2
        z_half = rng.standard_normal((half, n_steps))
        z = np.vstack([z_half, -z_half])
    else:
        z = rng.standard_normal((n_paths, n_steps))

    log_increments = drift + vol * z
    log_paths = np.cumsum(log_increments, axis=1)
    return S0 * np.exp(log_paths)


def _summarize_estimator(samples: np.ndarray, ci_z: float = 1.96) -> dict[str, float]:
    """Hitung harga, standard error, dan CI 95% dari sampel independen."""
    mean = float(samples.mean())
    se = float(samples.std(ddof=1) / math.sqrt(samples.size))
    return {
        "price": mean,
        "se": se,
        "ci_low": mean - ci_z * se,
        "ci_high": mean + ci_z * se,
    }


def price_european_mc(
    paths: np.ndarray,
    K: float,
    r: float,
    T: float,
) -> dict[str, float]:
    """Estimasi harga call Eropa dari jalur GBM (untuk sanity check)."""
    payoffs = np.maximum(paths[:, -1] - K, 0.0) * math.exp(-r * T)
    return _summarize_estimator(payoffs)


def price_asian_naive(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator | None = 42,
    paths: np.ndarray | None = None,
) -> dict[str, float]:
    """Estimasi naive Monte Carlo untuk arithmetic Asian call."""
    if paths is None:
        paths = simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=seed)
    arith_avg = paths.mean(axis=1)
    payoffs = np.maximum(arith_avg - K, 0.0) * math.exp(-r * T)
    return _summarize_estimator(payoffs)


def price_asian_antithetic(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator | None = 42,
) -> dict[str, float]:
    """Estimasi antithetic variates untuk arithmetic Asian call (n_paths total = 2 * pasangan)."""
    paths = simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=seed, antithetic=True)
    arith_avg = paths.mean(axis=1)
    payoffs = np.maximum(arith_avg - K, 0.0) * math.exp(-r * T)
    half = n_paths // 2
    pair_means = 0.5 * (payoffs[:half] + payoffs[half:])
    return _summarize_estimator(pair_means)


def price_asian_control_variate(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator | None = 42,
    paths: np.ndarray | None = None,
) -> dict[str, float]:
    """Estimasi control variate memakai opsi Asia geometrik diskret."""
    if paths is None:
        paths = simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=seed)
    disc = math.exp(-r * T)
    arith_avg = paths.mean(axis=1)
    geom_avg = np.exp(np.log(paths).mean(axis=1))

    x = np.maximum(arith_avg - K, 0.0) * disc
    y = np.maximum(geom_avg - K, 0.0) * disc
    y_exact = geometric_asian_call(S0, K, r, sigma, T, n_steps)

    y_centered = y - y.mean()
    beta = float(np.dot(x - x.mean(), y_centered) / np.dot(y_centered, y_centered))
    z = x - beta * (y - y_exact)

    res = _summarize_estimator(z)
    res["beta"] = beta
    res["y_exact"] = y_exact
    return res


def delta_asian_pathwise(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    seed: int | np.random.Generator | None = 42,
    paths: np.ndarray | None = None,
) -> dict[str, float]:
    """Delta arithmetic Asian call dengan metode pathwise derivative."""
    if paths is None:
        paths = simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=seed)
    arith_avg = paths.mean(axis=1)
    samples = math.exp(-r * T) * np.where(arith_avg > K, arith_avg / S0, 0.0)
    res = _summarize_estimator(samples)
    return {"delta": res["price"], "se": res["se"], "ci_low": res["ci_low"], "ci_high": res["ci_high"]}


def delta_asian_bump(
    S0: float,
    K: float,
    r: float,
    sigma: float,
    T: float,
    n_steps: int,
    n_paths: int,
    h_rel: float = 0.01,
    seed: int | np.random.Generator | None = 42,
    paths: np.ndarray | None = None,
) -> dict[str, float]:
    """Delta arithmetic Asian call dengan central difference (CRN) dan bump relatif ``h_rel``."""
    if paths is None:
        paths = simulate_gbm_paths(S0, r, sigma, T, n_steps, n_paths, seed=seed)
    h = h_rel * S0
    arith_base = paths.mean(axis=1)
    arith_up = arith_base * ((S0 + h) / S0)
    arith_dn = arith_base * ((S0 - h) / S0)
    disc = math.exp(-r * T)
    samples = disc * (np.maximum(arith_up - K, 0.0) - np.maximum(arith_dn - K, 0.0)) / (2.0 * h)
    res = _summarize_estimator(samples)
    return {"delta": res["price"], "se": res["se"], "ci_low": res["ci_low"], "ci_high": res["ci_high"]}
