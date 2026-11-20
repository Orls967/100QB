"""Harga Black-Scholes, Greeks eksak, dan inversi implied volatility buatan sendiri."""

from __future__ import annotations

import numpy as np
from scipy.special import ndtr

_INV_SQRT_2PI = 1.0 / np.sqrt(2.0 * np.pi)


def _norm_pdf(x: np.ndarray | float) -> np.ndarray | float:
    return _INV_SQRT_2PI * np.exp(-0.5 * np.square(x))


def bs_call_price(
    S: np.ndarray | float,
    K: np.ndarray | float,
    T: np.ndarray | float,
    r: float,
    q: float,
    sigma: np.ndarray | float,
) -> np.ndarray | float:
    """Harga call Eropa Black-Scholes dengan dividend yield kontinu ``q``."""
    S_arr = np.asarray(S, dtype=np.float64)
    K_arr = np.asarray(K, dtype=np.float64)
    T_arr = np.asarray(T, dtype=np.float64)
    sig_arr = np.asarray(sigma, dtype=np.float64)

    scalar = (
        np.ndim(S) == 0
        and np.ndim(K) == 0
        and np.ndim(T) == 0
        and np.ndim(sigma) == 0
    )
    sqrt_T = np.sqrt(np.maximum(T_arr, 0.0))
    vol_sqrt_T = sig_arr * sqrt_T

    with np.errstate(divide="ignore", invalid="ignore"):
        d1 = (np.log(S_arr / K_arr) + (r - q + 0.5 * sig_arr * sig_arr) * T_arr) / vol_sqrt_T
        d2 = d1 - vol_sqrt_T
        price = S_arr * np.exp(-q * T_arr) * ndtr(d1) - K_arr * np.exp(-r * T_arr) * ndtr(d2)

    intrinsic = np.maximum(S_arr * np.exp(-q * T_arr) - K_arr * np.exp(-r * T_arr), 0.0)
    price = np.where((T_arr <= 0.0) | (vol_sqrt_T <= 0.0), intrinsic, price)
    return float(price) if scalar else price


def bs_call_from_total_variance(
    S0: float,
    K: np.ndarray | float,
    T: np.ndarray | float,
    r: float,
    q: float,
    w: np.ndarray | float,
) -> np.ndarray | float:
    """Harga call BS dari total variansi implisit ``w = sigma_imp**2 * T``."""
    T_arr = np.asarray(T, dtype=np.float64)
    w_arr = np.asarray(w, dtype=np.float64)
    sigma = np.sqrt(np.maximum(w_arr, 0.0) / np.maximum(T_arr, 1e-16))
    return bs_call_price(S0, K, T_arr, r, q, sigma)


def bs_vega(
    S: np.ndarray | float,
    K: np.ndarray | float,
    T: np.ndarray | float,
    r: float,
    q: float,
    sigma: np.ndarray | float,
) -> np.ndarray | float:
    """Vega Black-Scholes ``dC/dsigma``."""
    S_arr = np.asarray(S, dtype=np.float64)
    K_arr = np.asarray(K, dtype=np.float64)
    T_arr = np.asarray(T, dtype=np.float64)
    sig_arr = np.asarray(sigma, dtype=np.float64)
    sqrt_T = np.sqrt(np.maximum(T_arr, 1e-16))
    vol_sqrt_T = np.maximum(sig_arr * sqrt_T, 1e-16)
    d1 = (np.log(S_arr / K_arr) + (r - q + 0.5 * sig_arr * sig_arr) * T_arr) / vol_sqrt_T
    vega = S_arr * np.exp(-q * T_arr) * _norm_pdf(d1) * sqrt_T
    return float(vega) if np.ndim(vega) == 0 else vega


def bs_call_greeks(
    S: float,
    K: float,
    T: float,
    r: float,
    q: float,
    sigma: float,
) -> dict[str, float]:
    """Harga dan Greeks eksak Black-Scholes (delta, gamma, vega, theta)."""
    sqrt_T = np.sqrt(T)
    vol_sqrt_T = sigma * sqrt_T
    d1 = (np.log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / vol_sqrt_T
    d2 = d1 - vol_sqrt_T
    disc_q = np.exp(-q * T)
    disc_r = np.exp(-r * T)
    n_d1 = float(ndtr(d1))
    n_d2 = float(ndtr(d2))
    pdf_d1 = float(_norm_pdf(d1))

    price = S * disc_q * n_d1 - K * disc_r * n_d2
    delta = disc_q * n_d1
    gamma = disc_q * pdf_d1 / (S * vol_sqrt_T)
    vega = S * disc_q * pdf_d1 * sqrt_T
    theta = (
        -S * disc_q * pdf_d1 * sigma / (2.0 * sqrt_T)
        + q * S * disc_q * n_d1
        - r * K * disc_r * n_d2
    )
    return {
        "price": float(price),
        "delta": float(delta),
        "gamma": float(gamma),
        "vega": float(vega),
        "theta": float(theta),
    }


def implied_vol_call(
    price: float,
    S: float,
    K: float,
    T: float,
    r: float,
    q: float,
    tol: float = 1e-12,
    max_iter: int = 100,
) -> float:
    """Inversi implied volatility call dengan Newton-Raphson + bisection fallback."""
    disc_q = np.exp(-q * T)
    disc_r = np.exp(-r * T)
    lower_bound = max(S * disc_q - K * disc_r, 0.0)
    upper_bound = S * disc_q

    if price <= lower_bound + 1e-15:
        return 0.0
    if price >= upper_bound - 1e-15:
        raise ValueError(f"Harga {price} melewati batas atas no-arbitrage {upper_bound}")

    lo, hi = 1e-8, 5.0
    while bs_call_price(S, K, T, r, q, hi) < price:
        hi *= 2.0
        if hi > 50.0:
            raise ValueError("Implied volatility melebihi 50.0")

    # Tebakan awal Manaster-Koehler / Brenner-Subrahmanyam
    F = S * np.exp((r - q) * T)
    sigma = float(np.sqrt(2.0 * np.pi / T) * price / F)
    sigma = min(max(sigma, 0.05), 0.5 * hi)

    for _ in range(max_iter):
        c_val = bs_call_price(S, K, T, r, q, sigma)
        diff = c_val - price
        if abs(diff) < tol:
            return float(sigma)
        if diff < 0.0:
            lo = max(lo, sigma)
        else:
            hi = min(hi, sigma)

        vega = bs_vega(S, K, T, r, q, sigma)
        if vega > 1e-14:
            cand = sigma - diff / vega
            if lo < cand < hi:
                sigma = cand
                continue
        sigma = 0.5 * (lo + hi)
        if hi - lo < tol:
            return float(sigma)

    return float(0.5 * (lo + hi))
