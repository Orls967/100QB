"""Pricing Fourier (Gil-Pelaez & COS Fang-Oosterlee), Black-Scholes, dan Implied Vol."""

from __future__ import annotations

import numpy as np
from scipy import integrate
from scipy.optimize import brentq
from scipy.stats import norm

from .cf import HestonParams, heston_cf, heston_cumulants


def market_grid(
    s0: float = 100.0,
    r: float = 0.03,
    moneyness: tuple[float, ...] = (0.8, 0.9, 1.0, 1.1, 1.2),
    maturities: tuple[float, ...] = (0.25, 0.5, 1.0, 2.0),
) -> tuple[np.ndarray, np.ndarray]:
    """Bangkitkan grid ``(K, T)`` 20 opsi call Eropa dari moneyness ``K / (S0 * exp(rT))``."""
    m_arr = np.asarray(moneyness, dtype=float)
    t_arr = np.asarray(maturities, dtype=float)
    strikes = np.array([s0 * np.exp(r * t) * m_arr for t in t_arr], dtype=float)
    mats = np.repeat(t_arr[:, None], len(m_arr), axis=1)
    return strikes, mats


def bs_price(
    s0: float,
    k: np.ndarray | float,
    r: float,
    t: np.ndarray | float,
    sigma: np.ndarray | float,
    option_type: str = "call",
) -> np.ndarray | float:
    """Harga opsi Eropa Black-Scholes tervektorisasi."""
    k_arr = np.asarray(k, dtype=float)
    t_arr = np.asarray(t, dtype=float)
    sig_arr = np.asarray(sigma, dtype=float)
    sqrt_t = np.sqrt(t_arr)
    vol_t = sig_arr * sqrt_t
    d1 = (np.log(s0 / k_arr) + (r + 0.5 * sig_arr * sig_arr) * t_arr) / vol_t
    d2 = d1 - vol_t
    disc_k = k_arr * np.exp(-r * t_arr)
    if option_type == "call":
        res = s0 * norm.cdf(d1) - disc_k * norm.cdf(d2)
    else:
        res = disc_k * norm.cdf(-d2) - s0 * norm.cdf(-d1)
    return float(res) if np.ndim(k) == 0 and np.ndim(t) == 0 and np.ndim(sigma) == 0 else res


def implied_vol(
    price: np.ndarray | float,
    s0: float,
    k: np.ndarray | float,
    r: float,
    t: np.ndarray | float,
    option_type: str = "call",
) -> np.ndarray | float:
    """Inversi implied volatility Black-Scholes (Newton-Raphson + pengaman Brent)."""
    p_arr = np.atleast_1d(np.asarray(price, dtype=float))
    k_arr = np.broadcast_to(np.asarray(k, dtype=float), p_arr.shape)
    t_arr = np.broadcast_to(np.asarray(t, dtype=float), p_arr.shape)

    iv = np.empty_like(p_arr, dtype=float)
    for idx in np.ndindex(p_arr.shape):
        p_val = float(p_arr[idx])
        k_val = float(k_arr[idx])
        t_val = float(t_arr[idx])
        sig = 0.25
        converged = False
        for _ in range(30):
            model_p = float(bs_price(s0, k_val, r, t_val, sig, option_type=option_type))
            diff = model_p - p_val
            if abs(diff) < 1e-12:
                converged = True
                break
            d1 = (np.log(s0 / k_val) + (r + 0.5 * sig * sig) * t_val) / (sig * np.sqrt(t_val))
            vega = s0 * float(norm.pdf(d1)) * np.sqrt(t_val)
            if vega < 1e-12:
                break
            sig_next = sig - diff / vega
            if sig_next <= 1e-6 or sig_next > 5.0:
                break
            sig = sig_next
        if not converged:
            sig = float(
                brentq(
                    lambda s: float(bs_price(s0, k_val, r, t_val, s, option_type=option_type)) - p_val,
                    1e-7,
                    5.0,
                    xtol=1e-13,
                    rtol=1e-13,
                )
            )
        iv[idx] = sig

    if np.ndim(price) == 0 and np.ndim(k) == 0 and np.ndim(t) == 0:
        return float(iv.flat[0])
    return iv


def _gl_nodes_weights(u_max: float = 250.0, n_sub: int = 8, order: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """Node dan bobot Gauss-Legendre komposit pada ``[0, u_max]``."""
    x_ref, w_ref = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(0.0, u_max, n_sub + 1)
    nodes_list, weights_list = [], []
    for i in range(n_sub):
        a, b = edges[i], edges[i + 1]
        nodes_list.append(0.5 * (b - a) * x_ref + 0.5 * (a + b))
        weights_list.append(0.5 * (b - a) * w_ref)
    return np.concatenate(nodes_list), np.concatenate(weights_list)


def heston_gil_pelaez(
    s0: float,
    k: np.ndarray | float,
    r: float,
    t: np.ndarray | float,
    params: HestonParams,
    option_type: str = "call",
    u_max: float = 600.0,
    use_quad: bool = False,
) -> np.ndarray | float:
    """Harga opsi Eropa Heston via integral Gil-Pelaez (vektorisasi GL atau ``quad``)."""
    k_arr = np.atleast_1d(np.asarray(k, dtype=float))
    t_arr = np.broadcast_to(np.asarray(t, dtype=float), k_arr.shape)
    prices = np.empty_like(k_arr, dtype=float)

    if use_quad:
        for idx in np.ndindex(k_arr.shape):
            kv, tv = float(k_arr[idx]), float(t_arr[idx])
            ln_k = np.log(kv)
            phi_neg_i = s0 * np.exp(r * tv)

            def int1(u: float, lnk: float = ln_k, tau: float = tv, norm_f: float = phi_neg_i) -> float:
                f1 = heston_cf(u - 1j, tau, s0, r, params) / norm_f
                return float(np.real(np.exp(-1j * u * lnk) * f1 / (1j * u)))

            def int2(u: float, lnk: float = ln_k, tau: float = tv) -> float:
                f2 = heston_cf(u, tau, s0, r, params)
                return float(np.real(np.exp(-1j * u * lnk) * f2 / (1j * u)))

            p1 = 0.5 + integrate.quad(int1, 1e-12, u_max, epsabs=1e-11, epsrel=1e-11, limit=200)[0] / np.pi
            p2 = 0.5 + integrate.quad(int2, 1e-12, u_max, epsabs=1e-11, epsrel=1e-11, limit=200)[0] / np.pi
            call_val = s0 * p1 - kv * np.exp(-r * tv) * p2
            prices[idx] = call_val if option_type == "call" else kv * np.exp(-r * tv) * (1.0 - p2) - s0 * (1.0 - p1)
    else:
        u_nodes, w_nodes = _gl_nodes_weights(u_max=u_max, n_sub=16, order=64)
        unique_ts = np.unique(t_arr)
        for tv in unique_ts:
            mask = t_arr == tv
            ks_t = k_arr[mask]
            ln_ks = np.log(ks_t)[:, None]
            f2 = heston_cf(u_nodes, float(tv), s0, r, params)
            f1 = heston_cf(u_nodes - 1j, float(tv), s0, r, params) / (s0 * np.exp(r * float(tv)))
            exp_term = np.exp(-1j * ln_ks * u_nodes[None, :]) / (1j * u_nodes[None, :])
            p1 = 0.5 + np.real(exp_term * f1[None, :]) @ w_nodes / np.pi
            p2 = 0.5 + np.real(exp_term * f2[None, :]) @ w_nodes / np.pi
            if option_type == "call":
                prices[mask] = s0 * p1 - ks_t * np.exp(-r * float(tv)) * p2
            else:
                prices[mask] = ks_t * np.exp(-r * float(tv)) * (1.0 - p2) - s0 * (1.0 - p1)

    if np.ndim(k) == 0 and np.ndim(t) == 0:
        return float(prices.flat[0])
    return prices


def _cos_chi_psi(a: float, b: float, c: float, d: float, k_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Koefisien analitik ``chi_k`` dan ``psi_k`` metode COS Fang-Oosterlee pada ``[c, d]``."""
    omega = k_idx * np.pi / (b - a)
    arg_d = omega * (d - a)
    arg_c = omega * (c - a)
    exp_d = np.exp(d)
    exp_c = np.exp(c)
    chi = (
        np.cos(arg_d) * exp_d - np.cos(arg_c) * exp_c
        + omega * (np.sin(arg_d) * exp_d - np.sin(arg_c) * exp_c)
    ) / (1.0 + omega * omega)

    psi = np.empty_like(omega)
    psi[0] = d - c
    psi[1:] = (np.sin(arg_d[1:]) - np.sin(arg_c[1:])) / omega[1:]
    return chi, psi


def heston_cos(
    s0: float,
    k: np.ndarray | float,
    r: float,
    t: np.ndarray | float,
    params: HestonParams,
    n_cos: int = 256,
    l_trunc: float = 10.0,
    option_type: str = "call",
) -> np.ndarray | float:
    """Harga opsi Eropa Heston dengan metode COS (Fang-Oosterlee 2008), tervektorisasi."""
    k_arr = np.atleast_1d(np.asarray(k, dtype=float))
    t_arr = np.broadcast_to(np.asarray(t, dtype=float), k_arr.shape)
    prices = np.empty_like(k_arr, dtype=float)
    k_idx = np.arange(n_cos, dtype=float)

    for tv in np.unique(t_arr):
        mask = t_arr == tv
        ks_t = k_arr[mask]
        tau = float(tv)
        c1, c2, c4 = heston_cumulants(tau, s0, r, params)
        # Pusat cumulant di ln(S_T / K) agar [a, b] cocok dengan grid strike
        c1_x = c1 - float(np.mean(np.log(ks_t)))
        half_w = l_trunc * np.sqrt(c2 + np.sqrt(c4))
        a = c1_x - half_w
        b = c1_x + half_w

        omega = k_idx * np.pi / (b - a)
        phi_u = heston_cf(omega, tau, s0, r, params)
        x_k = np.log(s0 / ks_t)  # ln(S_0/K); c1_x sudah menggeser phi
        # Gunakan koordinat y = ln(S_T / K): phi_y(u) = phi(u) * exp(-i * u * ln(K))
        exp_shift = np.exp(1j * (-np.log(ks_t)[:, None] - a) * omega[None, :])
        re_term = np.real(phi_u[None, :] * exp_shift)

        if option_type == "call":
            chi, psi = _cos_chi_psi(a, b, 0.0, b, k_idx)
            vk = (2.0 / (b - a)) * (chi - psi)
        else:
            chi, psi = _cos_chi_psi(a, b, a, 0.0, k_idx)
            vk = (2.0 / (b - a)) * (-chi + psi)

        vk[0] *= 0.5
        prices[mask] = np.exp(-r * tau) * ks_t * (re_term @ vk)

    if np.ndim(k) == 0 and np.ndim(t) == 0:
        return float(prices.flat[0])
    return prices
