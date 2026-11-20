"""Permukaan volatilitas SSVI (Gatheral-Jacquier), turunan analitik, dan uji bebas-arbitrase."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SSVIParams:
    """Parameter tetap permukaan SSVI dan kurva pasar."""

    S0: float = 100.0
    r: float = 0.02
    q: float = 0.01
    rho: float = -0.7
    eta: float = 1.0
    gamma: float = 0.5
    theta_slope: float = 0.04
    theta_amp: float = 0.01
    theta_rate: float = 2.0


def forward_price(T: np.ndarray | float, S0: float = 100.0, r: float = 0.02, q: float = 0.01) -> np.ndarray | float:
    """Harga forward ``F(T) = S0 * exp((r - q) * T)``."""
    return S0 * np.exp((r - q) * np.asarray(T, dtype=np.float64))


def log_moneyness(
    K: np.ndarray | float,
    T: np.ndarray | float,
    S0: float = 100.0,
    r: float = 0.02,
    q: float = 0.01,
) -> np.ndarray | float:
    """Log-moneyness forward ``y = ln(K / F(T))``."""
    F = forward_price(T, S0=S0, r=r, q=q)
    return np.log(np.asarray(K, dtype=np.float64) / F)


def ssvi_theta(T: np.ndarray | float, params: SSVIParams = SSVIParams()) -> np.ndarray | float:
    """Variansi ATM total ``theta(T) = 0.04*T + 0.01*(1 - exp(-2*T))``."""
    T_arr = np.asarray(T, dtype=np.float64)
    val = params.theta_slope * T_arr + params.theta_amp * (1.0 - np.exp(-params.theta_rate * T_arr))
    return float(val) if np.ndim(T) == 0 else val


def ssvi_theta_prime(T: np.ndarray | float, params: SSVIParams = SSVIParams()) -> np.ndarray | float:
    """Turunan pertama ``theta'(T)``."""
    T_arr = np.asarray(T, dtype=np.float64)
    val = params.theta_slope + params.theta_amp * params.theta_rate * np.exp(-params.theta_rate * T_arr)
    return float(val) if np.ndim(T) == 0 else val


def ssvi_phi(theta: np.ndarray | float, params: SSVIParams = SSVIParams()) -> np.ndarray | float:
    """Fungsi power-law ``phi(theta) = eta / (theta^gamma * (1 + theta)^(1 - gamma))``."""
    th = np.asarray(theta, dtype=np.float64)
    val = params.eta / (np.power(th, params.gamma) * np.power(1.0 + th, 1.0 - params.gamma))
    return float(val) if np.ndim(theta) == 0 else val


def ssvi_phi_prime(theta: np.ndarray | float, params: SSVIParams = SSVIParams()) -> np.ndarray | float:
    """Turunan analitik ``phi'(theta) = -phi(theta) * (gamma + theta) / (theta * (1 + theta))``."""
    th = np.asarray(theta, dtype=np.float64)
    phi = ssvi_phi(th, params)
    val = -phi * (params.gamma + th) / (th * (1.0 + th))
    return float(val) if np.ndim(theta) == 0 else val


def ssvi_total_variance(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
    theta_override: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """Total variansi implisit SSVI ``w(y, theta(T))``."""
    y_arr = np.asarray(y, dtype=np.float64)
    th = np.asarray(ssvi_theta(T, params) if theta_override is None else theta_override, dtype=np.float64)
    phi = ssvi_phi(th, params)
    z = phi * y_arr + params.rho
    R = np.sqrt(z * z + (1.0 - params.rho * params.rho))
    w = 0.5 * th * (1.0 + params.rho * phi * y_arr + R)
    return float(w) if np.ndim(w) == 0 else w


def ssvi_dw_dy(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
    theta_override: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """Turunan analitik pertama ``dw/dy``."""
    y_arr = np.asarray(y, dtype=np.float64)
    th = np.asarray(ssvi_theta(T, params) if theta_override is None else theta_override, dtype=np.float64)
    phi = ssvi_phi(th, params)
    z = phi * y_arr + params.rho
    R = np.sqrt(z * z + (1.0 - params.rho * params.rho))
    val = 0.5 * th * phi * (params.rho + z / R)
    return float(val) if np.ndim(val) == 0 else val


def ssvi_d2w_dy2(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
    theta_override: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """Turunan analitik kedua ``d^2w/dy^2``."""
    y_arr = np.asarray(y, dtype=np.float64)
    th = np.asarray(ssvi_theta(T, params) if theta_override is None else theta_override, dtype=np.float64)
    phi = ssvi_phi(th, params)
    z = phi * y_arr + params.rho
    R = np.sqrt(z * z + (1.0 - params.rho * params.rho))
    val = 0.5 * th * phi * phi * (1.0 - params.rho * params.rho) / (R * R * R)
    return float(val) if np.ndim(val) == 0 else val


def ssvi_dw_dtheta(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
    theta_override: np.ndarray | float | None = None,
) -> np.ndarray | float:
    """Turunan analitik ``dw/dtheta`` pada ``y`` tetap."""
    y_arr = np.asarray(y, dtype=np.float64)
    th = np.asarray(ssvi_theta(T, params) if theta_override is None else theta_override, dtype=np.float64)
    phi = ssvi_phi(th, params)
    phi_p = ssvi_phi_prime(th, params)
    z = phi * y_arr + params.rho
    R = np.sqrt(z * z + (1.0 - params.rho * params.rho))
    A = 1.0 + params.rho * phi * y_arr + R
    val = 0.5 * A + 0.5 * th * phi_p * y_arr * (params.rho + z / R)
    return float(val) if np.ndim(val) == 0 else val


def ssvi_dw_dT(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
) -> np.ndarray | float:
    """Turunan analitik ``dw/dT = (dw/dtheta) * theta'(T)`` pada ``y`` tetap."""
    dw_dth = ssvi_dw_dtheta(y, T, params)
    th_p = ssvi_theta_prime(T, params)
    val = dw_dth * th_p
    return float(val) if np.ndim(val) == 0 else val


def ssvi_butterfly_g(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
) -> np.ndarray | float:
    """Fungsi densitas kupu-kupu Gatheral ``g(y, T)``."""
    w = ssvi_total_variance(y, T, params)
    wy = ssvi_dw_dy(y, T, params)
    wyy = ssvi_d2w_dy2(y, T, params)
    y_arr = np.asarray(y, dtype=np.float64)
    term1 = np.square(1.0 - (y_arr * wy) / (2.0 * w))
    term2 = 0.25 * np.square(wy) * (1.0 / w + 0.25)
    g = term1 - term2 + 0.5 * wyy
    return float(g) if np.ndim(g) == 0 else g


def check_butterfly_arbitrage(
    y_grid: np.ndarray,
    T_grid: np.ndarray,
    params: SSVIParams = SSVIParams(),
) -> dict[str, float | bool]:
    """Cek bebas-arbitrase kupu-kupu ``g(y, T) >= 0`` pada grid ``y_grid x T_grid``."""
    Y, TT = np.meshgrid(y_grid, T_grid, indexing="ij")
    G = ssvi_butterfly_g(Y, TT, params)
    min_g = float(np.min(G))
    idx = np.unravel_index(int(np.argmin(G)), G.shape)
    return {
        "min_g": min_g,
        "min_y": float(y_grid[idx[0]]),
        "min_T": float(T_grid[idx[1]]),
        "arbitrage_free": bool(min_g >= 0.0),
    }


def check_calendar_arbitrage(
    y_grid: np.ndarray,
    T_grid: np.ndarray,
    params: SSVIParams = SSVIParams(),
) -> dict[str, float | bool]:
    """Cek bebas-arbitrase kalender ``dw/dT >= 0`` pada grid ``y_grid x T_grid``."""
    Y, TT = np.meshgrid(y_grid, T_grid, indexing="ij")
    dw_dT = ssvi_dw_dT(Y, TT, params)
    min_dw_dT = float(np.min(dw_dT))
    idx = np.unravel_index(int(np.argmin(dw_dT)), dw_dT.shape)
    return {
        "min_dw_dT": min_dw_dT,
        "min_y": float(y_grid[idx[0]]),
        "min_T": float(T_grid[idx[1]]),
        "arbitrage_free": bool(min_dw_dT >= 0.0),
    }


def check_gatheral_jacquier_sufficient(
    T_grid: np.ndarray,
    params: SSVIParams = SSVIParams(),
) -> dict[str, float | bool]:
    """Syarat cukup Gatheral-Jacquier (2014) Theorem 4.1 & 4.2 + power-law Corollary 4.1."""
    th = ssvi_theta(T_grid, params)
    phi = ssvi_phi(th, params)
    abs_rho = abs(params.rho)
    cond_lin = float(np.max(th * phi * (1.0 + abs_rho)))
    cond_quad = float(np.max(th * phi * phi * (1.0 + abs_rho)))
    eta_lin = float(params.eta * (1.0 + abs_rho))
    eta_quad = float(params.eta * params.eta * (1.0 + abs_rho))
    return {
        "max_theta_phi_1_plus_abs_rho": cond_lin,
        "max_theta_phi2_1_plus_abs_rho": cond_quad,
        "eta_1_plus_abs_rho": eta_lin,
        "eta2_1_plus_abs_rho": eta_quad,
        "theorem_4_2_satisfied": bool(cond_lin <= 4.0 and cond_quad <= 4.0),
        "powerlaw_eta_le_4": bool(eta_lin <= 4.0),
        "powerlaw_strict_le_2": bool(eta_lin <= 2.0 and eta_quad <= 4.0),
    }
