"""Volatilitas lokal Dupire dari permukaan SSVI dan verifikasi beda hingga orde 4."""

from __future__ import annotations

import numpy as np

from .bs import bs_call_from_total_variance
from .ssvi import (
    SSVIParams,
    forward_price,
    log_moneyness,
    ssvi_butterfly_g,
    ssvi_dw_dT,
    ssvi_total_variance,
)


def dupire_local_var_ssvi(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
) -> np.ndarray | float:
    """Variansi lokal Dupire analitik ``sigma_loc^2(y, T) = (dw/dT) / g(y, T)``."""
    dw_dT = ssvi_dw_dT(y, T, params)
    g = ssvi_butterfly_g(y, T, params)
    val = dw_dT / g
    return float(val) if np.ndim(val) == 0 else val


def dupire_local_vol_ssvi(
    y: np.ndarray | float,
    T: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
) -> np.ndarray | float:
    """Volatilitas lokal Dupire ``sigma_loc(y, T) = sqrt(dupire_local_var_ssvi)``."""
    var = dupire_local_var_ssvi(y, T, params)
    val = np.sqrt(np.maximum(var, 0.0))
    return float(val) if np.ndim(val) == 0 else val


def dupire_local_vol_from_spot(
    S: np.ndarray | float,
    t: np.ndarray | float,
    params: SSVIParams = SSVIParams(),
    t_min: float = 1e-4,
) -> np.ndarray | float:
    """Pemetaan ``sigma_loc(S, t)`` ke koordinat SSVI: ``K = S``, ``T = max(t, t_min)``, ``y = ln(S / F(T))``."""
    t_eff = np.maximum(np.asarray(t, dtype=np.float64), t_min)
    y = log_moneyness(S, t_eff, S0=params.S0, r=params.r, q=params.q)
    return dupire_local_vol_ssvi(y, t_eff, params)


def dupire_local_var_from_prices_fd4(
    K: float,
    T: float,
    params: SSVIParams = SSVIParams(),
    dK: float = 0.05,
    dT: float = 5e-4,
) -> float:
    """Hitung ``sigma_loc^2(K, T)`` dari harga BS-implisit lewat definisi Dupire (beda hingga orde 4).

    Rumus Dupire:
        dC/dT = 0.5 * sigma_loc^2 * K^2 * d2C/dK2 - (r - q) * K * dC/dK - q * C
    """
    def call_at(k_val: float, t_val: float) -> float:
        y_val = log_moneyness(k_val, t_val, S0=params.S0, r=params.r, q=params.q)
        w_val = ssvi_total_variance(y_val, t_val, params)
        return float(bs_call_from_total_variance(params.S0, k_val, t_val, params.r, params.q, w_val))

    C0 = call_at(K, T)

    # Stencil beda hingga pusat orde 4 untuk turunan pertama & kedua
    c_Tm2 = call_at(K, T - 2.0 * dT)
    c_Tm1 = call_at(K, T - dT)
    c_Tp1 = call_at(K, T + dT)
    c_Tp2 = call_at(K, T + 2.0 * dT)
    dC_dT = (-c_Tp2 + 8.0 * c_Tp1 - 8.0 * c_Tm1 + c_Tm2) / (12.0 * dT)

    c_Km2 = call_at(K - 2.0 * dK, T)
    c_Km1 = call_at(K - dK, T)
    c_Kp1 = call_at(K + dK, T)
    c_Kp2 = call_at(K + 2.0 * dK, T)
    dC_dK = (-c_Kp2 + 8.0 * c_Kp1 - 8.0 * c_Km1 + c_Km2) / (12.0 * dK)
    d2C_dK2 = (-c_Kp2 + 16.0 * c_Kp1 - 30.0 * C0 + 16.0 * c_Km1 - c_Km2) / (12.0 * dK * dK)

    num = dC_dT + (params.r - params.q) * K * dC_dK + params.q * C0
    den = 0.5 * K * K * d2C_dK2
    return float(num / den)
