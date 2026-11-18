"""Fungsi keuangan dasar: returns, log returns, volatilitas, Sharpe, max drawdown.

Semua fungsi menerima ``pandas.Series``/``DataFrame`` atau array NumPy.
Konvensi: ``periods_per_year=252`` untuk data harian (hari bursa).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

TRADING_DAYS = 252


def _to_pandas(x):
    """Ubah array/list menjadi Series supaya operasi seragam."""
    if isinstance(x, (pd.Series, pd.DataFrame)):
        return x
    arr = np.asarray(x, dtype=float)
    return pd.DataFrame(arr) if arr.ndim == 2 else pd.Series(arr)


def simple_returns(prices):
    """Return sederhana: ``r_t = P_t / P_{t-1} - 1`` (baris pertama dibuang)."""
    p = _to_pandas(prices)
    return p.pct_change().iloc[1:]


def log_returns(prices):
    """Log return: ``r_t = ln(P_t / P_{t-1})`` (baris pertama dibuang)."""
    p = _to_pandas(prices)
    return np.log(p / p.shift(1)).iloc[1:]


def volatility(returns, periods_per_year: int = TRADING_DAYS, annualize: bool = True):
    """Volatilitas = standar deviasi sampel return (ddof=1).

    Bila ``annualize=True`` hasil dikali ``sqrt(periods_per_year)``.
    """
    r = _to_pandas(returns)
    vol = r.std(ddof=1)
    return vol * np.sqrt(periods_per_year) if annualize else vol


def sharpe_ratio(returns, risk_free: float = 0.0, periods_per_year: int = TRADING_DAYS):
    """Sharpe ratio tahunan.

    ``risk_free`` adalah suku bunga bebas risiko TAHUNAN; dikonversi ke per-periode
    lalu ``Sharpe = mean(excess) / std(excess) * sqrt(periods_per_year)``.
    Mengembalikan ``nan`` bila standar deviasi nol.
    """
    r = _to_pandas(returns)
    excess = r - risk_free / periods_per_year
    std = excess.std(ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        sr = excess.mean() / std * np.sqrt(periods_per_year)
    flat = np.isclose(std, 0.0, atol=1e-14)
    if np.isscalar(sr):
        return float("nan") if (flat or not np.isfinite(sr)) else float(sr)
    return sr.where(~flat)


def drawdown(prices):
    """Seri drawdown: ``P_t / max(P_0..P_t) - 1`` (nilai <= 0)."""
    p = _to_pandas(prices)
    return p / p.cummax() - 1.0


def max_drawdown(prices):
    """Max drawdown: penurunan terdalam dari puncak sebelumnya (nilai negatif atau 0)."""
    dd = drawdown(prices).min()
    return float(dd) if np.isscalar(dd) else dd
