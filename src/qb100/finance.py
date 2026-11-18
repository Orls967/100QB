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
