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
