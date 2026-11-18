"""Helper ambil & cache data harga dari yfinance ke ``data/cache``.

Data di-cache sebagai CSV supaya notebook tetap bisa jalan offline setelah
pengambilan pertama dan tidak membebani API berulang kali.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

from . import PROJECT_ROOT

CACHE_DIR = PROJECT_ROOT / os.environ.get("QB100_CACHE_DIR", "data/cache")


def cache_path(
    ticker: str,
    start: str | None = None,
    end: str | None = None,
    interval: str = "1d",
    cache_dir: str | Path | None = None,
) -> Path:
    """Bangun path file cache yang deterministik untuk kombinasi parameter."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", ticker.upper())
    parts = [safe, start or "min", end or "max", interval]
    folder = Path(cache_dir) if cache_dir is not None else CACHE_DIR
    return folder / ("_".join(parts) + ".csv")


def get_prices(
    ticker: str,
    start: str | None = None,
    end: str | None = None,
    interval: str = "1d",
    refresh: bool = False,
    cache_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Ambil data harga OHLCV satu ticker, pakai cache CSV bila tersedia.

    Parameters
    ----------
    ticker : simbol yfinance, mis. ``"BBCA.JK"`` atau ``"AAPL"``.
    start, end : tanggal ``YYYY-MM-DD`` (opsional).
    interval : interval data yfinance (default harian ``"1d"``).
    refresh : ``True`` untuk memaksa unduh ulang walau cache ada.
    cache_dir : folder cache alternatif (default ``data/cache``).

    Returns
    -------
    DataFrame ber-index tanggal dengan kolom Open/High/Low/Close/Volume.
    """
    path = cache_path(ticker, start, end, interval, cache_dir)
    if path.exists() and not refresh:
        return pd.read_csv(path, index_col=0, parse_dates=True)

    import yfinance as yf  # import lokal: hanya dibutuhkan saat unduh

    df = yf.download(
        ticker,
        start=start,
        end=end,
        interval=interval,
        auto_adjust=True,
        progress=False,
    )
    if df is None or df.empty:
        raise ValueError(f"Tidak ada data untuk ticker {ticker!r}.")
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = df.columns.get_level_values(0)

    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path)
    return df


def get_close(
    tickers: str | list[str],
    start: str | None = None,
    end: str | None = None,
    refresh: bool = False,
) -> pd.DataFrame:
    """Ambil harga penutupan (Close) satu atau beberapa ticker dalam satu DataFrame."""
    if isinstance(tickers, str):
        tickers = [tickers]
    closes = {t: get_prices(t, start, end, refresh=refresh)["Close"] for t in tickers}
    return pd.DataFrame(closes).dropna(how="all")
