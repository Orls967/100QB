"""qb100 — helper package untuk project 100QB (100 Days of Quant + Kebun Simulasi).

Modul:
- ``plotting``   : style grafik seragam dan fungsi simpan figure.
- ``data``       : ambil & cache data harga (yfinance) ke ``data/cache``.
- ``finance``    : returns, log returns, volatilitas, Sharpe, max drawdown.
- ``simulation`` : random walk, Monte Carlo generik, seed handling.
"""

from pathlib import Path

__version__ = "0.1.0"

#: Root repo (dua level di atas ``src/qb100``).
PROJECT_ROOT = Path(__file__).resolve().parents[2]

__all__ = ["PROJECT_ROOT", "__version__"]
