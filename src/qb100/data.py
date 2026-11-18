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
