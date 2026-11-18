"""Helper simulasi: seed handling, random walk, dan Monte Carlo generik."""

from __future__ import annotations

import os
from typing import Callable

import numpy as np

DEFAULT_SEED = 42


def get_rng(seed: int | None | np.random.Generator = None) -> np.random.Generator:
    """Kembalikan ``numpy.random.Generator`` yang reproducible.

    - ``Generator`` → dikembalikan apa adanya.
    - ``int``       → generator dengan seed tersebut.
    - ``None``      → pakai env ``QB100_SEED`` bila ada, kalau tidak ``DEFAULT_SEED``.
    """
    if isinstance(seed, np.random.Generator):
        return seed
    if seed is None:
        env = os.environ.get("QB100_SEED", "").strip()
        seed = int(env) if env else DEFAULT_SEED
    return np.random.default_rng(seed)
