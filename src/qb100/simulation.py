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


def random_walk(
    n_steps: int,
    n_paths: int = 1,
    step_std: float = 1.0,
    drift: float = 0.0,
    start: float = 0.0,
    seed: int | None | np.random.Generator = None,
) -> np.ndarray:
    """Random walk Gaussian.

    Returns
    -------
    Array bentuk ``(n_steps + 1, n_paths)``; baris pertama = ``start``.
    """
    if n_steps < 0 or n_paths < 1:
        raise ValueError("n_steps harus >= 0 dan n_paths harus >= 1")
    rng = get_rng(seed)
    steps = rng.normal(loc=drift, scale=step_std, size=(n_steps, n_paths))
    paths = np.empty((n_steps + 1, n_paths))
    paths[0] = start
    paths[1:] = start + np.cumsum(steps, axis=0)
    return paths
