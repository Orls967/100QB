"""Generator grid tak-seragam (sinh / Tavella-Randall) untuk solver PDE."""

from __future__ import annotations

import numpy as np


def sinh_grid(
    x_min: float,
    x_max: float,
    x_center: float,
    n_points: int,
    concentration: float = 0.15,
) -> np.ndarray:
    """Bangkitkan grid 1D tak-seragam yang terkonsentrasi mulus di sekitar ``x_center``.

    Menggunakan pemetaan ``x(xi) = x_center + c * sinh(c1 * (1 - xi) + c2 * xi)``
    dengan ``xi in [0, 1]`` seragam, sehingga beda hingga 2nd-order tetap berorde 2.
    """
    if n_points < 5:
        raise ValueError("n_points minimal 5")
    c = concentration * (x_max - x_min)
    c1 = np.arcsinh((x_min - x_center) / c)
    c2 = np.arcsinh((x_max - x_center) / c)
    xi = np.linspace(0.0, 1.0, n_points)
    grid = x_center + c * np.sinh(c1 * (1.0 - xi) + c2 * xi)
    grid[0] = x_min
    grid[-1] = x_max
    return grid


def nonuniform_fd_weights(grid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Hitung koefisien beda hingga pusat orde-2 pada grid tak-seragam untuk titik interior ``1..-2``.

    Mengembalikan ``(a1, b1, c1, a2, b2, c2)`` sehingga:
        f'(x_i)  ~= a1 * f_{i-1} + b1 * f_i + c1 * f_{i+1}
        f''(x_i) ~= a2 * f_{i-1} + b2 * f_i + c2 * f_{i+1}
    """
    hm = grid[1:-1] - grid[:-2]
    hp = grid[2:] - grid[1:-1]
    hsum = hm + hp

    a1 = -hp / (hm * hsum)
    b1 = (hp - hm) / (hm * hp)
    c1 = hm / (hp * hsum)

    a2 = 2.0 / (hm * hsum)
    b2 = -2.0 / (hm * hp)
    c2 = 2.0 / (hp * hsum)
    return a1, b1, c1, a2, b2, c2


def cubic_interp_1d(x_grid: np.ndarray, y_vals: np.ndarray, x_target: float) -> float:
    """Interpolasi kubik lokal (Lagrange 4-titik) pada grid 1D terurut."""
    n = len(x_grid)
    idx = int(np.searchsorted(x_grid, x_target))
    i0 = max(0, min(idx - 2, n - 4))
    xs = x_grid[i0 : i0 + 4]
    ys = y_vals[i0 : i0 + 4]
    total = 0.0
    for j in range(4):
        term = ys[j]
        for k in range(4):
            if k != j:
                term *= (x_target - xs[k]) / (xs[j] - xs[k])
        total += term
    return float(total)
