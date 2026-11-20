"""Solver PDE forward Dupire untuk permukaan call ``C(K, T)`` dalam satu sapuan."""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

from .grids import cubic_interp_1d, nonuniform_fd_weights, sinh_grid
from .pde_backward import _solve_tridiag

LocalVolFn = Callable[[np.ndarray, float], np.ndarray]


def _build_dupire_operator(
    K_int: np.ndarray,
    sigma_int: np.ndarray,
    r: float,
    q: float,
    w1: tuple[np.ndarray, np.ndarray, np.ndarray],
    w2: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Operator spasial Dupire: ``0.5 * sigma^2 * K^2 * d2C/dK2 - (r - q) * K * dC/dK - q * C``."""
    a1, b1, c1 = w1
    a2, b2, c2 = w2
    diff = 0.5 * sigma_int * sigma_int * K_int * K_int
    drift = -(r - q) * K_int

    lower = diff * a2 + drift * a1
    diag = diff * b2 + drift * b1 - q
    upper = diff * c2 + drift * c1
    return lower, diag, upper


def solve_forward_dupire_surface(
    S0: float,
    T_targets: Sequence[float],
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    n_space: int = 301,
    n_time_per_year: int = 120,
    K_min: float = 0.0,
    K_max: float = 450.0,
    concentration: float = 0.10,
    rannacher_steps: int = 4,
    time_power: float = 1.0,
) -> dict[str, np.ndarray | dict[float, np.ndarray]]:
    """Selesaikan PDE forward Dupire untuk seluruh permukaan ``C(K, T_j)`` dalam satu sapuan.

    Menggunakan grid tak-seragam ``sinh`` di ``K`` (terkonsentrasi di ``S0``) dan
    Crank-Nicolson dengan 4 langkah awal Rannacher di ``T = 0``.
    """
    T_sorted = sorted(float(t) for t in T_targets)
    T_max = T_sorted[-1]

    K_grid = sinh_grid(K_min, K_max, S0, n_space, concentration=concentration)
    K_int = K_grid[1:-1]

    a1, b1, c1, a2, b2, c2 = nonuniform_fd_weights(K_grid)
    w1 = (a1, b1, c1)
    w2 = (a2, b2, c2)

    # Bangun grid waktu yang mencakup semua T_targets persis
    n_total = max(int(np.ceil(T_max * n_time_per_year)), 8)
    s_grid = np.linspace(0.0, 1.0, n_total + 1)
    base_times = T_max * np.power(s_grid, time_power)
    all_times = np.unique(np.concatenate([base_times, np.asarray(T_sorted, dtype=np.float64)]))
    all_times.sort()

    # Kondisi awal pada T = 0: C(K, 0) = max(S0 - K, 0) dengan cell-averaging di sekitar kink S0
    C = np.maximum(S0 - K_grid, 0.0)
    k_lo = 0.5 * (K_grid[:-2] + K_grid[1:-1])
    k_hi = 0.5 * (K_grid[1:-1] + K_grid[2:])
    mask_kink = (k_lo < S0) & (S0 < k_hi)
    C[1:-1] = np.where(mask_kink, 0.5 * np.square(S0 - k_lo) / (k_hi - k_lo), C[1:-1])

    def bc_low(t_val: float) -> float:
        return float(S0 * np.exp(-q * t_val) - K_grid[0] * np.exp(-r * t_val))

    def bc_high(t_val: float) -> float:
        return 0.0

    snapshots: dict[float, np.ndarray] = {}
    n_ran_full = min(rannacher_steps // 2, len(all_times) - 1)

    t_curr = 0.0
    for idx in range(n_ran_full):
        t_end = float(all_times[idx + 1])
        half_dt = 0.5 * (t_end - t_curr)
        for _ in range(2):
            t_next = t_curr + half_dt
            sig_next = local_vol_fn(K_int, t_next)
            lo, di, up = _build_dupire_operator(K_int, sig_next, r, q, w1, w2)

            c_low_next = bc_low(t_next)
            c_high_next = bc_high(t_next)

            lhs_lo = -half_dt * lo
            lhs_di = 1.0 - half_dt * di
            lhs_up = -half_dt * up

            rhs = C[1:-1].copy()
            rhs[0] -= lhs_lo[0] * c_low_next
            rhs[-1] -= lhs_up[-1] * c_high_next

            C_int = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
            C[0] = c_low_next
            C[1:-1] = C_int
            C[-1] = c_high_next
            t_curr = t_next

        if t_end in T_sorted:
            snapshots[t_end] = C.copy()

    sig_curr = local_vol_fn(K_int, max(t_curr, 1e-6))
    lo_curr, di_curr, up_curr = _build_dupire_operator(K_int, sig_curr, r, q, w1, w2)

    for idx in range(n_ran_full, len(all_times) - 1):
        t_next = float(all_times[idx + 1])
        dt = t_next - t_curr
        half_dt = 0.5 * dt

        sig_next = local_vol_fn(K_int, t_next)
        lo_next, di_next, up_next = _build_dupire_operator(K_int, sig_next, r, q, w1, w2)

        c_low_next = bc_low(t_next)
        c_high_next = bc_high(t_next)

        rhs = C[1:-1] + half_dt * (lo_curr * C[:-2] + di_curr * C[1:-1] + up_curr * C[2:])
        lhs_lo = -half_dt * lo_next
        lhs_di = 1.0 - half_dt * di_next
        lhs_up = -half_dt * up_next

        rhs[0] -= lhs_lo[0] * c_low_next
        rhs[-1] -= lhs_up[-1] * c_high_next

        C_int = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
        C[0] = c_low_next
        C[1:-1] = C_int
        C[-1] = c_high_next

        t_curr = t_next
        lo_curr, di_curr, up_curr = lo_next, di_next, up_next

        for target in T_sorted:
            if abs(t_curr - target) < 1e-13:
                snapshots[target] = C.copy()

    return {
        "K": K_grid,
        "snapshots": snapshots,
    }


def eval_forward_surface_at_strikes(
    surface: dict[str, np.ndarray | dict[float, np.ndarray]],
    T: float,
    strikes: Sequence[float],
) -> np.ndarray:
    """Interpolasi kubik harga call dari hasil sapuan forward Dupire pada ``T`` tertentu."""
    K_grid = surface["K"]
    snapshots = surface["snapshots"]
    # Cari key terdekat (presisi float)
    matched_t = min(snapshots.keys(), key=lambda k: abs(k - T))
    C_vals = snapshots[matched_t]
    return np.asarray([cubic_interp_1d(K_grid, C_vals, float(k)) for k in strikes], dtype=np.float64)
