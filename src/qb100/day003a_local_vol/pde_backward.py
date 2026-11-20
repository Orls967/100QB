"""Solver PDE backward Black-Scholes lokal di koordinat ``x = ln(S)`` (Crank-Nicolson + Rannacher)."""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.linalg import solve_banded

from .grids import cubic_interp_1d, nonuniform_fd_weights, sinh_grid


LocalVolFn = Callable[[np.ndarray, float], np.ndarray]


def _build_spatial_operator(
    x_int: np.ndarray,
    sigma_int: np.ndarray,
    r: float,
    q: float,
    w1: tuple[np.ndarray, np.ndarray, np.ndarray],
    w2: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bangun tridiagonal ``(lower, diag, upper)`` untuk operator spasial ``L`` pada titik interior."""
    a1, b1, c1 = w1
    a2, b2, c2 = w2
    var = sigma_int * sigma_int
    diff = 0.5 * var
    drift = r - q - diff

    lower = diff * a2 + drift * a1
    diag = diff * b2 + drift * b1 - r
    upper = diff * c2 + drift * c1
    return lower, diag, upper


def _solve_tridiag(
    lower: np.ndarray,
    diag: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    """Selesaikan sistem tridiagonal dengan ``scipy.linalg.solve_banded``."""
    n = len(diag)
    ab = np.zeros((3, n), dtype=np.float64)
    ab[0, 1:] = upper[:-1]
    ab[1, :] = diag
    ab[2, :-1] = lower[1:]
    return solve_banded((1, 1), ab, rhs, overwrite_ab=True, overwrite_b=True, check_finite=False)


def solve_backward_european_call(
    S0: float,
    K: float,
    T: float,
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    n_space: int = 201,
    n_time: int = 100,
    x_width_stds: float = 6.0,
    concentration: float = 0.12,
    rannacher_steps: int = 4,
    time_power: float = 1.0,
    return_grid: bool = False,
    x_bounds: tuple[float, float] | None = None,
) -> float | dict[str, np.ndarray | float]:
    """Hitung harga call Eropa dengan PDE backward di ``x = ln(S)``.

    Menggunakan grid tak-seragam ``sinh`` di sekitar ``ln(K)`` dan ``ln(S0)``,
    4 langkah awal Rannacher (implisit penuh ``dt/2``) lalu Crank-Nicolson.
    """
    if T <= 0.0:
        price = max(S0 - K, 0.0)
        return {"price": price} if return_grid else price

    ref_vol = 0.25
    x0 = float(np.log(S0))
    xk = float(np.log(K))
    half_span = max(x_width_stds * ref_vol * np.sqrt(T), 1.5) + abs(x0 - xk)
    x_min = xk - half_span
    x_max = xk + half_span
    if x_bounds is not None:  # grid tetap, tidak bergantung S0 (untuk beda hingga bersih)
        x_min, x_max = x_bounds

    x = sinh_grid(x_min, x_max, xk, n_space, concentration=concentration)
    S_grid = np.exp(x)
    S_int = S_grid[1:-1]

    a1, b1, c1, a2, b2, c2 = nonuniform_fd_weights(x)
    w1 = (a1, b1, c1)
    w2 = (a2, b2, c2)

    # Kondisi terminal pada tau = 0 (t = T): payoff call
    V = np.maximum(S_grid - K, 0.0)

    def bc_low(tau: float) -> float:
        return 0.0

    def bc_high(tau: float) -> float:
        return float(S_grid[-1] * np.exp(-q * tau) - K * np.exp(-r * tau))

    # Grid waktu tau in [0, T]: jika time_power > 1, rapatkan di dekat t = 0 (tau = T)
    s_grid = np.linspace(1.0, 0.0, n_time + 1)
    tau_grid = T * (1.0 - np.power(s_grid, time_power))
    n_ran_full = min(rannacher_steps // 2, n_time)

    # Langkah 1..n_ran_full: masing-masing 2 half-step implisit penuh (Rannacher)
    tau = 0.0
    for step_idx in range(n_ran_full):
        dt = float(tau_grid[step_idx + 1] - tau_grid[step_idx])
        half_dt = 0.5 * dt
        for _ in range(2):
            tau_next = tau + half_dt
            t_next = max(T - tau_next, 1e-6)
            sig_next = local_vol_fn(S_int, t_next)
            lo, di, up = _build_spatial_operator(x[1:-1], sig_next, r, q, w1, w2)

            v_low_next = bc_low(tau_next)
            v_high_next = bc_high(tau_next)

            lhs_lo = -half_dt * lo
            lhs_di = 1.0 - half_dt * di
            lhs_up = -half_dt * up

            rhs = V[1:-1].copy()
            rhs[0] -= lhs_lo[0] * v_low_next
            rhs[-1] -= lhs_up[-1] * v_high_next

            V_int = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
            V[0] = v_low_next
            V[1:-1] = V_int
            V[-1] = v_high_next
            tau = tau_next

    # Langkah sisanya: Crank-Nicolson
    t_curr = max(T - tau, 1e-6)
    sig_curr = local_vol_fn(S_int, t_curr)
    lo_curr, di_curr, up_curr = _build_spatial_operator(x[1:-1], sig_curr, r, q, w1, w2)

    for step_idx in range(n_ran_full, n_time):
        tau_next = float(tau_grid[step_idx + 1])
        dt = tau_next - tau
        t_next = max(T - tau_next, 1e-6)
        sig_next = local_vol_fn(S_int, t_next)
        lo_next, di_next, up_next = _build_spatial_operator(x[1:-1], sig_next, r, q, w1, w2)

        v_low_next = bc_low(tau_next)
        v_high_next = bc_high(tau_next)

        half_dt = 0.5 * dt
        # RHS: (I + 0.5*dt*L_curr) * V_curr
        rhs = (
            V[1:-1]
            + half_dt * (lo_curr * V[:-2] + di_curr * V[1:-1] + up_curr * V[2:])
        )
        lhs_lo = -half_dt * lo_next
        lhs_di = 1.0 - half_dt * di_next
        lhs_up = -half_dt * up_next

        rhs[0] -= lhs_lo[0] * v_low_next
        rhs[-1] -= lhs_up[-1] * v_high_next

        V_int = _solve_tridiag(lhs_lo, lhs_di, lhs_up, rhs)
        V[0] = v_low_next
        V[1:-1] = V_int
        V[-1] = v_high_next

        tau = tau_next
        lo_curr, di_curr, up_curr = lo_next, di_next, up_next

    price = cubic_interp_1d(x, V, x0)
    if return_grid:
        return {
            "price": price,
            "x": x,
            "S": S_grid,
            "V": V,
        }
    return price
