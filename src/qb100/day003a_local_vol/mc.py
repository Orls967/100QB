"""Monte Carlo Euler (log-space, antitetik, PCG64) dan uji weak order di bawah sigma_loc."""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

LocalVolFn = Callable[[np.ndarray, float], np.ndarray]

# Toleransi weak order yang dinyatakan SEBELUM menjalankan eksperimen:
# Pada N = 400.000 jalur antitetik dan {25, 50, 100, 200} langkah/tahun dengan CRN,
# kemiringan log-log galat diskretisasi diharapkan berada di [-1.35, -0.65] (target -1.0 +/- 0.35).
WEAK_ORDER_SLOPE_TARGET = -1.0
WEAK_ORDER_SLOPE_TOL = 0.35


def mc_euler_european_calls(
    S0: float,
    strikes: Sequence[float],
    T: float,
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    n_paths: int = 400_000,
    steps_per_year: int = 200,
    seed: int = 20261121,
    batch_pairs: int = 25_000,
) -> dict[str, np.ndarray]:
    """Hitung harga dan standard error call Eropa via log-Euler antitetik.

    Menggunakan ``numpy.random.Generator(PCG64(seed))``.
    """
    if n_paths % 2 != 0:
        raise ValueError("n_paths harus genap untuk antitetik")
    n_pairs = n_paths // 2
    n_steps = max(int(round(T * steps_per_year)), 1)
    dt = T / n_steps
    sqrt_dt = np.sqrt(dt)
    disc = np.exp(-r * T)
    K_arr = np.asarray(strikes, dtype=np.float64)

    rng = np.random.Generator(np.random.PCG64(seed))
    pair_payoffs = np.empty((n_pairs, len(K_arr)), dtype=np.float64)

    done = 0
    while done < n_pairs:
        m = min(batch_pairs, n_pairs - done)
        x_p = np.full(m, np.log(S0), dtype=np.float64)
        x_m = np.full(m, np.log(S0), dtype=np.float64)

        for step in range(n_steps):
            t = max(step * dt, 1e-6)
            z = rng.standard_normal(m)
            s_p = np.exp(x_p)
            s_m = np.exp(x_m)
            sig_p = local_vol_fn(s_p, t)
            sig_m = local_vol_fn(s_m, t)

            x_p += (r - q - 0.5 * sig_p * sig_p) * dt + sig_p * sqrt_dt * z
            x_m += (r - q - 0.5 * sig_m * sig_m) * dt - sig_m * sqrt_dt * z

        sT_p = np.exp(x_p)[:, None]
        sT_m = np.exp(x_m)[:, None]
        pay_p = np.maximum(sT_p - K_arr[None, :], 0.0)
        pay_m = np.maximum(sT_m - K_arr[None, :], 0.0)
        pair_payoffs[done : done + m, :] = 0.5 * disc * (pay_p + pay_m)
        done += m

    prices = np.mean(pair_payoffs, axis=0)
    se = np.std(pair_payoffs, axis=0, ddof=1) / np.sqrt(n_pairs)
    return {"strikes": K_arr, "prices": prices, "se": se}


def mc_euler_weak_convergence(
    S0: float,
    K: float,
    T: float,
    r: float,
    q: float,
    local_vol_fn: LocalVolFn,
    ref_price: float,
    steps_list: Sequence[int] = (25, 50, 100, 200),
    n_paths: int = 400_000,
    seed: int = 20261121,
    batch_pairs: int = 25_000,
) -> dict[str, object]:
    """Uji weak convergence order Euler pada ``steps_list`` menggunakan Common Random Numbers."""
    max_steps = max(steps_list)
    for s in steps_list:
        if max_steps % s != 0:
            raise ValueError("Semua steps_list harus membagi max(steps_list) untuk CRN")

    n_pairs = n_paths // 2
    disc = np.exp(-r * T)
    dt_fine = T / max_steps
    sqrt_dt_fine = np.sqrt(dt_fine)

    rng = np.random.Generator(np.random.PCG64(seed))
    pair_payoffs = np.empty((len(steps_list), n_pairs), dtype=np.float64)

    done = 0
    while done < n_pairs:
        m = min(batch_pairs, n_pairs - done)
        # Bangkitkan inkremen Brownian halus dW ~ N(0, dt_fine)
        dW_fine = rng.standard_normal((m, max_steps)) * sqrt_dt_fine

        for idx_s, n_steps in enumerate(steps_list):
            stride = max_steps // n_steps
            dt = T / n_steps
            # Jumlahkan blok inkremen Brownian untuk koarsening
            dW = dW_fine.reshape(m, n_steps, stride).sum(axis=2)

            x_p = np.full(m, np.log(S0), dtype=np.float64)
            x_m = np.full(m, np.log(S0), dtype=np.float64)
            for step in range(n_steps):
                t = max(step * dt, 1e-6)
                sig_p = local_vol_fn(np.exp(x_p), t)
                sig_m = local_vol_fn(np.exp(x_m), t)
                x_p += (r - q - 0.5 * sig_p * sig_p) * dt + sig_p * dW[:, step]
                x_m += (r - q - 0.5 * sig_m * sig_m) * dt - sig_m * dW[:, step]

            pay_p = np.maximum(np.exp(x_p) - K, 0.0)
            pay_m = np.maximum(np.exp(x_m) - K, 0.0)
            pair_payoffs[idx_s, done : done + m] = 0.5 * disc * (pay_p + pay_m)

        done += m

    prices = np.mean(pair_payoffs, axis=1)
    se = np.std(pair_payoffs, axis=1, ddof=1) / np.sqrt(n_pairs)
    abs_errors = np.abs(prices - ref_price)
    slope = float(np.polyfit(np.log(steps_list), np.log(np.maximum(abs_errors, 1e-12)), 1)[0])

    return {
        "steps": list(steps_list),
        "prices": prices,
        "se": se,
        "abs_errors": abs_errors,
        "slope": slope,
        "target_slope": WEAK_ORDER_SLOPE_TARGET,
        "tolerance": WEAK_ORDER_SLOPE_TOL,
        "within_tolerance": bool(abs(slope - WEAK_ORDER_SLOPE_TARGET) <= WEAK_ORDER_SLOPE_TOL),
    }
