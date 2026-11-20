"""Tes unit untuk Bagian 3: PDE backward & PDE forward Dupire."""

from __future__ import annotations

import numpy as np

from qb100.day003a_local_vol.bs import bs_call_from_total_variance, bs_call_price, implied_vol_call
from qb100.day003a_local_vol.dupire import dupire_local_vol_from_spot
from qb100.day003a_local_vol.pde_backward import solve_backward_european_call
from qb100.day003a_local_vol.pde_forward import (
    eval_forward_surface_at_strikes,
    solve_forward_dupire_surface,
)
from qb100.day003a_local_vol.ssvi import SSVIParams, forward_price, log_moneyness, ssvi_total_variance


def _flat_vol_fn(sigma: float):
    return lambda s, t: np.full_like(s, sigma, dtype=np.float64)


def test_flat_convergence_order_backward_and_forward():
    S0, K, T, r, q, sigma = 100.0, 105.0, 1.0, 0.02, 0.01, 0.2
    exact = bs_call_price(S0, K, T, r, q, sigma)
    vol_fn = _flat_vol_fn(sigma)

    resolutions = [50, 100, 200, 400]
    err_bwd = []
    err_fwd = []
    for n in resolutions:
        p_bwd = solve_backward_european_call(
            S0, K, T, r, q, vol_fn, n_space=n + 1, n_time=n
        )
        err_bwd.append(abs(p_bwd - exact))

        surf = solve_forward_dupire_surface(
            S0, [T], r, q, vol_fn, n_space=n + 1, n_time_per_year=n
        )
        p_fwd = float(eval_forward_surface_at_strikes(surf, T, [K])[0])
        err_fwd.append(abs(p_fwd - exact))

    log_n = np.log(resolutions)
    slope_bwd = float(np.polyfit(log_n, np.log(err_bwd), 1)[0])
    slope_fwd = float(np.polyfit(log_n, np.log(err_fwd), 1)[0])

    assert abs(slope_bwd - (-2.0)) < 0.15, f"Backward slope={slope_bwd}, errs={err_bwd}"
    assert abs(slope_fwd - (-2.0)) < 0.15, f"Forward slope={slope_fwd}, errs={err_fwd}"


def test_ssvi_pde_vs_bs_implied_and_inter_solver():
    params = SSVIParams()
    vol_fn = lambda s, t: dupire_local_vol_from_spot(s, t, params)

    T_vals = [0.25, 0.5, 1.0, 2.0]
    m_vals = [0.8, 0.9, 1.0, 1.1, 1.2]

    surf = solve_forward_dupire_surface(
        params.S0, T_vals, params.r, params.q, vol_fn, n_space=501, n_time_per_year=400
    )

    max_iv_err_bwd = 0.0
    max_iv_err_fwd = 0.0
    max_price_diff = 0.0

    for T in T_vals:
        F = forward_price(T, S0=params.S0, r=params.r, q=params.q)
        strikes = [m * F for m in m_vals]
        p_fwd_list = eval_forward_surface_at_strikes(surf, T, strikes)

        for m, K, p_fwd in zip(m_vals, strikes, p_fwd_list):
            y = log_moneyness(K, T, S0=params.S0, r=params.r, q=params.q)
            w = ssvi_total_variance(y, T, params)
            iv_true = float(np.sqrt(w / T))
            p_bwd = solve_backward_european_call(
                params.S0,
                K,
                T,
                params.r,
                params.q,
                vol_fn,
                n_space=501,
                n_time=max(int(400 * T), 200),
                time_power=2.0,
            )
            iv_bwd = implied_vol_call(p_bwd, params.S0, K, T, params.r, params.q)
            iv_fwd = implied_vol_call(float(p_fwd), params.S0, K, T, params.r, params.q)

            max_iv_err_bwd = max(max_iv_err_bwd, abs(iv_bwd - iv_true))
            max_iv_err_fwd = max(max_iv_err_fwd, abs(iv_fwd - iv_true))
            max_price_diff = max(max_price_diff, abs(p_bwd - float(p_fwd)))

    assert max_iv_err_bwd < 2e-4, f"max_iv_err_bwd={max_iv_err_bwd}"
    assert max_iv_err_fwd < 2e-4, f"max_iv_err_fwd={max_iv_err_fwd}"
    assert max_price_diff < 1e-3, f"max_price_diff={max_price_diff}"
