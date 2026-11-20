"""Greeks dari PDE backward (delta, gamma) dan sensitivitas vega-bucket terhadap theta(T_i)."""

from __future__ import annotations

import numpy as np

from .bs import bs_call_greeks
from .dupire import dupire_local_vol_from_spot
from .grids import cubic_interp_1d, nonuniform_fd_weights
from .pde_backward import solve_backward_european_call
from .ssvi import (
    SSVIParams,
    log_moneyness,
    ssvi_butterfly_g,
    ssvi_dw_dtheta,
    ssvi_theta,
    ssvi_total_variance,
)


def pde_delta_gamma(
    S0: float,
    K: float,
    T: float,
    r: float,
    q: float,
    local_vol_fn,
    n_space: int = 401,
    n_time: int = 200,
    fd_bump: float = 0.1,
    time_power: float = 1.0,
) -> dict[str, float]:
    """Hitung harga, delta, dan gamma langsung dari grid PDE serta lewat beda hingga S0 +/- h."""
    res = solve_backward_european_call(
        S0, K, T, r, q, local_vol_fn,
        n_space=n_space, n_time=n_time, time_power=time_power, return_grid=True,
    )
    x = res["x"]
    V = res["V"]
    price = float(res["price"])

    a1, b1, c1, a2, b2, c2 = nonuniform_fd_weights(x)
    dV_dx_int = a1 * V[:-2] + b1 * V[1:-1] + c1 * V[2:]
    d2V_dx2_int = a2 * V[:-2] + b2 * V[1:-1] + c2 * V[2:]

    x0 = float(np.log(S0))
    dV_dx = cubic_interp_1d(x[1:-1], dV_dx_int, x0)
    d2V_dx2 = cubic_interp_1d(x[1:-1], d2V_dx2_int, x0)

    delta_grid = dV_dx / S0
    gamma_grid = (d2V_dx2 - dV_dx) / (S0 * S0)

    # Grid TETAP (batas x ditentukan dari S0 dasar) agar bump S0 tidak menggeser grid:
    # beda hingga harga hanya membaca respons PDE yang sama, bukan derau grid yang berubah.
    xk = float(np.log(K))
    half_span = max(6.0 * 0.25 * np.sqrt(T), 1.5) + abs(x0 - xk)
    x_bounds = (xk - half_span, xk + half_span)
    p_up = float(solve_backward_european_call(
        S0 + fd_bump, K, T, r, q, local_vol_fn,
        n_space=n_space, n_time=n_time, time_power=time_power, x_bounds=x_bounds,
    ))
    p_dn = float(solve_backward_european_call(
        S0 - fd_bump, K, T, r, q, local_vol_fn,
        n_space=n_space, n_time=n_time, time_power=time_power, x_bounds=x_bounds,
    ))
    price_fixed = float(solve_backward_european_call(
        S0, K, T, r, q, local_vol_fn,
        n_space=n_space, n_time=n_time, time_power=time_power, x_bounds=x_bounds,
    ))
    delta_fd = (p_up - p_dn) / (2.0 * fd_bump)
    gamma_fd = (p_up - 2.0 * price_fixed + p_dn) / (fd_bump * fd_bump)

    return {
        "price": price,
        "delta_grid": float(delta_grid),
        "gamma_grid": float(gamma_grid),
        "delta_fd": float(delta_fd),
        "gamma_fd": float(gamma_fd),
    }


def _piecewise_theta_and_prime(
    t: float,
    buckets: np.ndarray,
    thetas: np.ndarray,
) -> tuple[float, float]:
    """Interpolasi linear sepotong-sepotong theta(t) melalui (0, 0) dan (buckets[i], thetas[i])."""
    ts = np.concatenate([[0.0], buckets])
    vs = np.concatenate([[0.0], thetas])
    idx = int(np.searchsorted(ts, t, side="right")) - 1
    idx = max(0, min(idx, len(buckets) - 1))
    dt = ts[idx + 1] - ts[idx]
    slope = (vs[idx + 1] - vs[idx]) / dt
    val = vs[idx] + slope * (t - ts[idx])
    return float(max(val, 1e-10)), float(max(slope, 1e-10))


def _local_vol_from_theta_buckets(
    S: np.ndarray,
    t: float,
    params: SSVIParams,
    buckets: np.ndarray,
    thetas: np.ndarray,
) -> np.ndarray:
    t_eff = max(float(t), 1e-4)
    th, th_p = _piecewise_theta_and_prime(t_eff, buckets, thetas)
    y = log_moneyness(S, t_eff, S0=params.S0, r=params.r, q=params.q)
    dw_dth = ssvi_dw_dtheta(y, t_eff, params, theta_override=th)
    # Hitung g dengan theta_override=th
    from .ssvi import ssvi_d2w_dy2, ssvi_dw_dy, ssvi_total_variance

    w = ssvi_total_variance(y, t_eff, params, theta_override=th)
    wy = ssvi_dw_dy(y, t_eff, params, theta_override=th)
    wyy = ssvi_d2w_dy2(y, t_eff, params, theta_override=th)
    g = np.square(1.0 - (y * wy) / (2.0 * w)) - 0.25 * np.square(wy) * (1.0 / w + 0.25) + 0.5 * wyy
    loc_var = np.maximum((dw_dth * th_p) / g, 1e-10)
    return np.sqrt(loc_var)


def vega_bucket_analysis(
    params: SSVIParams = SSVIParams(),
    T: float = 1.0,
    buckets: tuple[float, ...] = (0.25, 0.5, 0.75, 1.0),
    bump: float = 1e-4,
    n_space: int = 401,
    n_time: int = 240,
) -> dict[str, object]:
    """Hitung sensitivitas harga call ATM T=1 terhadap 4 bucket theta(T_i) vs pergeseran paralel."""
    b_arr = np.asarray(buckets, dtype=np.float64)
    th_base = np.asarray(ssvi_theta(b_arr, params), dtype=np.float64)
    K_atm = params.S0  # ATM spot (atau bisa forward; prompt menyebut "harga ATM T = 1")

    def price_with_thetas(th_vec: np.ndarray) -> float:
        vol_fn = lambda s, t: _local_vol_from_theta_buckets(s, t, params, b_arr, th_vec)
        return float(
            solve_backward_european_call(
                params.S0, K_atm, T, params.r, params.q, vol_fn,
                n_space=n_space, n_time=n_time, time_power=2.0,
            )
        )

    base_price = price_with_thetas(th_base)
    bucket_sens = []
    for i in range(len(b_arr)):
        th_bump = th_base.copy()
        th_bump[i] += bump
        p_bump = price_with_thetas(th_bump)
        bucket_sens.append((p_bump - base_price) / bump)

    sum_buckets = float(np.sum(bucket_sens))
    p_parallel = price_with_thetas(th_base + bump)
    parallel_sens = float((p_parallel - base_price) / bump)
    rel_diff = abs(sum_buckets - parallel_sens) / max(abs(parallel_sens), 1e-12)

    return {
        "buckets": [float(x) for x in b_arr],
        "base_price": base_price,
        "bucket_sensitivities": [float(x) for x in bucket_sens],
        "sum_bucket_sensitivities": sum_buckets,
        "parallel_sens": parallel_sens,
        "rel_diff": float(rel_diff),
        "linear_consistent": bool(rel_diff < 1e-3),
    }


def gamma_diagnosis_ssvi(params: SSVIParams = SSVIParams()) -> dict[str, object]:
    """Diagnosis gamma SSVI di S0: grid vs beda hingga, konvergensi, titik non-sudut, t_min, dan C_KK.

    Temuan yang diukur (bukan disetel): di S0 = params.S0 permukaan sigma_loc(S, t) singular pada
    t -> 0 (lebar lapisan batas ~ sqrt(theta(t))), sehingga gamma di S0 bergantung pada regularisasi
    ``t_min`` dan tidak konvergen; di titik lain grid dan beda hingga cocok.
    """
    from .bs import bs_call_from_total_variance
    from .pde_forward import solve_forward_dupire_surface

    r, q = params.r, params.q
    vol_fn = lambda s, t: dupire_local_vol_from_spot(s, t, params)

    conv = []
    for n_space, n_time, h in [(601, 400, 0.05), (601, 400, 0.5), (1201, 800, 0.5), (2401, 1600, 0.5), (2401, 1600, 1.0)]:
        res = pde_delta_gamma(params.S0, params.S0, 1.0, r, q, vol_fn, n_space=n_space, n_time=n_time,
                              fd_bump=h, time_power=2.0)
        conv.append({"n_space": n_space, "n_time": n_time, "fd_bump": h,
                     "delta_grid": res["delta_grid"], "gamma_grid": res["gamma_grid"],
                     "delta_fd": res["delta_fd"], "gamma_fd": res["gamma_fd"]})

    off_corner = []
    for s0 in (90.0, 95.0, 105.0, 110.0):
        res = pde_delta_gamma(s0, params.S0, 1.0, r, q, vol_fn, n_space=1201, n_time=800,
                              fd_bump=0.5, time_power=2.0)
        off_corner.append({"S0": s0, "gamma_grid": res["gamma_grid"], "gamma_fd": res["gamma_fd"],
                           "rel_diff": abs(res["gamma_grid"] - res["gamma_fd"]) / res["gamma_fd"]})

    tmin_sens = []
    for t_min in (1e-3, 1e-4, 1e-5, 1e-6):
        fn = lambda s, t, tm=t_min: dupire_local_vol_from_spot(s, t, params, t_min=tm)
        res = pde_delta_gamma(params.S0, params.S0, 1.0, r, q, fn, n_space=1201, n_time=800,
                              fd_bump=0.5, time_power=2.0)
        tmin_sens.append({"t_min": t_min, "price": res["price"], "gamma_grid": res["gamma_grid"],
                          "gamma_fd": res["gamma_fd"]})

    # Pembanding independen: gamma STRIKE (C_KK) dari solver forward vs rumus SSVI analitik (densitas)
    surf = solve_forward_dupire_surface(params.S0, [1.0], r, q, vol_fn, n_space=2401,
                                        n_time_per_year=1600, time_power=2.0)
    k_grid = surf["K"]
    c_grid = surf["snapshots"][1.0]
    a1, b1, c1, a2, b2, c2 = nonuniform_fd_weights(k_grid)
    ckk = a2 * c_grid[:-2] + b2 * c_grid[1:-1] + c2 * c_grid[2:]
    idx = int(np.argmin(np.abs(k_grid[1:-1] - params.S0)))
    k_node = float(k_grid[1:-1][idx])

    def call_k(k_val: float) -> float:
        y = log_moneyness(k_val, 1.0, params.S0, r, q)
        return float(bs_call_from_total_variance(params.S0, k_val, 1.0, r, q, ssvi_total_variance(y, 1.0, params)))

    hk = 0.5
    ckk_analytic = (-call_k(k_node + 2 * hk) + 16 * call_k(k_node + hk) - 30 * call_k(k_node)
                    + 16 * call_k(k_node - hk) - call_k(k_node - 2 * hk)) / (12 * hk * hk)
    return {
        "convergence_at_S0": conv,
        "off_corner_points": off_corner,
        "t_min_sensitivity": tmin_sens,
        "strike_gamma_forward_node_K": k_node,
        "strike_gamma_forward": float(ckk[idx]),
        "strike_gamma_ssvi_analytic": float(ckk_analytic),
    }
