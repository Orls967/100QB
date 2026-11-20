"""Runner lengkap day003a: menjalankan Bagian 1-4, mencatat waktu, menyimpan results.json & grafik."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from qb100 import PROJECT_ROOT

from .barrier import (
    mc_down_and_out_brownian_bridge,
    reiner_rubinstein_down_call,
    solve_backward_down_and_out_call,
)
from .bs import bs_call_from_total_variance, bs_call_greeks, bs_call_price, implied_vol_call
from .calibration import calibrate_ssvi_lm, generate_noisy_ssvi_surface, run_calibration_mc_study
from .dupire import (
    dupire_local_var_from_prices_fd4,
    dupire_local_var_ssvi,
    dupire_local_vol_from_spot,
    dupire_local_vol_ssvi,
)
from .greeks import gamma_diagnosis_ssvi, pde_delta_gamma, vega_bucket_analysis
from .mc import mc_euler_european_calls, mc_euler_weak_convergence
from .pde_backward import solve_backward_european_call
from .pde_forward import eval_forward_surface_at_strikes, solve_forward_dupire_surface
from .plotting import make_all_figures
from .ssvi import (
    SSVIParams,
    check_butterfly_arbitrage,
    check_calendar_arbitrage,
    check_gatheral_jacquier_sufficient,
    forward_price,
    log_moneyness,
    ssvi_butterfly_g,
    ssvi_d2w_dy2,
    ssvi_dw_dT,
    ssvi_dw_dy,
    ssvi_theta_prime,
    ssvi_total_variance,
)


def _format_17(obj):
    """Konversi ndarray/float ke tipe Python standar; float akan diserialisasi pada presisi penuh 17 digit."""
    if isinstance(obj, np.ndarray):
        return [_format_17(x) for x in obj.tolist()]
    if isinstance(obj, (np.floating, float)):
        return float(f"{float(obj):.17g}")
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, dict):
        return {str(k): _format_17(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [_format_17(v) for v in obj]
    return obj


def run_all_day003a(save_outputs: bool = True) -> dict:
    """Jalankan seluruh pipeline day003a dan kembalikan dictionary hasil lengkap."""
    params = SSVIParams()
    timings: dict[str, float] = {}

    # === BAGIAN 1: Permukaan SSVI ===
    t0 = time.perf_counter()
    y_grid = np.linspace(-1.5, 1.5, 61)
    T_grid = np.linspace(0.05, 2.0, 40)
    Y, TT = np.meshgrid(y_grid, T_grid, indexing="ij")

    hy, hy2, hT = 1e-5, 1.5e-3, 1e-6
    wy_ana = ssvi_dw_dy(Y, TT, params)
    wy_fd = (ssvi_total_variance(Y + hy, TT, params) - ssvi_total_variance(Y - hy, TT, params)) / (2.0 * hy)
    rel_wy = float(np.max(np.abs(wy_ana - wy_fd) / np.maximum(np.abs(wy_fd), 1e-12)))

    wyy_ana = ssvi_d2w_dy2(Y, TT, params)
    wyy_fd = (
        -ssvi_total_variance(Y + 2.0 * hy2, TT, params)
        + 16.0 * ssvi_total_variance(Y + hy2, TT, params)
        - 30.0 * ssvi_total_variance(Y, TT, params)
        + 16.0 * ssvi_total_variance(Y - hy2, TT, params)
        - ssvi_total_variance(Y - 2.0 * hy2, TT, params)
    ) / (12.0 * hy2 * hy2)
    rel_wyy = float(np.max(np.abs(wyy_ana - wyy_fd) / np.maximum(np.abs(wyy_fd), 1e-8)))

    wT_ana = ssvi_dw_dT(Y, TT, params)
    wT_fd = (ssvi_total_variance(Y, TT + hT, params) - ssvi_total_variance(Y, TT - hT, params)) / (2.0 * hT)
    rel_wT = float(np.max(np.abs(wT_ana - wT_fd) / np.maximum(np.abs(wT_fd), 1e-12)))

    bf = check_butterfly_arbitrage(y_grid, T_grid, params)
    cal = check_calendar_arbitrage(y_grid, T_grid, params)
    gj = check_gatheral_jacquier_sufficient(T_grid, params)

    bad_params = SSVIParams(eta=4.0, rho=-0.9)
    bad_bf = check_butterfly_arbitrage(np.linspace(-2.5, 2.5, 101), T_grid, bad_params)
    g_bad_slice = ssvi_butterfly_g(y_grid, 1.0, bad_params)

    # Monotonisitas & konveksitas call
    K_chk = np.linspace(40.0, 180.0, 71)
    max_dC = -np.inf
    min_d2C = np.inf
    for T_c in (0.1, 0.5, 1.0, 2.0):
        yc = log_moneyness(K_chk, T_c, S0=params.S0, r=params.r, q=params.q)
        wc = ssvi_total_variance(yc, T_c, params)
        Cc = bs_call_from_total_variance(params.S0, K_chk, T_c, params.r, params.q, wc)
        max_dC = max(max_dC, float(np.max(np.diff(Cc))))
        min_d2C = min(min_d2C, float(np.min(np.diff(Cc, n=2))))

    timings["part1_ssvi"] = time.perf_counter() - t0

    # === BAGIAN 2: Volatilitas Lokal Dupire ===
    t0 = time.perf_counter()
    params_flat = SSVIParams(eta=0.0)
    loc_var_eta0 = dupire_local_var_ssvi(Y, TT, params_flat)
    theta_p_eta0 = ssvi_theta_prime(TT, params_flat)
    eta0_max_err = float(np.max(np.abs(loc_var_eta0 - theta_p_eta0)))

    fd4_checks = []
    max_fd4_rel_err = 0.0
    for T_c in (0.25, 0.5, 1.0, 1.5):
        Fc = forward_price(T_c, S0=params.S0, r=params.r, q=params.q)
        for mc_val in (0.8, 0.9, 1.0, 1.1, 1.2):
            Kc = mc_val * Fc
            yc = log_moneyness(Kc, T_c, S0=params.S0, r=params.r, q=params.q)
            va = float(dupire_local_var_ssvi(yc, T_c, params))
            vf = dupire_local_var_from_prices_fd4(Kc, T_c, params)
            re = abs(vf - va) / va
            max_fd4_rel_err = max(max_fd4_rel_err, re)
            fd4_checks.append({"T": T_c, "m": mc_val, "var_ana": va, "var_fd4": vf, "rel_err": re})

    loc_var_grid = dupire_local_var_ssvi(Y, TT, params)
    sigma_loc_grid = dupire_local_vol_ssvi(Y, TT, params)
    w_grid = ssvi_total_variance(Y, TT, params)
    iv_grid = np.sqrt(w_grid / TT)
    g_grid = ssvi_butterfly_g(Y, TT, params)
    timings["part2_dupire"] = time.perf_counter() - t0

    # === BAGIAN 3: Dua Solver PDE & Monte Carlo Euler ===
    t0 = time.perf_counter()
    # 3.3 Uji kasus datar
    flat_sig = 0.2
    flat_vol_fn = lambda s, t: np.full_like(s, flat_sig, dtype=np.float64)
    exact_flat = float(bs_call_price(params.S0, 105.0, 1.0, params.r, params.q, flat_sig))
    resolutions = [50, 100, 200, 400]
    err_bwd_flat, err_fwd_flat = [], []
    for n in resolutions:
        pb = float(solve_backward_european_call(params.S0, 105.0, 1.0, params.r, params.q, flat_vol_fn, n_space=n + 1, n_time=n))
        err_bwd_flat.append(abs(pb - exact_flat))
        sf = solve_forward_dupire_surface(params.S0, [1.0], params.r, params.q, flat_vol_fn, n_space=n + 1, n_time_per_year=n)
        pf = float(eval_forward_surface_at_strikes(sf, 1.0, [105.0])[0])
        err_fwd_flat.append(abs(pf - exact_flat))
    slope_bwd = float(np.polyfit(np.log(resolutions), np.log(err_bwd_flat), 1)[0])
    slope_fwd = float(np.polyfit(np.log(resolutions), np.log(err_fwd_flat), 1)[0])

    # 3.4 Uji permukaan SSVI pada resolusi halus
    vol_fn_ssvi = lambda s, t: dupire_local_vol_from_spot(s, t, params)
    T_vals = [0.25, 0.5, 1.0, 2.0]
    m_vals = [0.8, 0.9, 1.0, 1.1, 1.2]

    surf_ssvi = solve_forward_dupire_surface(
        params.S0, T_vals, params.r, params.q, vol_fn_ssvi,
        n_space=5201, n_time_per_year=3600, time_power=2.0,
    )
    iv_err_bwd_bp = np.zeros((len(m_vals), len(T_vals)), dtype=np.float64)
    iv_err_fwd_bp = np.zeros((len(m_vals), len(T_vals)), dtype=np.float64)
    price_diff_mat = np.zeros((len(m_vals), len(T_vals)), dtype=np.float64)
    pde_prices_T1 = []

    for j, T_c in enumerate(T_vals):
        Fc = forward_price(T_c, S0=params.S0, r=params.r, q=params.q)
        strikes = [m * Fc for m in m_vals]
        p_fwds = eval_forward_surface_at_strikes(surf_ssvi, T_c, strikes)
        for i, (m_c, Kc, pf) in enumerate(zip(m_vals, strikes, p_fwds)):
            yc = log_moneyness(Kc, T_c, S0=params.S0, r=params.r, q=params.q)
            iv_true = float(np.sqrt(ssvi_total_variance(yc, T_c, params) / T_c))
            pb = float(
                solve_backward_european_call(
                    params.S0, Kc, T_c, params.r, params.q, vol_fn_ssvi,
                    n_space=5201, n_time=max(int(3600 * T_c), 1200), time_power=2.0,
                )
            )
            iv_b = implied_vol_call(pb, params.S0, Kc, T_c, params.r, params.q)
            iv_f = implied_vol_call(float(pf), params.S0, Kc, T_c, params.r, params.q)
            iv_err_bwd_bp[i, j] = abs(iv_b - iv_true) * 1e4
            iv_err_fwd_bp[i, j] = abs(iv_f - iv_true) * 1e4
            price_diff_mat[i, j] = abs(pb - float(pf))
            if abs(T_c - 1.0) < 1e-12:
                pde_prices_T1.append(pb)

    timings["part3_pde"] = time.perf_counter() - t0

    # 3.5 Monte Carlo Euler
    t0 = time.perf_counter()
    F1 = forward_price(1.0, S0=params.S0, r=params.r, q=params.q)
    strikes_T1 = [m * F1 for m in m_vals]
    mc_res = mc_euler_european_calls(
        params.S0, strikes_T1, 1.0, params.r, params.q, vol_fn_ssvi,
        n_paths=400_000, steps_per_year=200, seed=20261121,
    )
    z_scores_mc = np.abs(mc_res["prices"] - np.asarray(pde_prices_T1)) / mc_res["se"]

    K_weak = float(strikes_T1[4])  # K/F = 1.2 agar bias Euler terlihat jelas vs noise
    ref_weak = float(pde_prices_T1[4])
    weak_res = mc_euler_weak_convergence(
        params.S0, K_weak, 1.0, params.r, params.q, vol_fn_ssvi,
        ref_price=ref_weak, steps_list=(25, 50, 100, 200), n_paths=400_000, seed=20261121,
    )
    timings["part3_mc"] = time.perf_counter() - t0

    # Diagnosis (A1): bias Euler vs jumlah langkah, CRN, kelima strike T=1 (bukan penyetelan)
    t0 = time.perf_counter()
    bias_steps = (50, 100, 200, 400)
    bias_rows = []
    for i_k, k_val in enumerate(strikes_T1):
        wb = mc_euler_weak_convergence(
            params.S0, float(k_val), 1.0, params.r, params.q, vol_fn_ssvi,
            ref_price=float(pde_prices_T1[i_k]), steps_list=bias_steps, n_paths=400_000, seed=20261121,
        )
        bias_rows.append({
            "m": m_vals[i_k], "K": float(k_val), "pde_price": float(pde_prices_T1[i_k]),
            "steps": list(bias_steps), "mc_prices": wb["prices"].tolist(), "mc_se": wb["se"].tolist(),
            "signed_bias": wb["signed_errors"].tolist(),
            "z_signed": (wb["signed_errors"] / wb["se"]).tolist(),
        })
    timings["part3_mc_bias_diag"] = time.perf_counter() - t0

    # === BAGIAN 4: Barrier, Greeks, Kalibrasi ===
    t0 = time.perf_counter()
    H_bar, K_bar, T_bar = 85.0, 100.0, 1.0
    rr_flat = reiner_rubinstein_down_call(params.S0, K_bar, H_bar, T_bar, params.r, params.q, flat_sig)
    rr_tiny = reiner_rubinstein_down_call(params.S0, K_bar, 1e-6, T_bar, params.r, params.q, flat_sig)
    vanilla_flat = float(bs_call_price(params.S0, K_bar, T_bar, params.r, params.q, flat_sig))
    pde_bar_flat = solve_backward_down_and_out_call(
        params.S0, K_bar, H_bar, T_bar, params.r, params.q, flat_vol_fn, n_space=1601, n_time=800
    )

    pde_bar_ssvi = solve_backward_down_and_out_call(
        params.S0, K_bar, H_bar, T_bar, params.r, params.q, vol_fn_ssvi,
        n_space=801, n_time=600, time_power=2.0,
    )
    mc_bar_ssvi = mc_down_and_out_brownian_bridge(
        params.S0, K_bar, H_bar, T_bar, params.r, params.q, vol_fn_ssvi,
        n_paths=400_000, steps_per_year=200, seed=20261121,
    )
    bar_z = abs(mc_bar_ssvi["price"] - pde_bar_ssvi) / mc_bar_ssvi["se"]
    timings["part4_barrier"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    dg_flat = pde_delta_gamma(params.S0, 100.0, 1.0, params.r, params.q, flat_vol_fn, n_space=801, n_time=400, fd_bump=0.05)
    bs_g_flat = bs_call_greeks(params.S0, 100.0, 1.0, params.r, params.q, flat_sig)
    dg_ssvi = pde_delta_gamma(params.S0, 100.0, 1.0, params.r, params.q, vol_fn_ssvi, n_space=601, n_time=400, fd_bump=0.05, time_power=2.0)
    vb_res = vega_bucket_analysis(params=params, T=1.0, n_space=401, n_time=240)
    timings["part4_greeks"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    gamma_diag = gamma_diagnosis_ssvi(params)
    timings["part4_gamma_diag"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    noisy_mkt = generate_noisy_ssvi_surface(params=params, noise_std=5e-4, seed=20261121)
    fit_single = calibrate_ssvi_lm(noisy_mkt["iv_noisy"], noisy_mkt["y_vals"], noisy_mkt["T_grid"])
    mc_calib = run_calibration_mc_study(params=params, n_rep=50, noise_std=5e-4, base_seed=20261121)
    timings["part4_calibration"] = time.perf_counter() - t0

    timings["total"] = sum(timings.values())

    results = {
        "part1": {
            "derivatives_rel_err": {"dw_dy": rel_wy, "d2w_dy2": rel_wyy, "dw_dT": rel_wT},
            "butterfly": bf,
            "calendar": cal,
            "gatheral_jacquier_sufficient": gj,
            "counterfactual_bad": bad_bf,
            "call_monotonicity_max_dC": max_dC,
            "call_convexity_min_d2C": min_d2C,
        },
        "part2": {
            "eta_zero_max_abs_err": eta0_max_err,
            "dupire_fd4_max_rel_err": max_fd4_rel_err,
            "dupire_fd4_checks": fd4_checks,
            "min_sigma_loc": float(np.min(sigma_loc_grid)),
            "max_sigma_loc": float(np.max(sigma_loc_grid)),
            "min_sigma_loc_sq": float(np.min(loc_var_grid)),
        },
        "part3": {
            "flat_convergence": {
                "resolutions": resolutions,
                "backward_errors": err_bwd_flat,
                "forward_errors": err_fwd_flat,
                "backward_slope": slope_bwd,
                "forward_slope": slope_fwd,
            },
            "ssvi_pde": {
                "m_vals": m_vals,
                "T_vals": T_vals,
                "iv_err_backward_bp": iv_err_bwd_bp.tolist(),
                "iv_err_forward_bp": iv_err_fwd_bp.tolist(),
                "max_iv_err_backward_bp": float(np.max(iv_err_bwd_bp)),
                "max_iv_err_forward_bp": float(np.max(iv_err_fwd_bp)),
                "price_diff_matrix": price_diff_mat.tolist(),
                "max_price_diff": float(np.max(price_diff_mat)),
                "pde_prices_T1": pde_prices_T1,
            },
            "mc_euler": {
                "strikes_T1": strikes_T1,
                "pde_prices": pde_prices_T1,
                "mc_prices": mc_res["prices"].tolist(),
                "mc_se": mc_res["se"].tolist(),
                "z_scores": z_scores_mc.tolist(),
                "max_z_score": float(np.max(z_scores_mc)),
            },
            "mc_euler_bias_vs_steps": bias_rows,
            "mc_weak_order": {
                "K_test": K_weak,
                "ref_price": ref_weak,
                "steps": weak_res["steps"],
                "prices": weak_res["prices"].tolist(),
                "se": weak_res["se"].tolist(),
                "abs_errors": weak_res["abs_errors"].tolist(),
                "slope": weak_res["slope"],
                "target_slope": weak_res["target_slope"],
                "tolerance": weak_res["tolerance"],
                "within_tolerance": weak_res["within_tolerance"],
            },
        },
        "part4": {
            "barrier": {
                "flat_rr_down_out": rr_flat["down_out"],
                "flat_rr_down_in": rr_flat["down_in"],
                "flat_rr_parity_err": abs((rr_flat["down_out"] + rr_flat["down_in"]) - vanilla_flat),
                "flat_rr_limit_H0_err": abs(rr_tiny["down_out"] - vanilla_flat),
                "flat_pde_down_out": pde_bar_flat,
                "flat_pde_vs_rr_err": abs(pde_bar_flat - rr_flat["down_out"]),
                "ssvi_pde_down_out": pde_bar_ssvi,
                "ssvi_mc_bridge_price": mc_bar_ssvi["price"],
                "ssvi_mc_bridge_se": mc_bar_ssvi["se"],
                "ssvi_z_score": bar_z,
            },
            "greeks": {
                "flat_pde": dg_flat,
                "flat_bs_exact": bs_g_flat,
                "flat_delta_grid_err": abs(dg_flat["delta_grid"] - bs_g_flat["delta"]),
                "flat_gamma_grid_err": abs(dg_flat["gamma_grid"] - bs_g_flat["gamma"]),
                "flat_delta_fd_err": abs(dg_flat["delta_fd"] - bs_g_flat["delta"]),
                "flat_gamma_fd_err": abs(dg_flat["gamma_fd"] - bs_g_flat["gamma"]),
                "ssvi_pde": dg_ssvi,
                "ssvi_gamma_diagnosis": gamma_diag,
            },
            "vega_bucket": vb_res,
            "calibration": {
                "true_params": mc_calib["true_params"].tolist(),
                "single_fit_params": fit_single["beta"].tolist(),
                "single_fit_rms_vol": fit_single["rms_vol_error"],
                "single_fit_iterations": fit_single["iterations"],
                "correlation_matrix": fit_single["corr"].tolist(),
                "mc_r50_mean": mc_calib["mean"].tolist(),
                "mc_r50_sd": mc_calib["sd"].tolist(),
                "mc_r50_bias": mc_calib["bias"].tolist(),
                "mc_r50_bias_over_sd": mc_calib["bias_over_sd"].tolist(),
                "mc_r50_mean_rms": mc_calib["mean_rms"],
            },
        },
        "timings_seconds": timings,
        "_plot_data": {
            "surfaces": {
                "y_grid": y_grid,
                "T_grid": T_grid,
                "iv_grid": iv_grid,
                "sigma_loc_grid": sigma_loc_grid,
                "g_grid": g_grid,
                "g_bad_slice": g_bad_slice,
            },
            "calib_estimates": mc_calib["estimates"],
        },
    }

    if save_outputs:
        out_dir = PROJECT_ROOT / "outputs" / "day003a"
        out_dir.mkdir(parents=True, exist_ok=True)
        clean_json = _format_17(results)
        (out_dir / "results.json").write_text(
            json.dumps(clean_json, indent=2), encoding="utf-8"
        )
        fig_paths = make_all_figures(results)
        results["figure_paths"] = [str(p.relative_to(PROJECT_ROOT)) for p in fig_paths]

    return results
