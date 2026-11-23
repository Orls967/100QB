use std::path::Path;
use std::time::Instant;

use day005d_merton_jump_diffusion::bs::implied_vol_call;
use day005d_merton_jump_diffusion::calib::{
    build_smile_table, calibrate_merton_nelder_mead, run_multistart_calibration,
    run_noisy_calibration_study, NoisyCalibConfig,
};
use day005d_merton_jump_diffusion::jsonout::{obj, write_json_file, JsonVal};
use day005d_merton_jump_diffusion::merton::{
    merton_call_series, short_maturity_first_order_put_iv, simulate_merton_mc, MertonParams,
};
use day005d_merton_jump_diffusion::parallel::default_workers;
use day005d_merton_jump_diffusion::pide::{compare_pide_vs_exact, pide_grid_convergence};
use day005d_merton_jump_diffusion::svg::{write_line_chart, Series};
use day005d_merton_jump_diffusion::tail_risk::{compare_otm_puts, simulate_short_put_tail_risk};

const MAIN_SEED: u64 = 20261123;

fn criterion(status: &str, measured: JsonVal, threshold: &str, note: &str) -> JsonVal {
    obj([
        ("status", status.into()),
        ("measured", measured),
        ("threshold", threshold.into()),
        ("note", note.into()),
    ])
}

fn main() {
    let t_total = Instant::now();
    let n_workers = default_workers();
    let p = MertonParams::default();
    println!("=== Day 005d: Jump-Diffusion Merton (workers = {n_workers}) ===");

    // -------------------------------------------------------------------------
    // Bagian 1: Seri eksak Merton & verifikasi MC (R = 400_000)
    // -------------------------------------------------------------------------
    let t1 = Instant::now();
    let ks_mc = [80.0, 90.0, 100.0, 110.0, 120.0];
    let mc_res = simulate_merton_mc(&p, &ks_mc, 1.0, 400_000, MAIN_SEED, n_workers);
    let dt_p1 = t1.elapsed().as_secs_f64();

    let max_mc_call_z = mc_res
        .rows
        .iter()
        .map(|r| r.call_z.abs())
        .fold(0.0_f64, f64::max);
    let max_mc_put_z = mc_res
        .rows
        .iter()
        .map(|r| r.put_z.abs())
        .fold(0.0_f64, f64::max);
    println!(
        "[Bagian 1] selesai dalam {dt_p1:.2}s | martingale_z = {:.2}, max |call_z| = {max_mc_call_z:.2}, max |put_z| = {max_mc_put_z:.2}",
        mc_res.martingale_z
    );

    // -------------------------------------------------------------------------
    // Bagian 2: PIDE Crank-Nicolson vs Seri Eksak + Konvergensi Grid
    // -------------------------------------------------------------------------
    let t2 = Instant::now();
    let pide_rows = compare_pide_vs_exact(&p, &ks_mc, 1.0, 401, 200);
    let conv_rows = pide_grid_convergence(&p, &ks_mc, 1.0, &[(101, 50), (201, 100), (401, 200)]);
    let dt_p2 = t2.elapsed().as_secs_f64();

    let max_pide_abs_err = pide_rows.iter().map(|r| r.abs_err).fold(0.0_f64, f64::max);
    let max_pide_rel_err = pide_rows.iter().map(|r| r.rel_err).fold(0.0_f64, f64::max);
    let conv_ratio = conv_rows[0].max_abs_err / conv_rows[2].max_abs_err.max(1e-15);
    println!(
        "[Bagian 2] selesai dalam {dt_p2:.2}s | max PIDE abs_err = {max_pide_abs_err:.5} (rel = {max_pide_rel_err:.5}), conv 101->401 ratio = {conv_ratio:.2}x"
    );

    // -------------------------------------------------------------------------
    // Bagian 3: Implied Volatility Smile & Kalibrasi Nelder-Mead
    // -------------------------------------------------------------------------
    let t3 = Instant::now();
    let ks_smile = [70.0, 80.0, 90.0, 95.0, 100.0, 105.0, 110.0, 120.0, 130.0];
    let ts_smile = [0.25, 1.0, 2.0];
    let smile_rows = build_smile_table(&p, &ts_smile, &ks_smile);
    let p_init = MertonParams {
        sigma: 0.20,
        lambda: 0.45,
        mu_j: -0.04,
        delta_j: 0.18,
        ..p
    };
    let calib = calibrate_merton_nelder_mead(
        &p,
        &p_init,
        &[0.25, 1.0],
        &[80.0, 90.0, 100.0, 110.0, 120.0],
        800,
    );
    let ks_calib = [80.0, 90.0, 100.0, 110.0, 120.0];
    let noisy_1mat = run_noisy_calibration_study(
        &p,
        &p_init,
        &[0.25],
        &ks_calib,
        NoisyCalibConfig {
            iv_noise_std: 0.0010,
            r_reps: 200,
            main_seed: MAIN_SEED,
            seed_stream: 310,
            n_workers,
        },
    );
    let noisy_2mat = run_noisy_calibration_study(
        &p,
        &p_init,
        &[0.25, 1.0],
        &ks_calib,
        NoisyCalibConfig {
            iv_noise_std: 0.0010,
            r_reps: 200,
            main_seed: MAIN_SEED,
            seed_stream: 320,
            n_workers,
        },
    );
    let multistart = run_multistart_calibration(&p, &[0.25, 1.0], &ks_calib, 24, MAIN_SEED);

    let short_skew_rows: Vec<(f64, f64, f64, f64)> = [0.01, 0.02, 0.05]
        .iter()
        .map(|&t_s| {
            let px = merton_call_series(&p, 80.0, t_s);
            let iv_ex = implied_vol_call(px, p.s0, 80.0, p.r, p.q, t_s);
            let iv_ap = short_maturity_first_order_put_iv(&p, 80.0, t_s);
            let rel_err = (iv_ap - iv_ex).abs() / iv_ex;
            (t_s, iv_ex, iv_ap, rel_err)
        })
        .collect();
    let max_short_skew_rel_err = short_skew_rows.iter().map(|r| r.3).fold(0.0_f64, f64::max);
    let dt_p3 = t3.elapsed().as_secs_f64();

    let iv_short_left = smile_rows
        .iter()
        .find(|r| (r.t - 0.25).abs() < 1e-9 && (r.k - 80.0).abs() < 1e-9)
        .map(|r| r.implied_vol)
        .unwrap_or(0.0);
    let iv_short_atm = smile_rows
        .iter()
        .find(|r| (r.t - 0.25).abs() < 1e-9 && (r.k - 100.0).abs() < 1e-9)
        .map(|r| r.implied_vol)
        .unwrap_or(0.0);
    println!(
        "[Bagian 3] selesai dalam {dt_p3:.2}s | IV(T=0.25, K=80)={:.2}% > ATM={:.2}%, calib price_rmse={:.2e}, iv_rmse={:.3} bps",
        iv_short_left * 100.0,
        iv_short_atm * 100.0,
        calib.price_rmse,
        calib.iv_rmse_bps
    );

    // -------------------------------------------------------------------------
    // Bagian 4: Risiko Ekor (OTM Put & 99% VaR/ES Short Put 1 Bulan)
    // -------------------------------------------------------------------------
    let t4 = Instant::now();
    let otm_puts = compare_otm_puts(&p, &[70.0, 80.0, 90.0], 0.25);
    let tail_risk =
        simulate_short_put_tail_risk(&p, 90.0, 1.0 / 12.0, 200_000, MAIN_SEED, n_workers);
    let dt_p4 = t4.elapsed().as_secs_f64();
    println!(
        "[Bagian 4] selesai dalam {dt_p4:.2}s | K=70 put ratio = {:.2}x BS, 1m K=90 ES99 Merton = {:.3} vs BS-ATM = {:.3} vs BS-VarMatched = {:.3}",
        otm_puts[0].ratio_merton_to_bs,
        tail_risk.es99_merton,
        tail_risk.es99_bs_atm,
        tail_risk.es99_bs_var_matched
    );

    // Worker invariance check (1 vs n_workers)
    let inv1 = simulate_merton_mc(&p, &[100.0], 1.0, 2048, MAIN_SEED, 1);
    let inv_n = simulate_merton_mc(&p, &[100.0], 1.0, 2048, MAIN_SEED, n_workers);
    let worker_bit_diff =
        (inv1.rows[0].call_mc.to_bits() != inv_n.rows[0].call_mc.to_bits()) as usize;

    // -------------------------------------------------------------------------
    // SVG Figures
    // -------------------------------------------------------------------------
    let fig_dir = Path::new("outputs/figures");
    // 1. PIDE vs Exact Call Prices
    let pide_k: Vec<f64> = pide_rows.iter().map(|r| r.k).collect();
    let pide_ex: Vec<f64> = pide_rows.iter().map(|r| r.exact_call).collect();
    let pide_cn: Vec<f64> = pide_rows.iter().map(|r| r.pide_call).collect();
    write_line_chart(
        &fig_dir.join("day005d_pide_vs_exact.svg"),
        "Harga Call Merton (T=1): Seri Eksak vs PIDE Crank-Nicolson",
        "Strike K",
        "Harga Call",
        &[
            Series {
                label: "Seri Eksak",
                x: &pide_k,
                y: &pide_ex,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "PIDE CN (401x200)",
                x: &pide_k,
                y: &pide_cn,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .expect("gagal menulis SVG pide_vs_exact");

    // 2. PIDE Grid Convergence
    let gx: Vec<f64> = conv_rows.iter().map(|r| (r.n_x as f64).log10()).collect();
    let gy: Vec<f64> = conv_rows.iter().map(|r| r.max_abs_err.log10()).collect();
    write_line_chart(
        &fig_dir.join("day005d_pide_convergence.svg"),
        "Konvergensi Grid PIDE Crank-Nicolson (log10 Max Abs Error vs log10 N_x)",
        "log10(N_x)",
        "log10(Max Abs Error)",
        &[Series {
            label: "Max Abs Error",
            x: &gx,
            y: &gy,
            color: "#dc2626",
            dashed: false,
        }],
    )
    .expect("gagal menulis SVG pide_convergence");

    // 3. Implied Volatility Smile across maturities
    let smile_slice = |t_val: f64| -> Vec<f64> {
        smile_rows
            .iter()
            .filter(|r| (r.t - t_val).abs() < 1e-9)
            .map(|r| r.implied_vol * 100.0)
            .collect()
    };
    let iv_025 = smile_slice(0.25);
    let iv_100 = smile_slice(1.0);
    let iv_200 = smile_slice(2.0);
    write_line_chart(
        &fig_dir.join("day005d_iv_smile.svg"),
        "Implied Volatility Smile Merton (T = 0.25, 1.0, 2.0 Tahun)",
        "Strike K",
        "Implied Volatility (%)",
        &[
            Series {
                label: "T = 0.25y",
                x: &ks_smile,
                y: &iv_025,
                color: "#dc2626",
                dashed: false,
            },
            Series {
                label: "T = 1.00y",
                x: &ks_smile,
                y: &iv_100,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "T = 2.00y",
                x: &ks_smile,
                y: &iv_200,
                color: "#059669",
                dashed: true,
            },
        ],
    )
    .expect("gagal menulis SVG iv_smile");

    // 4. Deep OTM Put Merton vs BS ATM IV
    let ok_x: Vec<f64> = otm_puts.iter().map(|r| r.k).collect();
    let ok_m: Vec<f64> = otm_puts.iter().map(|r| r.merton_put).collect();
    let ok_b: Vec<f64> = otm_puts.iter().map(|r| r.bs_atm_put).collect();
    write_line_chart(
        &fig_dir.join("day005d_otm_put_tail.svg"),
        "Harga Put OTM (T = 0.25y): Merton vs Black-Scholes ATM IV",
        "Strike K",
        "Harga Put",
        &[
            Series {
                label: "Merton Jump-Diffusion",
                x: &ok_x,
                y: &ok_m,
                color: "#dc2626",
                dashed: false,
            },
            Series {
                label: "Black-Scholes (ATM IV)",
                x: &ok_x,
                y: &ok_b,
                color: "#2563eb",
                dashed: true,
            },
        ],
    )
    .expect("gagal menulis SVG otm_put_tail");

    let total_elapsed = t_total.elapsed().as_secs_f64();
    println!("Total waktu eksekusi: {total_elapsed:.2}s");

    // -------------------------------------------------------------------------
    // Criteria evaluation (a-priori thresholds)
    // -------------------------------------------------------------------------
    let c1_ok = mc_res.martingale_z.abs() <= 3.0 && max_mc_call_z <= 3.0 && max_mc_put_z <= 3.0;
    let c2_ok = max_pide_abs_err < 0.01 && conv_ratio > 8.0;
    let c3_ok = iv_short_left > iv_short_atm + 0.02 && max_short_skew_rel_err < 0.05;
    let c4_ok = calib.price_rmse < 1e-2
        && calib.iv_rmse_bps < 5.0
        && noisy_2mat.sigma_bias.abs() < 0.02
        && noisy_2mat.lambda_bias.abs() < 0.15
        && noisy_2mat.lambda_sd < noisy_1mat.lambda_sd
        && multistart.global_fraction >= 0.50;
    let c5_ok = otm_puts[0].ratio_merton_to_bs > 2.0
        && tail_risk.es99_merton > tail_risk.es99_bs_var_matched;
    let c6_ok = worker_bit_diff == 0;
    let c7_ok = total_elapsed < 240.0;

    let criteria = obj([
        (
            "c1_merton_series_vs_mc_3se",
            criterion(
                if c1_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("martingale_z", mc_res.martingale_z.into()),
                    ("max_call_z", max_mc_call_z.into()),
                    ("max_put_z", max_mc_put_z.into()),
                ]),
                "martingale_z, max_call_z, max_put_z <= 3.0 (R=400k)",
                "Verifikasi martingal dan kecocokan seri eksak Merton vs MC",
            ),
        ),
        (
            "c2_pide_cn_matches_exact_and_converges",
            criterion(
                if c2_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("max_pide_abs_err", max_pide_abs_err.into()),
                    ("max_pide_rel_err", max_pide_rel_err.into()),
                    ("conv_ratio_101_to_401", conv_ratio.into()),
                ]),
                "max_pide_abs_err < 0.01 dan conv_ratio_101_to_401 > 8.0",
                "Akurasi PIDE Crank-Nicolson + Rannacher dan konvergensi orde ~2",
            ),
        ),
        (
            "c3_short_maturity_jump_skew",
            criterion(
                if c3_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("iv_t025_k80", iv_short_left.into()),
                    ("iv_t025_k100", iv_short_atm.into()),
                    ("max_rel_err_t001_005", max_short_skew_rel_err.into()),
                ]),
                "IV(T=0.25, K=80) > IV(T=0.25, K=100) + 0.02 dan max_rel_err_t001_005 < 0.05",
                "Skew jangka pendek sesuai aproksimasi lompatan tunggal urutan pertama O(T)",
            ),
        ),
        (
            "c4_nelder_mead_calibration_recovery",
            criterion(
                if c4_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("price_rmse", calib.price_rmse.into()),
                    ("iv_rmse_bps", calib.iv_rmse_bps.into()),
                    ("iters", calib.iters.into()),
                    ("noisy_2mat_sigma_bias", noisy_2mat.sigma_bias.into()),
                    ("noisy_2mat_lambda_bias", noisy_2mat.lambda_bias.into()),
                    ("noisy_1mat_lambda_sd", noisy_1mat.lambda_sd.into()),
                    ("noisy_2mat_lambda_sd", noisy_2mat.lambda_sd.into()),
                    ("multistart_global_fraction", multistart.global_fraction.into()),
                ]),
                "price_rmse < 1e-2, iv_rmse_bps < 5.0, |sigma_bias| < 0.02, |lambda_bias| < 0.15, 2mat_lambda_sd < 1mat_lambda_sd, multistart_global_fraction >= 0.50",
                "Kalibrasi Nelder-Mead pada permukaan bebas-derau, kuotasi IV berderau R=200 (1-mat vs 2-mat), dan multi-start 24 titik",
            ),
        ),
        (
            "c5_tail_risk_otm_put_and_es99",
            criterion(
                if c5_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    (
                        "k70_put_ratio_merton_vs_bs",
                        otm_puts[0].ratio_merton_to_bs.into(),
                    ),
                    ("es99_merton", tail_risk.es99_merton.into()),
                    ("es99_bs_atm", tail_risk.es99_bs_atm.into()),
                    ("es99_bs_var_matched", tail_risk.es99_bs_var_matched.into()),
                ]),
                "ratio(K=70) > 2.0 dan ES99_merton > ES99_bs_var_matched",
                "Merton tetap memiliki ES99 lebih tinggi bahkan dibanding Black-Scholes variansi total sama",
            ),
        ),
        (
            "c6_worker_invariance_1_vs_n",
            criterion(
                if c6_ok { "tercapai" } else { "belum tercapai" },
                obj([("bit_diffs", worker_bit_diff.into())]),
                "bit_diffs == 0",
                "Reduksi blok tetap menjamin hasil bit-identik 1 vs N worker",
            ),
        ),
        (
            "c7_runtime_under_240s",
            criterion(
                if c7_ok { "tercapai" } else { "belum tercapai" },
                obj([("total_seconds", total_elapsed.into())]),
                "< 240.0s",
                "Anggaran waktu eksekusi runner penuh",
            ),
        ),
    ]);

    let root_json = obj([
        (
            "metadata",
            obj([
                ("day", 5_usize.into()),
                ("session", "d".into()),
                (
                    "topic_source",
                    "backup_roadmap_075_merton_jump_diffusion".into(),
                ),
                ("language", "rust".into()),
                ("main_seed", (MAIN_SEED as usize).into()),
                ("workers_used", n_workers.into()),
                ("run_completed_without_errors", true.into()),
                (
                    "elapsed_seconds",
                    obj([
                        ("part1_series_mc", dt_p1.into()),
                        ("part2_pide_cn", dt_p2.into()),
                        ("part3_smile_calib", dt_p3.into()),
                        ("part4_tail_risk", dt_p4.into()),
                        ("total", total_elapsed.into()),
                    ]),
                ),
            ]),
        ),
        ("criteria", criteria),
        (
            "part1_series_vs_mc",
            obj([
                ("martingale_mean", mc_res.martingale_mean.into()),
                ("martingale_se", mc_res.martingale_se.into()),
                ("martingale_z", mc_res.martingale_z.into()),
                (
                    "rows",
                    JsonVal::Arr(
                        mc_res
                            .rows
                            .iter()
                            .map(|r| {
                                obj([
                                    ("k", r.k.into()),
                                    ("call_exact", r.call_exact.into()),
                                    ("call_mc", r.call_mc.into()),
                                    ("call_se", r.call_se.into()),
                                    ("call_z", r.call_z.into()),
                                    ("put_exact", r.put_exact.into()),
                                    ("put_mc", r.put_mc.into()),
                                    ("put_se", r.put_se.into()),
                                    ("put_z", r.put_z.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "part2_pide_cn",
            obj([
                (
                    "rows",
                    JsonVal::Arr(
                        pide_rows
                            .iter()
                            .map(|r| {
                                obj([
                                    ("k", r.k.into()),
                                    ("exact_call", r.exact_call.into()),
                                    ("pide_call", r.pide_call.into()),
                                    ("abs_err", r.abs_err.into()),
                                    ("rel_err", r.rel_err.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "convergence",
                    JsonVal::Arr(
                        conv_rows
                            .iter()
                            .map(|r| {
                                obj([
                                    ("n_x", r.n_x.into()),
                                    ("n_t", r.n_t.into()),
                                    ("max_abs_err", r.max_abs_err.into()),
                                    ("atm_abs_err", r.atm_abs_err.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "part3_smile_and_calib",
            obj([
                (
                    "smile",
                    JsonVal::Arr(
                        smile_rows
                            .iter()
                            .map(|r| {
                                obj([
                                    ("t", r.t.into()),
                                    ("k", r.k.into()),
                                    ("call_price", r.call_price.into()),
                                    ("implied_vol", r.implied_vol.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "calibration",
                    obj([
                        ("fitted_sigma", calib.fitted_params.sigma.into()),
                        ("fitted_lambda", calib.fitted_params.lambda.into()),
                        ("fitted_mu_j", calib.fitted_params.mu_j.into()),
                        ("fitted_delta_j", calib.fitted_params.delta_j.into()),
                        ("price_rmse", calib.price_rmse.into()),
                        ("iv_rmse_bps", calib.iv_rmse_bps.into()),
                        ("iters", calib.iters.into()),
                    ]),
                ),
                (
                    "short_maturity_first_order_skew",
                    JsonVal::Arr(
                        short_skew_rows
                            .iter()
                            .map(|r| {
                                obj([
                                    ("t", r.0.into()),
                                    ("iv_exact_k80", r.1.into()),
                                    ("iv_first_order_k80", r.2.into()),
                                    ("rel_err", r.3.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "noisy_calib_1mat_t025",
                    obj([
                        ("r_reps", noisy_1mat.r_reps.into()),
                        ("sigma_bias", noisy_1mat.sigma_bias.into()),
                        ("sigma_sd", noisy_1mat.sigma_sd.into()),
                        ("lambda_bias", noisy_1mat.lambda_bias.into()),
                        ("lambda_sd", noisy_1mat.lambda_sd.into()),
                        ("mu_j_bias", noisy_1mat.mu_j_bias.into()),
                        ("mu_j_sd", noisy_1mat.mu_j_sd.into()),
                        ("delta_j_bias", noisy_1mat.delta_j_bias.into()),
                        ("delta_j_sd", noisy_1mat.delta_j_sd.into()),
                        ("corr_lambda_delta_j", noisy_1mat.corr_lambda_delta_j.into()),
                        ("corr_lambda_sigma", noisy_1mat.corr_lambda_sigma.into()),
                    ]),
                ),
                (
                    "noisy_calib_2mat_t025_100",
                    obj([
                        ("r_reps", noisy_2mat.r_reps.into()),
                        ("sigma_bias", noisy_2mat.sigma_bias.into()),
                        ("sigma_sd", noisy_2mat.sigma_sd.into()),
                        ("lambda_bias", noisy_2mat.lambda_bias.into()),
                        ("lambda_sd", noisy_2mat.lambda_sd.into()),
                        ("mu_j_bias", noisy_2mat.mu_j_bias.into()),
                        ("mu_j_sd", noisy_2mat.mu_j_sd.into()),
                        ("delta_j_bias", noisy_2mat.delta_j_bias.into()),
                        ("delta_j_sd", noisy_2mat.delta_j_sd.into()),
                        ("corr_lambda_delta_j", noisy_2mat.corr_lambda_delta_j.into()),
                        ("corr_lambda_sigma", noisy_2mat.corr_lambda_sigma.into()),
                    ]),
                ),
                (
                    "multistart_calib",
                    obj([
                        ("n_starts", multistart.n_starts.into()),
                        ("n_global_converged", multistart.n_global_converged.into()),
                        ("global_fraction", multistart.global_fraction.into()),
                        ("best_price_rmse", multistart.best_price_rmse.into()),
                        (
                            "local_optima_rmses",
                            JsonVal::Arr(
                                multistart
                                    .local_optima_rmses
                                    .iter()
                                    .map(|&v| v.into())
                                    .collect(),
                            ),
                        ),
                    ]),
                ),
            ]),
        ),
        (
            "part4_tail_risk",
            obj([
                (
                    "otm_puts",
                    JsonVal::Arr(
                        otm_puts
                            .iter()
                            .map(|r| {
                                obj([
                                    ("k", r.k.into()),
                                    ("t", r.t.into()),
                                    ("merton_put", r.merton_put.into()),
                                    ("bs_atm_put", r.bs_atm_put.into()),
                                    ("ratio_merton_to_bs", r.ratio_merton_to_bs.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "short_put_1m_k90",
                    obj([
                        ("bs_atm_vol", tail_risk.bs_atm_vol.into()),
                        ("bs_var_matched_vol", tail_risk.bs_var_matched_vol.into()),
                        ("var99_merton", tail_risk.var99_merton.into()),
                        ("es99_merton", tail_risk.es99_merton.into()),
                        ("var99_bs_atm", tail_risk.var99_bs_atm.into()),
                        ("es99_bs_atm", tail_risk.es99_bs_atm.into()),
                        (
                            "var99_bs_var_matched",
                            tail_risk.var99_bs_var_matched.into(),
                        ),
                        ("es99_bs_var_matched", tail_risk.es99_bs_var_matched.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(Path::new("outputs/day005d/results.json"), &root_json)
        .expect("gagal menulis outputs/day005d/results.json");
}
