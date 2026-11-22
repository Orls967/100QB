//! Unified CLI runner for day004h:
//! Part 1: Volterra covariance anchors, 2n x 2n Cholesky exact reference (R=200k, n=50), hybrid kappa=1 bias vs n.
//! Part 2: Variance & price martingale anchors, variance swap, pure MC vs conditional BS (VRF), Cholesky price check.
//! Part 3: Implied vol smile (13 strikes x 5 maturities, R=100k, n=250), ATM skew ~ T^{H-1/2}, Markov H=0.5 compare.
//! Part 4: Synthetic market (R=200k, n=500), Nelder-Mead calibration (R=20k, n=100 CRN), H-profile, 30-rep, T>=0.5 vs T<=0.1.

use std::path::PathBuf;
use std::time::Instant;

use day004h_rough_bergomi::calib::run_calibration_core;
use day004h_rough_bergomi::jsonout::{write_json_file, JsonVal};
use day004h_rough_bergomi::parallel;
use day004h_rough_bergomi::quad::GaussLegendre;
use day004h_rough_bergomi::rbergomi::{
    compare_hybrid_vs_cholesky_ref, run_call_estimators_comparison,
    run_exact_cholesky_martingale_check, run_martingale_anchors, RBergomiParams,
};
use day004h_rough_bergomi::rng::MAIN_SEED;
use day004h_rough_bergomi::smile::{
    compute_smile, run_atm_skew_scaling_study, run_markov_comparison,
};
use day004h_rough_bergomi::svg::{self, LineSeries, PALETTE};
use day004h_rough_bergomi::volterra::{
    exact_cov_closed_form, exact_cov_quadrature, exact_cross_cov, hybrid_bias_study_with_scheme,
    run_cholesky_reference_check, HybridScheme, HybridWeights,
};

fn find_repo_root() -> PathBuf {
    let mut dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    for _ in 0..8 {
        if dir.join("outputs").is_dir() && dir.join("projects").join("rust").is_dir() {
            return dir;
        }
        if !dir.pop() {
            break;
        }
    }
    PathBuf::from(".")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let quick = args.iter().any(|a| a == "--quick" || a == "-quick");
    let n_workers = parallel::default_workers();
    let root = find_repo_root();
    let out_json = if quick {
        std::env::temp_dir().join("day004h_quick_results.json")
    } else {
        root.join("outputs").join("day004h").join("results.json")
    };
    let fig_dir = if quick {
        std::env::temp_dir().join("day004h_quick_figures")
    } else {
        root.join("outputs").join("figures")
    };

    let t_total = Instant::now();
    let p = RBergomiParams::default_params();
    let r_part1 = if quick { 5_000 } else { 200_000 };
    let r_part2 = if quick { 5_000 } else { 100_000 };
    let r_part3 = if quick { 5_000 } else { 100_000 };

    // ==================== BAGIAN 1: Volterra & Hybrid ====================
    let t1 = Instant::now();
    let gl = GaussLegendre::new(64);
    let pairs = [
        (0.05, 0.05),
        (0.05, 0.5),
        (0.1, 0.9),
        (0.25, 0.27),
        (0.5, 1.0),
        (0.85, 0.9),
    ];
    let mut max_cov_quad_rel_err = 0.0_f64;
    let mut cov_pairs_json = Vec::new();
    for &(s, t) in &pairs {
        let cf = exact_cov_closed_form(s, t, p.h);
        let qd = exact_cov_quadrature(s, t, p.h, &gl);
        let rel = (cf - qd).abs() / cf.abs().max(1e-12);
        max_cov_quad_rel_err = max_cov_quad_rel_err.max(rel);
        cov_pairs_json.push(JsonVal::Obj(vec![
            ("s", s.into()),
            ("t", t.into()),
            ("closed_form", cf.into()),
            ("quadrature", qd.into()),
            ("rel_err", rel.into()),
        ]));
    }
    let cross_cov_sample = exact_cross_cov(0.5, 0.1, 0.12, p.h);

    let chol_ref =
        run_cholesky_reference_check(50, 1.0, p.h, r_part1, n_workers, MAIN_SEED, 0x1100);
    // For H=0.5, W_tilde^{0.5} = sum dB_j is an exact linear combination of dB_j (2n x 2n is rank n),
    // so we verify exact Brownian recovery via Cov(s,t) = min(s,t) and hybrid_bias_h05.
    let mut h05_max_cov_err = 0.0_f64;
    for &s in &[0.1_f64, 0.4, 0.9] {
        for &t in &[0.1_f64, 0.4, 0.9] {
            h05_max_cov_err =
                h05_max_cov_err.max((exact_cov_closed_form(s, t, 0.5) - s.min(t)).abs());
        }
    }

    let ns_bias = [50_usize, 100, 250, 500];
    let hyb_bias_re = hybrid_bias_study_with_scheme(
        0.1,
        1.0,
        &ns_bias,
        r_part1 / 2,
        n_workers,
        MAIN_SEED,
        HybridScheme::RightEndpoint,
    );
    let hyb_bias_h01 = hybrid_bias_study_with_scheme(
        0.1,
        1.0,
        &ns_bias,
        r_part1 / 2,
        n_workers,
        MAIN_SEED,
        HybridScheme::BlpOptimal,
    );
    let hyb_bias_vm = hybrid_bias_study_with_scheme(
        0.1,
        1.0,
        &ns_bias,
        r_part1 / 2,
        n_workers,
        MAIN_SEED,
        HybridScheme::VarianceMatching,
    );
    let hyb_bias_h05 = hybrid_bias_study_with_scheme(
        0.5,
        1.0,
        &ns_bias,
        r_part1 / 4,
        n_workers,
        MAIN_SEED,
        HybridScheme::BlpOptimal,
    );
    let part1_secs = t1.elapsed().as_secs_f64();

    // ==================== BAGIAN 2: Martingales & Pricing ====================
    let t2 = Instant::now();
    let mart_ts = [0.25, 0.5, 1.0];
    let mart_anchors = run_martingale_anchors(p, &mart_ts, 250, r_part2, n_workers, MAIN_SEED);
    let chol_mart = run_exact_cholesky_martingale_check(p, r_part1, n_workers, MAIN_SEED);
    let pricing_ts = [0.1, 0.25, 0.5, 1.0];
    let pricing_ks = [-0.2, -0.1, 0.0, 0.1, 0.2];
    let pricing_cmp = run_call_estimators_comparison(
        p,
        &pricing_ts,
        &pricing_ks,
        250,
        r_part2,
        n_workers,
        MAIN_SEED,
    );
    let ref_vs_hyb = compare_hybrid_vs_cholesky_ref(
        p,
        50,
        1.0,
        &[-0.2, 0.0, 0.2],
        r_part2 / 2,
        n_workers,
        MAIN_SEED,
    );
    let part2_secs = t2.elapsed().as_secs_f64();

    // ==================== BAGIAN 3: Smile & ATM Skew ====================
    let t3 = Instant::now();
    let smile_ts = [0.05, 0.1, 0.25, 0.5, 1.0];
    let smile_ks: Vec<f64> = (-6..=6).map(|i| i as f64 * 0.05).collect();
    let smile = compute_smile(
        p, &smile_ts, &smile_ks, 250, r_part3, n_workers, MAIN_SEED, 0x3000,
    );
    let skew_hs = [0.05, 0.1, 0.25, 0.5];
    let skew_study = run_atm_skew_scaling_study(
        p,
        &skew_hs,
        &smile_ts,
        250,
        r_part3 / 2,
        n_workers,
        MAIN_SEED,
    );
    let smile_t01: Vec<_> = smile
        .iter()
        .filter(|pt| (pt.t_horizon - 0.1).abs() < 1e-12)
        .cloned()
        .collect();
    let markov_cmp = run_markov_comparison(
        p,
        &smile_t01,
        &smile_ks,
        250,
        r_part3 / 2,
        n_workers,
        MAIN_SEED,
    );
    let part3_secs = t3.elapsed().as_secs_f64();

    // ==================== BAGIAN 4: Calibration & Identifiability ====================
    let t4 = Instant::now();
    let calib = run_calibration_core(p, n_workers, MAIN_SEED);
    let part4_secs = t4.elapsed().as_secs_f64();

    // ==================== BAGIAN 5: 7 SVG Figures ====================
    let t5 = Instant::now();
    // 1. Exact vs Hybrid covariance curve for t=0.5 against s in [0.05..1.0]
    let s_grid: Vec<f64> = (1..=25).map(|i| i as f64 * 0.04).collect();
    let cov_ex_curve: Vec<f64> = s_grid
        .iter()
        .map(|&s| exact_cov_closed_form(s, 0.5, p.h))
        .collect();
    let hw250 = HybridWeights::new(250, 1.0, p.h);
    let idx_05 = 124_usize; // (125)*0.004 = 0.5
    let cov_hyb_curve: Vec<f64> = s_grid
        .iter()
        .map(|&s| {
            let idx_s = ((s / hw250.dt).round() as usize).clamp(1, 250) - 1;
            let (i, j) = if idx_s <= idx_05 {
                (idx_s, idx_05)
            } else {
                (idx_05, idx_s)
            };
            if i == j {
                let mut acc = hw250.dt_pow_h * hw250.dt_pow_h;
                for m in 0..i {
                    acc += hw250.conv_weights[i + 1 - m] * hw250.conv_weights[i + 1 - m] * hw250.dt;
                }
                acc
            } else {
                let mut acc =
                    hw250.dt_pow_h * hw250.rho_h * hw250.sqrt_dt * hw250.conv_weights[j + 1 - i];
                for m in 0..i {
                    acc += hw250.conv_weights[i + 1 - m] * hw250.conv_weights[j + 1 - m] * hw250.dt;
                }
                acc
            }
        })
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_cov-exact-vs-hybrid.svg"),
        "Kovarians Volterra Cov(W_tilde^H_s, W_tilde^H_0.5): Eksak vs Hybrid (n=250)",
        "Waktu s",
        "Cov(W_s, W_0.5)",
        &[
            LineSeries {
                label: "Eksak (2F1)",
                x: &s_grid,
                y: &cov_ex_curve,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Hybrid kappa=1",
                x: &s_grid,
                y: &cov_hyb_curve,
                color: PALETTE[1],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 2. Hybrid variance bias vs n
    let ns_f64: Vec<f64> = ns_bias.iter().map(|&n| n as f64).collect();
    let bias_h01: Vec<f64> = hyb_bias_h01
        .iter()
        .map(|r| r.var_bias_theoretical.abs())
        .collect();
    let bias_h05: Vec<f64> = hyb_bias_h05
        .iter()
        .map(|r| r.var_bias_theoretical.abs())
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_hybrid-bias-vs-n.svg"),
        "Bias Teoretis Var(W_tilde^H_1) Skema Hybrid terhadap Jumlah Langkah n",
        "Jumlah Langkah Grid n",
        "|Var_hybrid - 1.0|",
        &[
            LineSeries {
                label: "H = 0.1",
                x: &ns_f64,
                y: &bias_h01,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "H = 0.5 (Brownian)",
                x: &ns_f64,
                y: &bias_h05,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 3. Smile across maturities (T = 0.05, 0.25, 1.0)
    let iv_005: Vec<f64> = smile
        .iter()
        .filter(|pt| (pt.t_horizon - 0.05).abs() < 1e-12)
        .map(|pt| pt.iv)
        .collect();
    let iv_025: Vec<f64> = smile
        .iter()
        .filter(|pt| (pt.t_horizon - 0.25).abs() < 1e-12)
        .map(|pt| pt.iv)
        .collect();
    let iv_100: Vec<f64> = smile
        .iter()
        .filter(|pt| (pt.t_horizon - 1.0).abs() < 1e-12)
        .map(|pt| pt.iv)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_smile-term-structure.svg"),
        "Smile Volatilitas Implisit rough Bergomi (H=0.1, eta=1.9, rho=-0.7)",
        "Log-moneyness k = ln(K/S0)",
        "Volatilitas Implisit",
        &[
            LineSeries {
                label: "T = 0.05",
                x: &smile_ks,
                y: &iv_005,
                color: PALETTE[1],
                dashed: false,
            },
            LineSeries {
                label: "T = 0.25",
                x: &smile_ks,
                y: &iv_025,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "T = 1.00",
                x: &smile_ks,
                y: &iv_100,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 4. Skew vs T log-log
    let log_ts_vec: Vec<f64> = smile_ts.iter().map(|t| t.ln()).collect();
    let log_sk_005: Vec<f64> = skew_study[0].atm_skews.iter().map(|s| s.ln()).collect();
    let log_sk_010: Vec<f64> = skew_study[1].atm_skews.iter().map(|s| s.ln()).collect();
    let log_sk_025: Vec<f64> = skew_study[2].atm_skews.iter().map(|s| s.ln()).collect();
    let log_sk_050: Vec<f64> = skew_study[3].atm_skews.iter().map(|s| s.ln()).collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_skew-loglog.svg"),
        "Skala Skew ATM log psi(T) vs log T untuk H in {0.05, 0.1, 0.25, 0.5}",
        "ln(T)",
        "ln(psi(T))",
        &[
            LineSeries {
                label: "H=0.05",
                x: &log_ts_vec,
                y: &log_sk_005,
                color: PALETTE[1],
                dashed: false,
            },
            LineSeries {
                label: "H=0.10",
                x: &log_ts_vec,
                y: &log_sk_010,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "H=0.25",
                x: &log_ts_vec,
                y: &log_sk_025,
                color: PALETTE[3],
                dashed: false,
            },
            LineSeries {
                label: "H=0.50 (Markov)",
                x: &log_ts_vec,
                y: &log_sk_050,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 5. VRF across strikes at T=0.1 and T=1.0
    let vrf_01: Vec<f64> = pricing_cmp
        .rows
        .iter()
        .filter(|r| (r.t_horizon - 0.1).abs() < 1e-12)
        .map(|r| r.vrf)
        .collect();
    let vrf_10: Vec<f64> = pricing_cmp
        .rows
        .iter()
        .filter(|r| (r.t_horizon - 1.0).abs() < 1e-12)
        .map(|r| r.vrf)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_vrf-conditional-bs.svg"),
        "Faktor Reduksi Variansi (VRF) Estimator Bersyarat Black-Scholes",
        "Log-moneyness k",
        "VRF = Var(MC) / Var(Cond-BS)",
        &[
            LineSeries {
                label: "T = 0.1",
                x: &pricing_ks,
                y: &vrf_01,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "T = 1.0",
                x: &pricing_ks,
                y: &vrf_10,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 6. H-profile RMSE
    let prof_h: Vec<f64> = calib.h_profile.iter().map(|r| r.h).collect();
    let prof_rmse: Vec<f64> = calib.h_profile.iter().map(|r| r.rmse * 100.0).collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_calib-h-profile.svg"),
        "Profil Objektif Kalibrasi RMSE (%) terhadap Hurst H",
        "Parameter Hurst H",
        "IV RMSE (%)",
        &[LineSeries {
            label: "Profile RMSE(H)",
            x: &prof_h,
            y: &prof_rmse,
            color: PALETTE[0],
            dashed: false,
        }],
    )
    .unwrap();

    // 7. Rough vs Markov smile at T=0.1
    let rb_01_ivs: Vec<f64> = markov_cmp.rbergomi_smile.iter().map(|pt| pt.iv).collect();
    let mk_01_ivs: Vec<f64> = markov_cmp.markov_smile.iter().map(|pt| pt.iv).collect();
    svg::write_line_chart(
        &fig_dir.join("day004h_rough-vs-markov-t01.svg"),
        "Perbandingan Smile T=0.1: rough Bergomi (H=0.1) vs Markov (H=0.5)",
        "Log-moneyness k",
        "Volatilitas Implisit",
        &[
            LineSeries {
                label: "rough Bergomi (H=0.1)",
                x: &smile_ks,
                y: &rb_01_ivs,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Markov (H=0.5, ATM matched)",
                x: &smile_ks,
                y: &mk_01_ivs,
                color: PALETTE[1],
                dashed: true,
            },
        ],
    )
    .unwrap();

    let part5_secs = t5.elapsed().as_secs_f64();
    let total_secs = t_total.elapsed().as_secs_f64();

    // ==================== Criteria Evaluation ====================
    let c1_cov_quad = max_cov_quad_rel_err < 1e-8;
    // A priori Bonferroni threshold for m = n*(n+1)/2 = 1275 W-block covariance pairs at two-sided alpha = 0.05:
    // Phi^{-1}(1 - 0.025/1275) = 4.112 -> 4.12
    let z_bonf_c2 = 4.12_f64;
    let c2_chol_ref = chol_ref.max_chol_recon_err < 1e-12 && chol_ref.max_cov_z_score <= z_bonf_c2;
    let c3_mart_anchors = chol_mart.z_exp_w <= 3.0 && mart_anchors.iter().all(|m| m.z_s_t <= 3.0);
    let c4_var_swap = mart_anchors.iter().all(|m| m.z_var_swap <= 3.0);
    let max_z_cond_vs_mc = pricing_cmp
        .rows
        .iter()
        .map(|r| r.z_diff)
        .fold(0.0_f64, f64::max);
    let min_vrf = pricing_cmp
        .rows
        .iter()
        .map(|r| r.vrf)
        .fold(f64::INFINITY, f64::min);
    let max_parity_err = pricing_cmp
        .rows
        .iter()
        .map(|r| r.put_call_parity_err)
        .fold(0.0_f64, f64::max);
    let max_z_ref_hyb = ref_vs_hyb.iter().map(|r| r.z_comb).fold(0.0_f64, f64::max);
    let c5_cond_bs = max_z_cond_vs_mc <= 3.5 && min_vrf > 1.5 && max_parity_err < 1e-12;
    let c6_ref_price = max_z_ref_hyb <= 3.0;
    // Skew scaling: H=0.1 slope within +/- 0.15 of -0.40, and H=0.5 slope near 0 (|slope| < 0.08)
    let c7_skew_scaling =
        (skew_study[1].slope - (-0.4)).abs() < 0.15 && skew_study[3].slope.abs() < 0.08;
    // Calibration identifiability pre-set threshold: |bias_H| < 0.03
    let c8_calib_bias = calib.replication.bias_h.abs() < 0.03;
    let c9_budget = total_secs < 300.0 || quick;

    let crit_str = |ok: bool| if ok { "tercapai" } else { "belum tercapai" };

    let doc = JsonVal::Obj(vec![
        (
            "metadata",
            JsonVal::Obj(vec![
                ("session", "day004h".into()),
                ("date_wita", "2026-11-22".into()),
                ("language", "rust".into()),
                ("quick_mode", quick.into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as usize).into()),
                (
                    "elapsed_seconds",
                    JsonVal::Obj(vec![
                        ("part1_volterra_hybrid", part1_secs.into()),
                        ("part2_martingale_pricing", part2_secs.into()),
                        ("part3_smile_skew", part3_secs.into()),
                        ("part4_calibration", part4_secs.into()),
                        ("part5_figures", part5_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part1_volterra",
            JsonVal::Obj(vec![
                ("max_cov_quad_rel_err", max_cov_quad_rel_err.into()),
                ("cross_cov_sample_05_01_012", cross_cov_sample.into()),
                ("cov_pairs", cov_pairs_json.into()),
                (
                    "cholesky_ref_h01",
                    JsonVal::Obj(vec![
                        ("n", chol_ref.n.into()),
                        ("n_paths", chol_ref.n_paths.into()),
                        ("max_chol_recon_err", chol_ref.max_chol_recon_err.into()),
                        ("max_cov_z_score", chol_ref.max_cov_z_score.into()),
                        ("var_w_terminal_emp", chol_ref.var_w_terminal_emp.into()),
                        ("var_w_terminal_exact", chol_ref.var_w_terminal_exact.into()),
                    ]),
                ),
                (
                    "brownian_limit_h05",
                    JsonVal::Obj(vec![("max_cov_err_vs_min", h05_max_cov_err.into())]),
                ),
                (
                    "hybrid_bias_right_endpoint_h01",
                    hyb_bias_re
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("n", r.n.into()),
                                ("var_exact", r.var_w_exact.into()),
                                ("var_hybrid_th", r.var_w_hybrid_theoretical.into()),
                                ("var_hybrid_mc", r.var_w_hybrid_mc.into()),
                                ("bias_th", r.var_bias_theoretical.into()),
                                ("max_cov_abs_err", r.max_cov_abs_err_theoretical.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "hybrid_bias_h01",
                    hyb_bias_h01
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("n", r.n.into()),
                                ("var_exact", r.var_w_exact.into()),
                                ("var_hybrid_th", r.var_w_hybrid_theoretical.into()),
                                ("var_hybrid_mc", r.var_w_hybrid_mc.into()),
                                ("bias_th", r.var_bias_theoretical.into()),
                                ("max_cov_abs_err", r.max_cov_abs_err_theoretical.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "hybrid_bias_variance_matching_h01",
                    hyb_bias_vm
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("n", r.n.into()),
                                ("var_exact", r.var_w_exact.into()),
                                ("var_hybrid_th", r.var_w_hybrid_theoretical.into()),
                                ("var_hybrid_mc", r.var_w_hybrid_mc.into()),
                                ("bias_th", r.var_bias_theoretical.into()),
                                ("max_cov_abs_err", r.max_cov_abs_err_theoretical.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "hybrid_bias_h05",
                    hyb_bias_h05
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("n", r.n.into()),
                                ("var_exact", r.var_w_exact.into()),
                                ("var_hybrid_th", r.var_w_hybrid_theoretical.into()),
                                ("var_hybrid_mc", r.var_w_hybrid_mc.into()),
                                ("bias_th", r.var_bias_theoretical.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part2_martingale_and_pricing",
            JsonVal::Obj(vec![
                (
                    "cholesky_martingale_check",
                    JsonVal::Obj(vec![
                        ("mean_exp_w", chol_mart.mean_exp_w.into()),
                        ("se_exp_w", chol_mart.se_exp_w.into()),
                        ("z_exp_w", chol_mart.z_exp_w.into()),
                        ("boot_q025", chol_mart.boot_q025.into()),
                        ("boot_q500", chol_mart.boot_q500.into()),
                        ("boot_q975", chol_mart.boot_q975.into()),
                    ]),
                ),
                (
                    "martingale_anchors",
                    mart_anchors
                        .iter()
                        .map(|m| {
                            JsonVal::Obj(vec![
                                ("t", m.t_horizon.into()),
                                ("mean_exp_w", m.mean_exp_w.into()),
                                ("se_exp_w", m.se_exp_w.into()),
                                ("z_exp_w", m.z_exp_w.into()),
                                ("mean_int_v", m.mean_int_v.into()),
                                ("se_int_v", m.se_int_v.into()),
                                ("exact_var_swap", m.exact_var_swap.into()),
                                ("z_var_swap", m.z_var_swap.into()),
                                ("mean_s_t", m.mean_s_t.into()),
                                ("se_s_t", m.se_s_t.into()),
                                ("z_s_t", m.z_s_t.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("mc_secs", pricing_cmp.mc_secs.into()),
                ("cond_secs", pricing_cmp.cond_secs.into()),
                (
                    "call_estimators",
                    pricing_cmp
                        .rows
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("t", r.t_horizon.into()),
                                ("log_k", r.log_k.into()),
                                ("mc_call", r.mc_call.into()),
                                ("mc_se", r.mc_se.into()),
                                ("cond_call", r.cond_call.into()),
                                ("cond_se", r.cond_se.into()),
                                ("z_diff", r.z_diff.into()),
                                ("vrf", r.vrf.into()),
                                ("put_call_parity_err", r.put_call_parity_err.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "ref_vs_hybrid_n50_t1",
                    ref_vs_hyb
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("log_k", r.log_k.into()),
                                ("ref_call", r.ref_call.into()),
                                ("ref_se", r.ref_se.into()),
                                ("hyb_call", r.hyb_call.into()),
                                ("hyb_se", r.hyb_se.into()),
                                ("z_comb", r.z_comb.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part3_smile_and_skew",
            JsonVal::Obj(vec![
                (
                    "atm_skew_scaling",
                    skew_study
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("h", r.h.into()),
                                ("theoretical_slope", r.theoretical_slope.into()),
                                ("measured_slope", r.slope.into()),
                                ("se_slope", r.se_slope.into()),
                                ("r2", r.r2.into()),
                                ("atm_skews", r.atm_skews.clone().into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "markov_compare_t01",
                    JsonVal::Obj(vec![
                        ("eta_markov", markov_cmp.eta_markov.into()),
                        ("rbergomi_atm_iv", markov_cmp.rbergomi_atm_iv.into()),
                        ("markov_atm_iv", markov_cmp.markov_atm_iv.into()),
                        ("rbergomi_skew", markov_cmp.rbergomi_skew.into()),
                        ("markov_skew", markov_cmp.markov_skew.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part4_calibration",
            JsonVal::Obj(vec![
                (
                    "main_fit",
                    JsonVal::Obj(vec![
                        ("h_hat", calib.main_fit.h_hat.into()),
                        ("eta_hat", calib.main_fit.eta_hat.into()),
                        ("rho_hat", calib.main_fit.rho_hat.into()),
                        ("rmse_opt", calib.main_fit.rmse_opt.into()),
                        ("rmse_at_true", calib.main_fit.rmse_at_true.into()),
                        ("evals", calib.main_fit.evals.into()),
                    ]),
                ),
                ("corr_h_eta", calib.corr_h_eta.into()),
                (
                    "h_profile",
                    calib
                        .h_profile
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("h", r.h.into()),
                                ("eta_opt", r.eta_opt.into()),
                                ("rho_opt", r.rho_opt.into()),
                                ("rmse", r.rmse.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "replication_30",
                    JsonVal::Obj(vec![
                        ("n_rep", calib.replication.n_rep.into()),
                        ("bias_h", calib.replication.bias_h.into()),
                        ("sd_h", calib.replication.sd_h.into()),
                        ("bias_eta", calib.replication.bias_eta.into()),
                        ("sd_eta", calib.replication.sd_eta.into()),
                        ("bias_rho", calib.replication.bias_rho.into()),
                        ("sd_rho", calib.replication.sd_rho.into()),
                    ]),
                ),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            JsonVal::Obj(vec![
                (
                    "c1_volterra_cov_quad_within_1e8",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c1_cov_quad).into()),
                        ("max_rel_err", max_cov_quad_rel_err.into()),
                        ("threshold", 1e-8_f64.into()),
                    ]),
                ),
                (
                    "c2_cholesky_exact_cov_bonferroni",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c2_chol_ref).into()),
                        ("max_chol_recon_err", chol_ref.max_chol_recon_err.into()),
                        ("max_cov_z_score", chol_ref.max_cov_z_score.into()),
                        ("bonferroni_z_threshold", z_bonf_c2.into()),
                    ]),
                ),
                (
                    "c3_martingales_within_3se",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c3_mart_anchors).into()),
                        ("chol_z_exp_w", chol_mart.z_exp_w.into()),
                        (
                            "hybrid_max_z_exp_w",
                            mart_anchors
                                .iter()
                                .map(|m| m.z_exp_w)
                                .fold(0.0_f64, f64::max)
                                .into(),
                        ),
                        (
                            "max_z_s_t",
                            mart_anchors
                                .iter()
                                .map(|m| m.z_s_t)
                                .fold(0.0_f64, f64::max)
                                .into(),
                        ),
                    ]),
                ),
                (
                    "c4_variance_swap_within_3se",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c4_var_swap).into()),
                        (
                            "max_z_var_swap",
                            mart_anchors
                                .iter()
                                .map(|m| m.z_var_swap)
                                .fold(0.0_f64, f64::max)
                                .into(),
                        ),
                    ]),
                ),
                (
                    "c5_cond_bs_matches_mc_and_vrf",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c5_cond_bs).into()),
                        ("max_z_diff", max_z_cond_vs_mc.into()),
                        ("min_vrf", min_vrf.into()),
                        ("max_parity_err", max_parity_err.into()),
                    ]),
                ),
                (
                    "c6_hybrid_matches_cholesky_price_3se",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c6_ref_price).into()),
                        ("max_z_comb", max_z_ref_hyb.into()),
                    ]),
                ),
                (
                    "c7_atm_skew_power_law_and_markov_zero",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c7_skew_scaling).into()),
                        ("slope_h01", skew_study[1].slope.into()),
                        ("se_slope_h01", skew_study[1].se_slope.into()),
                        ("slope_h05", skew_study[3].slope.into()),
                        ("se_slope_h05", skew_study[3].se_slope.into()),
                    ]),
                ),
                (
                    "c8_calib_bias_h_under_003",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c8_calib_bias).into()),
                        ("bias_h", calib.replication.bias_h.into()),
                        ("sd_h", calib.replication.sd_h.into()),
                        ("threshold_abs_bias_h", 0.03_f64.into()),
                    ]),
                ),
                (
                    "c9_runtime_under_300s",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c9_budget).into()),
                        ("total_seconds", total_secs.into()),
                        ("budget_seconds", 300.0_f64.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &doc).expect("failed to write results.json");
    println!(
        "day004h selesai dalam {:.2}s | workers={n_workers} | H_hat={:.4} (bias={:.4}) | slope(H=0.1)={:.4} | c1..c9={}/{}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        calib.main_fit.h_hat,
        calib.replication.bias_h,
        skew_study[1].slope,
        c1_cov_quad,
        c2_chol_ref,
        c3_mart_anchors,
        c4_var_swap,
        c5_cond_bs,
        c6_ref_price,
        c7_skew_scaling,
        c8_calib_bias,
        c9_budget
    );
}
