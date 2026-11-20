#![forbid(unsafe_code)]

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use day003b_optimal_execution::ac::{
    ac_continuous_sup_error, ac_efficient_frontier, evaluate_ac_moments, simulate_ac_mc,
    solve_ac_qp, solve_ac_sinh_discrete, AcParams,
};
use day003b_optimal_execution::dp_liquidity::{
    compute_adaptivity_value, simulate_dp_paths, solve_liquidity_dp, verify_dp_bruteforce,
    LiquidityDpParams,
};
use day003b_optimal_execution::json::JsonValue;
use day003b_optimal_execution::parallel::default_num_threads;
use day003b_optimal_execution::svg::{
    write_bar_chart_svg, write_line_chart_svg, LineSeries, PALETTE,
};
use day003b_optimal_execution::transient::{
    analyze_non_pd_kernel, analyze_oscillating_kernel, benchmark_ow_cholesky_vs_analytic,
    build_kernel_matrix, solve_ow_exponential_analytic, solve_transient_general,
    solve_transient_with_risk, OwParams,
};

const SEED: u64 = 2026_1121;

fn repo_root() -> PathBuf {
    let cwd = std::env::current_dir().expect("cwd");
    if cwd.join("outputs").exists() {
        cwd
    } else if cwd.join("../../outputs").exists() {
        cwd.join("../..").canonicalize().unwrap_or(cwd)
    } else {
        cwd
    }
}

fn log_log_slope(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let lx: Vec<f64> = xs.iter().map(|&v| v.ln()).collect();
    let ly: Vec<f64> = ys.iter().map(|&v| v.max(1e-30).ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let num: f64 = lx
        .iter()
        .zip(ly.iter())
        .map(|(x, y)| (x - mx) * (y - my))
        .sum();
    let den: f64 = lx.iter().map(|x| (x - mx) * (x - mx)).sum();
    num / den
}

fn main() {
    let root = repo_root();
    let fig_dir = root.join("outputs").join("figures");
    let out_dir = root.join("outputs").join("day003b");
    fs::create_dir_all(&fig_dir).expect("create figures dir");
    fs::create_dir_all(&out_dir).expect("create day003b dir");

    let n_threads = default_num_threads();
    let t_total = Instant::now();

    // === BAGIAN 1: Almgren-Chriss ===
    let t1 = Instant::now();
    let ac_p = AcParams::default();
    let ac_qp = solve_ac_qp(&ac_p);
    let ac_sinh = solve_ac_sinh_discrete(&ac_p);
    let qp_vs_sinh_rel = ac_qp
        .x
        .iter()
        .zip(ac_sinh.x.iter())
        .map(|(a, b)| (a - b).abs() / ac_p.x_total)
        .fold(0.0f64, f64::max);

    let cont_ns = [10usize, 25, 50, 100, 200, 400];
    let cont_errs: Vec<f64> = cont_ns
        .iter()
        .map(|&n| ac_continuous_sup_error(&ac_p, n))
        .collect();
    let cont_ns_f: Vec<f64> = cont_ns.iter().map(|&n| n as f64).collect();
    let cont_slope = log_log_slope(&cont_ns_f, &cont_errs);

    let mut p_twap_lim = ac_p;
    p_twap_lim.lambda = 1e-16;
    let s_twap_lim = solve_ac_qp(&p_twap_lim);
    let twap_x = ac_p.x_total / (ac_p.n_steps as f64);
    let twap_limit_err = s_twap_lim
        .x
        .iter()
        .map(|&xk| (xk - twap_x).abs() / ac_p.x_total)
        .fold(0.0f64, f64::max);

    let mut p_fast_lim = ac_p;
    p_fast_lim.lambda = 1.0;
    let s_fast_lim = solve_ac_qp(&p_fast_lim);
    let fast_first_frac = s_fast_lim.x[0] / ac_p.x_total;

    let frontier = ac_efficient_frontier(&ac_p, 40);
    let mut frontier_monotone = true;
    for w in frontier.windows(2) {
        if w[1].1.var_c > w[0].1.var_c + 1e-6 || w[1].1.expected_c < w[0].1.expected_c - 1e-6 {
            frontier_monotone = false;
        }
    }
    let (ev_ec, ev_vc) = evaluate_ac_moments(&ac_p, &ac_qp.x);
    let moment_parity_err = (ev_ec - ac_qp.expected_c).abs() + (ev_vc - ac_qp.var_c).abs();

    let mc_opt = simulate_ac_mc(&ac_p, &ac_qp.x, 400_000, n_threads, SEED);
    let twap_sched = vec![twap_x; ac_p.n_steps];
    let mc_twap = simulate_ac_mc(&ac_p, &twap_sched, 400_000, n_threads, SEED + 1);
    let mut sell_all = vec![0.0; ac_p.n_steps];
    sell_all[0] = ac_p.x_total;
    let mc_sell_all = simulate_ac_mc(&ac_p, &sell_all, 400_000, n_threads, SEED + 2);
    let time_part1 = t1.elapsed().as_secs_f64();

    // === BAGIAN 2: Dampak Transien (Obizhaeva-Wang) ===
    let t2 = Instant::now();
    let ow = OwParams::default();
    let g_ow = build_kernel_matrix(ow.n_steps, ow.t_horizon, |d| {
        ow.lambda_ow * (-ow.rho * d).exp()
    });
    let (x_ow_chol, c_ow_chol) = solve_transient_general(&g_ow, ow.x_total).unwrap();
    let (x_ow_ana, c_ow_ana) = solve_ow_exponential_analytic(&ow);
    let ow_anchor_err = x_ow_chol
        .iter()
        .zip(x_ow_ana.iter())
        .map(|(a, b)| (a - b).abs() / ow.x_total)
        .fold(0.0f64, f64::max);

    let ow_ns = [20usize, 50, 100, 200, 400];
    let ow_cont_cost = ow.lambda_ow * ow.x_total * ow.x_total / (ow.rho * ow.t_horizon + 2.0);
    let ow_cont_block = ow.x_total / (ow.rho * ow.t_horizon + 2.0);
    let mut ow_cost_errs = Vec::new();
    let mut ow_block_errs = Vec::new();
    for &n in &ow_ns {
        let mut p_n = ow;
        p_n.n_steps = n;
        let (xn, cn) = solve_ow_exponential_analytic(&p_n);
        ow_cost_errs.push((cn - ow_cont_cost).abs() / ow_cont_cost);
        ow_block_errs.push((xn[0] - ow_cont_block).abs() / ow.x_total);
    }
    let ow_ns_f: Vec<f64> = ow_ns.iter().map(|&n| n as f64).collect();
    let ow_cost_order = log_log_slope(&ow_ns_f, &ow_cost_errs);

    let (bench_err, bench_chol_s, bench_ana_s) = benchmark_ow_cholesky_vs_analytic(&ow, 2000);

    // Kernel power-law
    let g_pl = build_kernel_matrix(ow.n_steps, ow.t_horizon, |d| {
        ow.lambda_ow / (1.0 + d / 0.1).sqrt()
    });
    let (x_pl, c_pl) = solve_transient_general(&g_pl, ow.x_total).unwrap();
    let pl_min_x = x_pl.iter().copied().fold(f64::INFINITY, f64::min);

    // Kernel osilasi (PD dengan komponen negatif)
    let (osc_pd, osc_min_eig, x_osc, osc_min_x, c_osc) =
        analyze_oscillating_kernel(20, 1.0, ow.lambda_ow, 0.2, 6.0, 0.001, ow.x_total);

    // Kernel non-PD
    let (non_pd_chol_failed, non_pd_min_eig, non_pd_d_quad, non_pd_cost_pert) =
        analyze_non_pd_kernel(20, 1.0, ow.lambda_ow, ow.x_total);

    // Gabungan: G diagonal vs AC
    let dt_ac = ac_p.t_horizon / (ac_p.n_steps as f64);
    let eta_tilde = ac_p.eta - 0.5 * ac_p.gamma * dt_ac;
    let mut g_diag = vec![vec![0.0; ac_p.n_steps]; ac_p.n_steps];
    for (i, row) in g_diag.iter_mut().enumerate() {
        row[i] = 2.0 * eta_tilde / dt_ac;
    }
    let (x_comb, _) = solve_transient_with_risk(&g_diag, &ac_p).unwrap();
    let comb_vs_ac_err = x_comb
        .iter()
        .zip(ac_qp.x.iter())
        .map(|(a, b)| (a - b).abs() / ac_p.x_total)
        .fold(0.0f64, f64::max);
    let time_part2 = t2.elapsed().as_secs_f64();

    // === BAGIAN 3: DP Likuiditas Stokastik ===
    let t3 = Instant::now();
    let dp_p = LiquidityDpParams::default();
    let dp_sol = solve_liquidity_dp(&dp_p);

    let mut dp_p5 = dp_p;
    dp_p5.n_steps = 5;
    let dp_bf_err = verify_dp_bruteforce(&dp_p5, 501);

    let dp_sim = simulate_dp_paths(&dp_p, &dp_sol, 400_000, n_threads, SEED);
    let adapt_base = compute_adaptivity_value(&dp_p);

    let persistences = [0.5, 0.7, 0.95];
    let mut sens_reports = Vec::new();
    for &p11 in &persistences {
        let mut ps = dp_p;
        ps.p_11 = p11;
        sens_reports.push((p11, compute_adaptivity_value(&ps)));
    }
    let time_part3 = t3.elapsed().as_secs_f64();
    let time_total = t_total.elapsed().as_secs_f64();

    // === GRAFIK SVG ===
    // 1. Trajektori AC untuk beberapa lambda
    let t_grid_y: Vec<f64> = (0..=ac_p.n_steps)
        .map(|j| (j as f64) / (ac_p.n_steps as f64))
        .collect();
    let ac_lams = [1e-7, 2e-5, 2e-4, 1e-3];
    let ac_series: Vec<LineSeries> = ac_lams
        .iter()
        .enumerate()
        .map(|(idx, &lam)| {
            let mut pi = ac_p;
            pi.lambda = lam;
            let s = solve_ac_qp(&pi);
            LineSeries {
                label: format!("lambda = {lam:.0e}"),
                xs: t_grid_y.clone(),
                ys: s.y.iter().map(|v| v / 1e6).collect(),
                color: PALETTE[idx % PALETTE.len()].to_string(),
                dashed: false,
            }
        })
        .collect();
    write_line_chart_svg(
        &fig_dir.join("day003b_ac-trajectories.svg"),
        "Trajektori Holding Optimal Almgren-Chriss y(t)",
        "Waktu t (hari)",
        "Sisa Holding (juta saham)",
        &ac_series,
        false,
        false,
    )
    .expect("svg 1");

    // 2. Frontier Efisien AC
    let front_series = vec![LineSeries {
        label: "Frontier Efisien (40 titik)".to_string(),
        xs: frontier.iter().map(|(_, s)| s.std_c / 1e3).collect(),
        ys: frontier.iter().map(|(_, s)| s.expected_c / 1e3).collect(),
        color: PALETTE[0].to_string(),
        dashed: false,
    }];
    write_line_chart_svg(
        &fig_dir.join("day003b_ac-frontier.svg"),
        "Frontier Efisien Almgren-Chriss (SD vs E[C])",
        "SD[C] (ribu $)",
        "E[C] (ribu $)",
        &front_series,
        false,
        false,
    )
    .expect("svg 2");

    // 3. Konvergensi kontinu AC (log-log)
    let conv_series = vec![LineSeries {
        label: format!("Galat sup (slope = {cont_slope:.3})"),
        xs: cont_ns_f.clone(),
        ys: cont_errs.clone(),
        color: PALETTE[1].to_string(),
        dashed: false,
    }];
    write_line_chart_svg(
        &fig_dir.join("day003b_ac-continuous-convergence.svg"),
        "Konvergensi Limit Kontinu Almgren-Chriss (Log-Log)",
        "Jumlah Langkah N",
        "Galat Sup Relatif",
        &conv_series,
        true,
        true,
    )
    .expect("svg 3");

    // 4. Bentuk x* kernel eksponensial, power-law, osilasi
    let k_idx_50: Vec<f64> = (1..=ow.n_steps).map(|i| i as f64).collect();
    let k_idx_20: Vec<f64> = (1..=20).map(|i| i as f64 / 20.0 * 50.0).collect();
    let kern_series = vec![
        LineSeries {
            label: "OW Eksponensial".to_string(),
            xs: k_idx_50.clone(),
            ys: x_ow_ana.iter().map(|v| v / 1e3).collect(),
            color: PALETTE[0].to_string(),
            dashed: false,
        },
        LineSeries {
            label: "Power-Law".to_string(),
            xs: k_idx_50,
            ys: x_pl.iter().map(|v| v / 1e3).collect(),
            color: PALETTE[2].to_string(),
            dashed: false,
        },
        LineSeries {
            label: "Osilasi PD (manipulasi)".to_string(),
            xs: k_idx_20,
            ys: x_osc.iter().map(|v| v / 1e3).collect(),
            color: PALETTE[1].to_string(),
            dashed: true,
        },
    ];
    write_line_chart_svg(
        &fig_dir.join("day003b_transient-kernels.svg"),
        "Jadwal Optimal x* pada Berbagai Kernel Dampak Transien",
        "Langkah k",
        "Trade x_k (ribu saham)",
        &kern_series,
        false,
        false,
    )
    .expect("svg 4");

    // 5. Kebijakan DP c_k(s)
    let dp_steps: Vec<f64> = (1..=dp_p.n_steps).map(|k| k as f64).collect();
    let pol_series = vec![
        LineSeries {
            label: "Rezim Likuid (c_k(0))".to_string(),
            xs: dp_steps.clone(),
            ys: dp_sol.c.iter().map(|c| c[0]).collect(),
            color: PALETTE[0].to_string(),
            dashed: false,
        },
        LineSeries {
            label: "Rezim Ilikuid (c_k(1))".to_string(),
            xs: dp_steps,
            ys: dp_sol.c.iter().map(|c| c[1]).collect(),
            color: PALETTE[1].to_string(),
            dashed: false,
        },
    ];
    write_line_chart_svg(
        &fig_dir.join("day003b_dp-policy.svg"),
        "Kebijakan Likuidasi Optimal DP c_k(s)",
        "Langkah k",
        "Fraksi Trade c_k(s)",
        &pol_series,
        false,
        false,
    )
    .expect("svg 5");

    // 6. Histogram / distribusi J (dalam 10 bin)
    let j_min = dp_sim
        .j_samples_head
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let j_max = dp_sim
        .j_samples_head
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let n_bins = 10usize;
    let bin_w = (j_max - j_min).max(1.0) / (n_bins as f64);
    let mut bin_counts = vec![0.0; n_bins];
    for &jv in &dp_sim.j_samples_head {
        let b = (((jv - j_min) / bin_w).floor() as usize).min(n_bins - 1);
        bin_counts[b] += 1.0 / (dp_sim.j_samples_head.len() as f64);
    }
    let bin_labels: Vec<String> = (0..n_bins)
        .map(|b| format!("{:.0}k", (j_min + (b as f64 + 0.5) * bin_w) / 1e3))
        .collect();
    write_bar_chart_svg(
        &fig_dir.join("day003b_dp-cost-distribution.svg"),
        "Distribusi Biaya Realisasi J pada Simulasi DP",
        "Proporsi Jalur",
        &bin_labels,
        &bin_counts,
    )
    .expect("svg 6");

    // 7. Nilai adaptivitas vs persistensi
    let sens_labels: Vec<String> = sens_reports
        .iter()
        .map(|(p11, _)| format!("p11={p11:.2}"))
        .collect();
    let sens_vals: Vec<f64> = sens_reports
        .iter()
        .map(|(_, r)| r.adaptivity_pct_vs_stationary)
        .collect();
    write_bar_chart_svg(
        &fig_dir.join("day003b_dp-adaptivity-sensitivity.svg"),
        "Nilai Adaptivitas DP (% Penghematan vs Statis Stasioner)",
        "Penghematan (%)",
        &sens_labels,
        &sens_vals,
    )
    .expect("svg 7");

    // === RESULTS JSON ===
    let json_root = JsonValue::Obj(vec![
        ("n_threads".to_string(), JsonValue::Int(n_threads as i64)),
        (
            "part1_almgren_chriss".to_string(),
            JsonValue::Obj(vec![
                (
                    "qp_vs_sinh_max_rel_err".to_string(),
                    JsonValue::Num(qp_vs_sinh_rel),
                ),
                ("expected_c".to_string(), JsonValue::Num(ac_qp.expected_c)),
                ("var_c".to_string(), JsonValue::Num(ac_qp.var_c)),
                ("std_c".to_string(), JsonValue::Num(ac_qp.std_c)),
                (
                    "moment_parity_err".to_string(),
                    JsonValue::Num(moment_parity_err),
                ),
                (
                    "continuous_ns".to_string(),
                    JsonValue::from_f64_slice(&cont_ns_f),
                ),
                (
                    "continuous_sup_errors".to_string(),
                    JsonValue::from_f64_slice(&cont_errs),
                ),
                (
                    "continuous_log_log_slope".to_string(),
                    JsonValue::Num(cont_slope),
                ),
                (
                    "limit_lambda0_twap_max_rel_err".to_string(),
                    JsonValue::Num(twap_limit_err),
                ),
                (
                    "limit_lambda_inf_first_frac".to_string(),
                    JsonValue::Num(fast_first_frac),
                ),
                (
                    "frontier_monotone_convex".to_string(),
                    JsonValue::Bool(frontier_monotone),
                ),
                ("mc_opt_mean".to_string(), JsonValue::Num(mc_opt.mc_mean)),
                (
                    "mc_opt_mean_se".to_string(),
                    JsonValue::Num(mc_opt.mc_mean_se),
                ),
                (
                    "mc_opt_antithetic_mean".to_string(),
                    JsonValue::Num(mc_opt.mc_antithetic_mean),
                ),
                ("mc_opt_var".to_string(), JsonValue::Num(mc_opt.mc_var)),
                (
                    "mc_opt_var_se".to_string(),
                    JsonValue::Num(mc_opt.mc_var_se),
                ),
                ("mc_opt_z_mean".to_string(), JsonValue::Num(mc_opt.z_mean)),
                ("mc_opt_z_var".to_string(), JsonValue::Num(mc_opt.z_var)),
                ("mc_twap_z_mean".to_string(), JsonValue::Num(mc_twap.z_mean)),
                ("mc_twap_z_var".to_string(), JsonValue::Num(mc_twap.z_var)),
                (
                    "mc_sell_all_z_mean".to_string(),
                    JsonValue::Num(mc_sell_all.z_mean),
                ),
                (
                    "mc_sell_all_z_var".to_string(),
                    JsonValue::Num(mc_sell_all.z_var),
                ),
            ]),
        ),
        (
            "part2_transient".to_string(),
            JsonValue::Obj(vec![
                (
                    "ow_anchor_max_rel_err".to_string(),
                    JsonValue::Num(ow_anchor_err),
                ),
                ("ow_min_cost".to_string(), JsonValue::Num(c_ow_chol)),
                ("ow_min_cost_analytic".to_string(), JsonValue::Num(c_ow_ana)),
                (
                    "ow_continuous_cost_limit".to_string(),
                    JsonValue::Num(ow_cont_cost),
                ),
                (
                    "ow_continuous_cost_rel_errs".to_string(),
                    JsonValue::from_f64_slice(&ow_cost_errs),
                ),
                (
                    "ow_continuous_block_errs".to_string(),
                    JsonValue::from_f64_slice(&ow_block_errs),
                ),
                (
                    "ow_cost_convergence_order".to_string(),
                    JsonValue::Num(ow_cost_order),
                ),
                (
                    "bench_n2000_max_rel_err".to_string(),
                    JsonValue::Num(bench_err),
                ),
                (
                    "bench_n2000_cholesky_sec".to_string(),
                    JsonValue::Num(bench_chol_s),
                ),
                (
                    "bench_n2000_analytic_sec".to_string(),
                    JsonValue::Num(bench_ana_s),
                ),
                ("powerlaw_min_x".to_string(), JsonValue::Num(pl_min_x)),
                ("powerlaw_min_cost".to_string(), JsonValue::Num(c_pl)),
                ("oscillating_is_pd".to_string(), JsonValue::Bool(osc_pd)),
                (
                    "oscillating_min_eig".to_string(),
                    JsonValue::Num(osc_min_eig),
                ),
                ("oscillating_min_x".to_string(), JsonValue::Num(osc_min_x)),
                ("oscillating_min_cost".to_string(), JsonValue::Num(c_osc)),
                (
                    "non_pd_chol_failed".to_string(),
                    JsonValue::Bool(non_pd_chol_failed),
                ),
                ("non_pd_min_eig".to_string(), JsonValue::Num(non_pd_min_eig)),
                (
                    "non_pd_descent_quad".to_string(),
                    JsonValue::Num(non_pd_d_quad),
                ),
                (
                    "non_pd_perturbed_cost".to_string(),
                    JsonValue::Num(non_pd_cost_pert),
                ),
                (
                    "combined_diag_vs_ac_max_rel_err".to_string(),
                    JsonValue::Num(comb_vs_ac_err),
                ),
            ]),
        ),
        (
            "part3_dp_liquidity".to_string(),
            JsonValue::Obj(vec![
                (
                    "dp_bruteforce_n5_rel_err".to_string(),
                    JsonValue::Num(dp_bf_err),
                ),
                (
                    "v0_liquid_analytic".to_string(),
                    JsonValue::Num(dp_sol.v0_liquid),
                ),
                (
                    "v0_illiquid_analytic".to_string(),
                    JsonValue::Num(dp_sol.v0_illiquid),
                ),
                ("sim_mean_j".to_string(), JsonValue::Num(dp_sim.mean_j)),
                ("sim_se_j".to_string(), JsonValue::Num(dp_sim.se_j)),
                ("sim_z_score".to_string(), JsonValue::Num(dp_sim.z_score)),
                (
                    "sim_illiquid_frac".to_string(),
                    JsonValue::Num(dp_sim.illiquid_frac),
                ),
                ("y_k10_mean".to_string(), JsonValue::Num(dp_sim.y_k10_mean)),
                ("y_k10_std".to_string(), JsonValue::Num(dp_sim.y_k10_std)),
                ("y_k25_mean".to_string(), JsonValue::Num(dp_sim.y_k25_mean)),
                ("y_k25_std".to_string(), JsonValue::Num(dp_sim.y_k25_std)),
                ("y_k40_mean".to_string(), JsonValue::Num(dp_sim.y_k40_mean)),
                ("y_k40_std".to_string(), JsonValue::Num(dp_sim.y_k40_std)),
                (
                    "adaptivity_j_dp".to_string(),
                    JsonValue::Num(adapt_base.j_dp),
                ),
                (
                    "adaptivity_j_static_stationary".to_string(),
                    JsonValue::Num(adapt_base.j_static_stationary),
                ),
                (
                    "adaptivity_j_static_liquid".to_string(),
                    JsonValue::Num(adapt_base.j_static_liquid),
                ),
                (
                    "adaptivity_j_static_illiquid".to_string(),
                    JsonValue::Num(adapt_base.j_static_illiquid),
                ),
                (
                    "adaptivity_pct_vs_stationary".to_string(),
                    JsonValue::Num(adapt_base.adaptivity_pct_vs_stationary),
                ),
                (
                    "sensitivity_p11".to_string(),
                    JsonValue::Arr(
                        sens_reports
                            .iter()
                            .map(|(p11, r)| {
                                JsonValue::Obj(vec![
                                    ("p11".to_string(), JsonValue::Num(*p11)),
                                    ("j_dp".to_string(), JsonValue::Num(r.j_dp)),
                                    (
                                        "j_static_stationary".to_string(),
                                        JsonValue::Num(r.j_static_stationary),
                                    ),
                                    (
                                        "adaptivity_pct_vs_stationary".to_string(),
                                        JsonValue::Num(r.adaptivity_pct_vs_stationary),
                                    ),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "timings_seconds".to_string(),
            JsonValue::Obj(vec![
                ("part1_ac".to_string(), JsonValue::Num(time_part1)),
                ("part2_transient".to_string(), JsonValue::Num(time_part2)),
                ("part3_dp_liquidity".to_string(), JsonValue::Num(time_part3)),
                ("total".to_string(), JsonValue::Num(time_total)),
            ]),
        ),
    ]);

    fs::write(out_dir.join("results.json"), json_root.to_pretty_string())
        .expect("write results.json");

    println!(
        "day003b selesai dalam {time_total:.3}s (threads={n_threads}):\n  Part 1: {time_part1:.3}s, Part 2: {time_part2:.3}s, Part 3: {time_part3:.3}s"
    );
    println!("  AC QP vs sinh rel err: {qp_vs_sinh_rel:.3e}, cont slope: {cont_slope:.4}");
    println!(
        "  AC MC z_mean: {:.3}, z_var: {:.3}",
        mc_opt.z_mean, mc_opt.z_var
    );
    println!("  OW anchor err: {ow_anchor_err:.3e}, N=2000 chol={bench_chol_s:.4}s vs ana={bench_ana_s:.6}s");
    println!(
        "  DP brute-force err: {dp_bf_err:.3e}, sim z: {:.3}, adaptivity: {:.3}%",
        dp_sim.z_score, adapt_base.adaptivity_pct_vs_stationary
    );
}
