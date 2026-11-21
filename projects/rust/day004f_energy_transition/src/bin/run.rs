//! Unified CLI runner for day004f: runs Panels A-F, verifies all anchors & parity,
//! writes outputs/day004f/results.json and 7 SVG figures in outputs/figures/.

use std::path::{Path, PathBuf};
use std::time::Instant;

use day004f_energy_transition::afns::AfnsParams;
use day004f_energy_transition::dp::{trace_example_day, PriceGrid};
use day004f_energy_transition::integration;
use day004f_energy_transition::jsonout::{write_json_file, JsonVal};
use day004f_energy_transition::ns::NsParams;
use day004f_energy_transition::panel_a;
use day004f_energy_transition::panel_b;
use day004f_energy_transition::panel_c;
use day004f_energy_transition::panel_d;
use day004f_energy_transition::panel_e;
use day004f_energy_transition::parallel;
use day004f_energy_transition::risk::STANDARD_KRD_NODES;
use day004f_energy_transition::rng;
use day004f_energy_transition::svg::{self, LineSeries, PALETTE};

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
        std::env::temp_dir().join("day004f_quick_results.json")
    } else {
        root.join("outputs").join("day004f").join("results.json")
    };
    let fig_dir = if quick {
        std::env::temp_dir().join("day004f_quick_figures")
    } else {
        root.join("outputs").join("figures")
    };

    let t_total = Instant::now();
    println!("[day004f] Running Panels A-F (quick={quick}, workers={n_workers})...");

    let pa = panel_a::run_panel_a(n_workers, quick);
    println!(
        "[day004f] Panel A (day004a) done in {:.2}s",
        pa.elapsed_seconds
    );

    let pb = panel_b::run_panel_b(n_workers, quick);
    println!(
        "[day004f] Panel B (day004b) done in {:.2}s",
        pb.elapsed_seconds
    );

    let pc = panel_c::run_panel_c(n_workers, quick);
    println!(
        "[day004f] Panel C (day004c) done in {:.2}s",
        pc.elapsed_seconds
    );

    let pd = panel_d::run_panel_d(n_workers, quick);
    println!(
        "[day004f] Panel D (day004d) done in {:.2}s",
        pd.elapsed_seconds
    );

    let pe = panel_e::run_panel_e(n_workers, quick);
    println!(
        "[day004f] Panel E (day004e) done in {:.2}s",
        pe.elapsed_seconds
    );

    let pf = integration::run_integration(
        &pe.res_dp33.daily_profits,
        &pb.baseline.replications,
        &pc.jackknife_20y.full_fits,
        &pd.outer_fits,
        quick,
    );
    println!(
        "[day004f] Panel F (Integration) done in {:.2}s",
        pf.elapsed_seconds
    );

    // Acceptance checks
    let c_anchors = pa.max_integral_err_ns < 1e-12
        && pa.max_bootstrap_df_err < 1e-12
        && pa.bond_risk.mod_duration_ytm_rel_err < 1e-6
        && pa.bond_risk.convexity_ytm_rel_err < 1e-6
        && pa.krd.sum_vs_par_closed_err < 1e-10
        && pa.afns_riccati_max_a_err < 1e-10
        && pa.afns_phi_max_err < 1e-12
        && pa.afns_q_vl_max_err < 1e-12
        && pa.anchor_t5.log_lik_abs_diff < 1e-9
        && pa.anchor_t5.max_smooth_mean_abs_diff < 1e-9
        && pa.anchor_t5.max_smooth_cov_abs_diff < 1e-9;

    let c_pca = pa.pca_acceptance.aggregate_passed;
    let c_jk_alg = pc.algebraic_c_over_h_max_err < 1e-12;
    let c_power = pe.stationary_residual < 1e-13 && pe.jump_verify_ok && pe.step_verify_ok;
    let c_battery_dp = pe.bf_abs_diff < 1e-11
        && pe.policy_17.max_bellman_residual < 1e-11
        && pe.policy_33.max_bellman_residual < 1e-11
        && pe.policy_exp38.max_bellman_residual < 1e-11
        && pe.res_dp33.soc_violation_count == 0
        && pe.res_dp33.pf_dom_violations == 0
        && pe.res_dp17.pf_dom_violations == 0
        && pe.res_threshold.pf_dom_violations == 0;
    let c_integration = pf.duration_rel_err < 1e-6
        && pf.convexity_rel_err < 1e-6
        && pf.krd.sum_vs_par_closed_err < 1e-10
        && pf.hedge.residual_node1.abs() < 1e-9
        && pf.hedge.residual_node2.abs() < 1e-9;
    let total_seconds = t_total.elapsed().as_secs_f64();
    let c_budget = total_seconds < 600.0 || quick;
    let all_technical_passed =
        c_anchors && c_pca && c_jk_alg && c_power && c_battery_dp && c_integration && c_budget;

    write_all_figures(&fig_dir, &pa, &pb, &pc, &pd, &pe, &pf).expect("failed to write SVG figures");

    let json_doc = build_results_json(
        quick,
        n_workers,
        total_seconds,
        &pa,
        &pb,
        &pc,
        &pd,
        &pe,
        &pf,
        c_anchors,
        c_pca,
        c_jk_alg,
        c_power,
        c_battery_dp,
        c_integration,
        c_budget,
        all_technical_passed,
    );
    write_json_file(&out_json, &json_doc).expect("failed to write results.json");

    println!(
        "day004f selesai dalam {:.2}s | workers={n_workers} | all_technical_passed={all_technical_passed} | NPV_base={:.2} EUR",
        total_seconds, pf.base_npv_pre_capex_eur
    );
}

fn write_all_figures(
    fig_dir: &Path,
    _pa: &panel_a::PanelAResult,
    pb: &panel_b::PanelBResult,
    pc: &panel_c::PanelCResult,
    pd: &panel_d::PanelDResult,
    pe: &panel_e::PanelEResult,
    pf: &integration::IntegrationResult,
) -> std::io::Result<()> {
    let taus: Vec<f64> = (1..=60).map(|i| i as f64 * 0.5).collect();
    let ns_p = NsParams::default_ns();
    let afns_p = AfnsParams::default_params();
    let y_ns: Vec<f64> = taus.iter().map(|&t| ns_p.yield_at(t) * 100.0).collect();
    let f_ns: Vec<f64> = taus.iter().map(|&t| ns_p.forward(t) * 100.0).collect();
    let y_afns: Vec<f64> = taus
        .iter()
        .map(|&t| afns_p.yield_at(afns_p.theta_p, t) * 100.0)
        .collect();
    let adj_bps: Vec<f64> = taus
        .iter()
        .map(|&t| afns_p.yield_adjustment(t).0 * 10_000.0)
        .collect();

    // 1. NS/NSS & AFNS
    svg::write_line_chart(
        &fig_dir.join("day004f_ns-afns-curves-and-adjustment.svg"),
        "Kurva Yield NS vs AFNS (%) dan Koreksi Konveksitas",
        "Tenor (Tahun)",
        "Yield (%) / Adj (bps/100)",
        &[
            LineSeries {
                label: "NS Yield (%)",
                x: &taus,
                y: &y_ns,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "NS Forward (%)",
                x: &taus,
                y: &f_ns,
                color: PALETTE[1],
                dashed: true,
            },
            LineSeries {
                label: "AFNS Yield (%)",
                x: &taus,
                y: &y_afns,
                color: PALETTE[2],
                dashed: false,
            },
            LineSeries {
                label: "AFNS -A/tau (bps/10)",
                x: &taus,
                y: &adj_bps.iter().map(|x| x * 0.1).collect::<Vec<_>>(),
                color: PALETTE[3],
                dashed: true,
            },
        ],
    )?;

    // 2. Horizon diagnostics (kp11, kp22, kp33 bias vs H)
    let hs: Vec<f64> = pb.kp_trends[0].horizons.to_vec();
    let b1: Vec<f64> = pb.kp_trends[0].biases.to_vec();
    let b2: Vec<f64> = pb.kp_trends[1].biases.to_vec();
    let b3: Vec<f64> = pb.kp_trends[2].biases.to_vec();
    svg::write_line_chart(
        &fig_dir.join("day004f_kalman-pca-horizon-diagnostics.svg"),
        "Peluruhan Bias Finite-Sample KP_ii terhadap Horizon H (10y-80y)",
        "Horizon H (Tahun)",
        "Bias MLE (thn^-1)",
        &[
            LineSeries {
                label: "kp11 bias",
                x: &hs,
                y: &b1,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "kp22 bias",
                x: &hs,
                y: &b2,
                color: PALETTE[1],
                dashed: false,
            },
            LineSeries {
                label: "kp33 bias",
                x: &hs,
                y: &b3,
                color: PALETTE[2],
                dashed: false,
            },
        ],
    )?;

    // 3. Interval coverage & bias correction
    let cov_labels = [
        "kp11 Wald",
        "kp11 Pct",
        "kp11 Basic",
        "kp11 Stud",
        "kp11 LogB",
        "th1 Wald",
        "th1 Stud",
    ];
    let cov_vals = [
        pd.params[1].wald.coverage_95 * 100.0,
        pd.params[1].percentile.coverage_95 * 100.0,
        pd.params[1].basic.coverage_95 * 100.0,
        pd.params[1].studentized.coverage_95 * 100.0,
        pd.params[1].log_basic.coverage_95 * 100.0,
        pd.params[4].wald.coverage_95 * 100.0,
        pd.params[4].studentized.coverage_95 * 100.0,
    ];
    let _ = pc;
    svg::write_bar_chart(
        &fig_dir.join("day004f_bias-correction-and-interval-coverage.svg"),
        "Cakupan Empiris Interval 95% per Metode (H=20y, R=100, B=100)",
        "Cakupan (%)",
        &cov_labels,
        &cov_vals,
        PALETTE[0],
    )?;

    // 4. Electricity sample prices 7d
    let mut px = Vec::new();
    let mut py = Vec::new();
    for (d, day) in pe.sample_window_7d.iter().enumerate() {
        for h in 0..24 {
            px.push((d * 24 + h) as f64);
            py.push(day.prices[h]);
        }
    }
    svg::write_line_chart(
        &fig_dir.join("day004f_electricity-prices-regimes-jumps.svg"),
        "Lintasan Harga Listrik Sintetis 7 Hari (Rezim, Harga Negatif & Jump)",
        "Jam Pasca Burn-In (0..167)",
        "Harga Spot (EUR/MWh)",
        &[LineSeries {
            label: "Spot Price",
            x: &px,
            y: &py,
            color: PALETTE[0],
            dashed: false,
        }],
    )?;

    // 5. Battery dispatch SoC
    let ex = &pe.eval_days[pe.example_day_idx];
    let bat = pe.policy_33.battery;
    let tr_dp33 = trace_example_day(ex, bat, Some(&pe.policy_33), false, false);
    let tr_th = trace_example_day(ex, bat, None, true, false);
    let tr_pf = trace_example_day(ex, bat, None, false, true);
    let hrs25: Vec<f64> = (0..=24).map(|h| h as f64).collect();
    svg::write_line_chart(
        &fig_dir.join("day004f_battery-dispatch-and-strategies.svg"),
        "Jejak SoC Baterai pada Hari Contoh (DP-33 vs Threshold vs Perfect Foresight)",
        "Jam (0..24)",
        "SoC (MWh)",
        &[
            LineSeries {
                label: "DP-33",
                x: &hrs25,
                y: &tr_dp33.soc_mwh,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Threshold",
                x: &hrs25,
                y: &tr_th.soc_mwh,
                color: PALETTE[3],
                dashed: true,
            },
            LineSeries {
                label: "Perfect Foresight",
                x: &hrs25,
                y: &tr_pf.soc_mwh,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )?;

    // 6. Expanded-tail & grid sensitivity
    let grid_labels = [
        "DP-17",
        "DP-33",
        "Expanded-38",
        "Threshold",
        "Perfect Foresight",
    ];
    let grid_vals = [
        pe.res_dp17.mean_profit,
        pe.res_dp33.mean_profit,
        pe.res_exp38.mean_profit,
        pe.res_threshold.mean_profit,
        pe.res_pf.mean_profit,
    ];
    svg::write_bar_chart(
        &fig_dir.join("day004f_expanded-tail-grid-sensitivity.svg"),
        "Sensitivitas Profit Harian terhadap Resolusi Grid & Audit Expanded-Tail [-1000, 5000]",
        "Mean Profit Harian (EUR)",
        &grid_labels,
        &grid_vals,
        PALETTE[2],
    )?;

    // 7. Integration Panel F: Project KRD profile
    let krd_labels = ["1y", "2y", "3y", "5y", "7y", "10y", "20y", "30y"];
    let _ = PriceGrid::expanded_tail_38();
    svg::write_bar_chart(
        &fig_dir.join("day004f_battery-afns-npv-krd.svg"),
        &format!(
            "Panel F: Profil KRD Proyek Baterai 10y (NPV Pre-CAPEX = {:.0} EUR, Dur = {:.2}y)",
            pf.base_npv_pre_capex_eur, pf.eff_duration_closed
        ),
        "Key-Rate Duration (Tahun)",
        &krd_labels,
        &pf.krd.krd_closed,
        PALETTE[4],
    )?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_results_json(
    quick: bool,
    n_workers: usize,
    total_seconds: f64,
    pa: &panel_a::PanelAResult,
    pb: &panel_b::PanelBResult,
    pc: &panel_c::PanelCResult,
    pd: &panel_d::PanelDResult,
    pe: &panel_e::PanelEResult,
    pf: &integration::IntegrationResult,
    c_anchors: bool,
    c_pca: bool,
    c_jk_alg: bool,
    c_power: bool,
    c_battery_dp: bool,
    c_integration: bool,
    c_budget: bool,
    all_technical_passed: bool,
) -> JsonVal {
    let status = |ok: bool| if ok { "tercapai" } else { "belum tercapai" };

    JsonVal::Obj(vec![
        ("metadata", JsonVal::Obj(vec![
            ("session", "day004f".into()),
            ("date_wita", "2026-11-22".into()),
            ("language", "rust".into()),
            ("quick_mode", quick.into()),
            ("workers_used", n_workers.into()),
            ("seeds", JsonVal::Obj(vec![
                ("day004a", (rng::SEED_DAY004A as usize).into()),
                ("day004b", (rng::SEED_DAY004B as usize).into()),
                ("day004c", (rng::SEED_DAY004C as usize).into()),
                ("day004d", (rng::SEED_DAY004D as usize).into()),
                ("day004e_main", (rng::SEED_DAY004E_MAIN as usize).into()),
                ("day004e_kernel", (rng::SEED_DAY004E_KERNEL as usize).into()),
                ("day004e_eval", (rng::SEED_DAY004E_EVAL as usize).into()),
                ("day004e_verify", (rng::SEED_DAY004E_VERIFY as usize).into()),
                ("day004f_integration", (rng::SEED_DAY004F_INTEGRATION as usize).into()),
            ])),
            ("elapsed_seconds", JsonVal::Obj(vec![
                ("panel_a_day004a", pa.elapsed_seconds.into()),
                ("panel_b_day004b", pb.elapsed_seconds.into()),
                ("panel_c_day004c", pc.elapsed_seconds.into()),
                ("panel_d_day004d", pd.elapsed_seconds.into()),
                ("panel_e_day004e", pe.elapsed_seconds.into()),
                ("panel_f_integration", pf.elapsed_seconds.into()),
                ("total", total_seconds.into()),
            ])),
        ])),
        ("panel_a_day004a", JsonVal::Obj(vec![
            ("max_integral_err_ns", pa.max_integral_err_ns.into()),
            ("max_integral_err_nss", pa.max_integral_err_nss.into()),
            ("max_bootstrap_df_err", pa.max_bootstrap_df_err.into()),
            ("max_bootstrap_yield_err", pa.max_bootstrap_yield_err.into()),
            ("interp_loglin_neg_forward_mean", pa.noise_exp.log_linear_p.neg_forward_mean_count.into()),
            ("interp_loglin_forward_tv_mean", pa.noise_exp.log_linear_p.forward_tv_mean.into()),
            ("interp_loglin_yield_rmse_bps", (pa.noise_exp.log_linear_p.yield_rmse_mean * 10_000.0).into()),
            ("interp_spline_neg_forward_mean", pa.noise_exp.cubic_spline.neg_forward_mean_count.into()),
            ("interp_spline_forward_tv_mean", pa.noise_exp.cubic_spline.forward_tv_mean.into()),
            ("interp_spline_yield_rmse_bps", (pa.noise_exp.cubic_spline.yield_rmse_mean * 10_000.0).into()),
            ("bond_10y_4pct_mod_dur_closed", pa.bond_risk.mod_duration_ytm_closed.into()),
            ("bond_10y_4pct_mod_dur_rel_err", pa.bond_risk.mod_duration_ytm_rel_err.into()),
            ("bond_10y_4pct_conv_closed", pa.bond_risk.convexity_ytm_closed.into()),
            ("bond_10y_4pct_conv_rel_err", pa.bond_risk.convexity_ytm_rel_err.into()),
            ("krd_sum_vs_parallel_closed_err", pa.krd.sum_vs_par_closed_err.into()),
            ("swap_10y_par_rate", pa.swap.par_swap_rate.into()),
            ("swap_10y_pv01_per_1m", pa.swap.annuity_pv01_per_1m.into()),
            ("two_bond_hedge_res_5y", pa.hedge.residual_node1.into()),
            ("two_bond_hedge_res_10y", pa.hedge.residual_node2.into()),
            ("afns_riccati_max_a_err", pa.afns_riccati_max_a_err.into()),
            ("afns_adj_10y_bps", pa.afns_adj_10y_bps.into()),
            ("afns_adj_30y_bps", pa.afns_adj_30y_bps.into()),
            ("afns_phi_max_err", pa.afns_phi_max_err.into()),
            ("afns_q_van_loan_max_err", pa.afns_q_vl_max_err.into()),
            ("afns_mc_transition_max_z", pa.afns_mc_transition.max_z_score.into()),
            ("afns_tp_10y_closed_bps", pa.afns_tp_10y_closed_bps.into()),
            ("afns_tp_10y_mc_bps", pa.afns_tp_10y_mc_bps.into()),
            ("afns_tp_10y_z", pa.afns_tp_10y_z.into()),
            ("kalman_anchor_t5_ll_diff", pa.anchor_t5.log_lik_abs_diff.into()),
            ("kalman_anchor_t5_mean_diff", pa.anchor_t5.max_smooth_mean_abs_diff.into()),
            ("kalman_anchor_t5_cov_diff", pa.anchor_t5.max_smooth_cov_abs_diff.into()),
            ("single_mle_est_ll", pa.single_mle.est_log_lik.into()),
            ("single_mle_true_ll", pa.single_mle.true_log_lik.into()),
            ("dns_ll_drop", pa.dns_ll_drop.into()),
            ("dns_lambda_bias", pa.dns_lambda_bias.into()),
            ("mle_r50_ll_ge_true", pa.mle_ll_ge_true_count.into()),
            ("mle_r50_params", JsonVal::Arr(pa.mle_param_summaries.iter().map(|p| JsonVal::Obj(vec![
                ("name", p.name.into()),
                ("true", p.true_val.into()),
                ("bias", p.bias.into()),
                ("bias_mcse", p.bias_mcse.into()),
                ("sd", p.sd.into()),
                ("rmse", p.rmse.into()),
                ("coverage_95", p.coverage_95.into()),
            ])).collect())),
            ("pca_top3_cum_var_theo", pa.pca_summary.theo_top3_cum_var.into()),
            ("pca_top3_cos_mean", pa.pca_summary.emp_mean_top3_cos.to_vec().into()),
            ("pca_emp_sign_changes", pa.pca_summary.emp_sign_changes.to_vec().into()),
            ("pca_noise_bulk_z", pa.pca_summary.noise_bulk_z_score.into()),
            ("pca_aggregate_passed", pa.pca_acceptance.aggregate_passed.into()),
            ("pca_diag_all10_within_3se", pa.pca_acceptance.diagnostic_all10_sorted_within_3se.into()),
        ])),
        ("panel_b_day004b", JsonVal::Obj(vec![
            ("baseline_kp11_bias", pb.baseline.params[1].bias.into()),
            ("baseline_theta1_cov95", pb.baseline.params[4].coverage_95.into()),
            ("longspan_horizons", pb.longspan.iter().map(|s| JsonVal::Obj(vec![
                ("name", s.spec.name.clone().into()),
                ("horizon_years", s.spec.horizon_years.into()),
                ("kp11_bias", s.params[1].bias.into()),
                ("kp11_bias_mcse", s.params[1].bias_mcse.into()),
                ("kp22_bias", s.params[2].bias.into()),
                ("kp33_bias", s.params[3].bias.into()),
                ("theta1_cov95", s.params[4].coverage_95.into()),
                ("grad_tol_count", s.grad_tol_count.into()),
                ("line_search_failed_count", s.line_search_failed_count.into()),
                ("max_iter_count", s.max_iter_count.into()),
                ("bound_hit_count", s.bound_hit_count.into()),
            ])).collect::<Vec<_>>().into()),
            ("infill_20y", pb.infill.iter().map(|s| JsonVal::Obj(vec![
                ("name", s.spec.name.clone().into()),
                ("t_steps", s.spec.t_steps.into()),
                ("kp11_rmse", s.params[1].rmse.into()),
                ("theta1_rmse", s.params[4].rmse.into()),
                ("lambda_rmse", s.params[0].rmse.into()),
                ("sigma_eps_rmse", s.params[10].rmse.into()),
            ])).collect::<Vec<_>>().into()),
            ("kp_inv_h_trends", pb.kp_trends.iter().map(|t| JsonVal::Obj(vec![
                ("param", t.param_name.into()),
                ("ols_slope", t.ols_slope.into()),
                ("ols_slope_mcse", t.ols_slope_mcse.into()),
                ("wls_slope", t.wls_slope.into()),
                ("wls_slope_mcse", t.wls_slope_mcse.into()),
                ("bias_times_h", t.bias_times_h.to_vec().into()),
            ])).collect::<Vec<_>>().into()),
        ])),
        ("panel_c_day004c", JsonVal::Obj(vec![
            ("algebraic_c_over_h_err", pc.algebraic_c_over_h_max_err.into()),
            ("jackknife_20y_kp11", JsonVal::Obj(vec![
                ("raw_bias", pc.jackknife_20y.params[1].raw_bias.into()),
                ("jk_bias", pc.jackknife_20y.params[1].jk_bias.into()),
                ("raw_rmse", pc.jackknife_20y.params[1].raw_rmse.into()),
                ("jk_rmse", pc.jackknife_20y.params[1].jk_rmse.into()),
                ("non_physical_count", pc.jackknife_20y.params[1].non_physical_count.into()),
            ])),
            ("jackknife_80y_kp11", JsonVal::Obj(vec![
                ("raw_bias", pc.jackknife_80y.params[1].raw_bias.into()),
                ("jk_bias", pc.jackknife_80y.params[1].jk_bias.into()),
                ("raw_rmse", pc.jackknife_80y.params[1].raw_rmse.into()),
                ("jk_rmse", pc.jackknife_80y.params[1].jk_rmse.into()),
                ("non_physical_count", pc.jackknife_80y.params[1].non_physical_count.into()),
            ])),
            ("bootstrap_20y_total_refits", pc.bootstrap_20y.total_refits.into()),
            ("bootstrap_20y_finite_refits", pc.bootstrap_20y.finite_refits.into()),
            ("bootstrap_20y_kp11_pct_cov95", pc.bootstrap_20y.params[1].pct_cov_95.into()),
            ("bootstrap_20y_kp11_boot_corr_rmse", pc.bootstrap_20y.params[1].boot_corr_rmse.into()),
            ("bootstrap_20y_kp11_non_physical", pc.bootstrap_20y.params[1].non_physical_count.into()),
            ("bootstrap_20y_theta1_pct_cov95", pc.bootstrap_20y.params[4].pct_cov_95.into()),
        ])),
        ("panel_d_day004d", JsonVal::Obj(vec![
            ("total_inner_refits", pd.total_inner_refits.into()),
            ("inner_finite", pd.inner_opt_diag.finite_count.into()),
            ("kp11_intervals", JsonVal::Obj(vec![
                ("wald_cov", pd.params[1].wald.coverage_95.into()),
                ("pct_cov", pd.params[1].percentile.coverage_95.into()),
                ("basic_cov", pd.params[1].basic.coverage_95.into()),
                ("basic_non_phys", pd.params[1].basic.non_physical_count.into()),
                ("basic_phys_valid_cov", pd.params[1].basic.physically_valid_covered_count.into()),
                ("stud_cov", pd.params[1].studentized.coverage_95.into()),
                ("log_basic_cov", pd.params[1].log_basic.coverage_95.into()),
                ("log_basic_non_phys", pd.params[1].log_basic.non_physical_count.into()),
            ])),
            ("theta1_intervals", JsonVal::Obj(vec![
                ("wald_cov", pd.params[4].wald.coverage_95.into()),
                ("wald_med_width", pd.params[4].wald.median_width.into()),
                ("pct_cov", pd.params[4].percentile.coverage_95.into()),
                ("basic_cov", pd.params[4].basic.coverage_95.into()),
                ("stud_cov", pd.params[4].studentized.coverage_95.into()),
                ("stud_med_width", pd.params[4].studentized.median_width.into()),
                ("stud_iqr_width", pd.params[4].studentized.iqr_width.into()),
                ("stud_mean_width", pd.params[4].studentized.mean_width.into()),
            ])),
            ("stratification_kp11_median", pd.stratification.median_outer_kp11.into()),
            ("strat_low_kp11_theta1_stud_cov", pd.stratification.low_kp11_theta1_stud_cov.into()),
            ("strat_high_kp11_theta1_stud_cov", pd.stratification.high_kp11_theta1_stud_cov.into()),
        ])),
        ("panel_e_day004e", JsonVal::Obj(vec![
            ("stationary_pi", pe.stationary_pi.to_vec().into()),
            ("stationary_residual", pe.stationary_residual.into()),
            ("neg_price_pct", pe.path_diag.neg_price_pct.into()),
            ("total_nonzero_jumps", pe.path_diag.total_nonzero_jumps.into()),
            ("dp17_mean_profit", pe.res_dp17.mean_profit.into()),
            ("dp17_mcse", pe.res_dp17.mcse_profit.into()),
            ("dp17_upper_tail_mass", pe.kernel_17.upper_tail_mass.into()),
            ("dp17_lower_tail_mass", pe.kernel_17.lower_tail_mass.into()),
            ("dp33_mean_profit", pe.res_dp33.mean_profit.into()),
            ("dp33_mcse", pe.res_dp33.mcse_profit.into()),
            ("dp33_upper_tail_mass", pe.kernel_33.upper_tail_mass.into()),
            ("dp33_lower_tail_mass", pe.kernel_33.lower_tail_mass.into()),
            ("exp38_mean_profit", pe.res_exp38.mean_profit.into()),
            ("exp38_upper_tail_mass", pe.kernel_exp38.upper_tail_mass.into()),
            ("exp38_lower_tail_mass", pe.kernel_exp38.lower_tail_mass.into()),
            ("exp38_eval_agree_vs_dp33", pe.eval_agree_exp38_vs_33.into()),
            ("exp38_interior_agree_vs_dp33", pe.interior_agree_exp38_vs_33.into()),
            ("threshold_mean_profit", pe.res_threshold.mean_profit.into()),
            ("pf_mean_profit", pe.res_pf.mean_profit.into()),
            ("dp33_minus_thresh", pe.diff_33_vs_thresh.mean_diff.into()),
            ("dp33_minus_dp17", pe.diff_33_vs_17.mean_diff.into()),
            ("exp38_minus_dp33", pe.diff_exp38_vs_33.mean_diff.into()),
        ])),
        ("panel_f_integration", JsonVal::Obj(vec![
            ("disclaimer", "Studi valuasi sintetis/risk-free-discounted expected cashflow (bebas korelasi antara simulasi harga listrik dan estimasi kurva AFNS); bukan harga pasar arbitrage-free.".into()),
            ("project_years", pf.project_years.into()),
            ("mean_annual_cf_eur", pf.mean_annual_cf_eur.into()),
            ("sd_annual_cf_eur", pf.sd_annual_cf_eur.into()),
            ("base_npv_pre_capex_eur", pf.base_npv_pre_capex_eur.into()),
            ("break_even_capex_eur", pf.break_even_capex_eur.into()),
            ("resampled_npv_mean_eur", pf.resampled_npv_mean_eur.into()),
            ("resampled_npv_sd_eur", pf.resampled_npv_sd_eur.into()),
            ("resampled_npv_q05_eur", pf.resampled_npv_q05_eur.into()),
            ("resampled_npv_median_eur", pf.resampled_npv_median_eur.into()),
            ("resampled_npv_q95_eur", pf.resampled_npv_q95_eur.into()),
            ("eff_duration_closed", pf.eff_duration_closed.into()),
            ("eff_duration_fd", pf.eff_duration_fd.into()),
            ("duration_rel_err", pf.duration_rel_err.into()),
            ("eff_convexity_closed", pf.eff_convexity_closed.into()),
            ("convexity_rel_err", pf.convexity_rel_err.into()),
            ("parallel_pv01_eur", pf.parallel_pv01_eur.into()),
            ("krd_nodes", STANDARD_KRD_NODES.to_vec().into()),
            ("krd_closed", pf.krd.krd_closed.clone().into()),
            ("krd_sum_vs_par_err", pf.krd.sum_vs_par_closed_err.into()),
            ("hedge_weight_5y_4pct", pf.hedge.weight_bond_5y_4pct.into()),
            ("hedge_weight_10y_2pct", pf.hedge.weight_bond_10y_2pct.into()),
            ("hedge_residual_5y", pf.hedge.residual_node1.into()),
            ("hedge_residual_10y", pf.hedge.residual_node2.into()),
            ("afns_param_sensitivity", pf.afns_param_sensitivities.iter().map(|s| JsonVal::Obj(vec![
                ("source", s.source.into()),
                ("total_fits", s.total_fits.into()),
                ("finite_fits", s.finite_fits.into()),
                ("non_converged_count", s.non_converged_count.into()),
                ("bound_hit_count", s.bound_hit_count.into()),
                ("npv_mean_eur", s.npv_mean.into()),
                ("npv_sd_eur", s.npv_sd.into()),
                ("npv_q05_eur", s.npv_q05.into()),
                ("npv_median_eur", s.npv_median.into()),
                ("npv_q95_eur", s.npv_q95.into()),
            ])).collect::<Vec<_>>().into()),
            ("illustrative_capex_scenarios", pf.capex_scenarios.iter().map(|c| JsonVal::Obj(vec![
                ("capex_eur_assumed", c.capex_eur.into()),
                ("annual_om_eur_assumed", c.annual_om_eur.into()),
                ("post_capex_npv_eur", c.post_capex_npv_eur.into()),
            ])).collect::<Vec<_>>().into()),
        ])),
        ("acceptance_criteria", JsonVal::Obj(vec![
            ("c1_anchors_passed", c_anchors.into()),
            ("c1_status", status(c_anchors).into()),
            ("c2_pca_aggregate_passed", c_pca.into()),
            ("c2_status", status(c_pca).into()),
            ("c3_jackknife_algebra_passed", c_jk_alg.into()),
            ("c3_status", status(c_jk_alg).into()),
            ("c4_power_model_verified", c_power.into()),
            ("c4_status", status(c_power).into()),
            ("c5_battery_dp_verified", c_battery_dp.into()),
            ("c5_status", status(c_battery_dp).into()),
            ("c6_integration_anchors_passed", c_integration.into()),
            ("c6_status", status(c_integration).into()),
            ("c7_within_10min_budget", c_budget.into()),
            ("c7_status", status(c_budget).into()),
            ("all_technical_passed", all_technical_passed.into()),
        ])),
    ])
}
