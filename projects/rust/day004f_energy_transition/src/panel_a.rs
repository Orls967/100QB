//! Panel A: complete port of day004a (NS/NSS, bootstrap, interpolation, bond/KRD/swap/hedge,
//! AFNS/Riccati/P-dynamics/term premium, Kalman/RTS/MLE/R=50/DNS misspec, PCA T=5000 R=100).

use crate::afns::{
    simulate_stationary_path, simulate_stationary_path_opt, verify_long_run_tp_mc,
    verify_transition_mc, AfnsParams,
};
use crate::bootstrap;
use crate::interp::{self, NoiseExperimentSummary};
use crate::kalman::{
    estimate_mle, fast_log_likelihood, filter_and_smooth, perturb_initial_guess,
    verify_joint_gaussian_anchor, JointGaussianAnchor, MleResult, ModelSpec, DEFAULT_SIGMA_EPS,
    PARAM_NAMES, STANDARD_TENORS, WEEKLY_DT,
};
use crate::ns::{self, NsParams};
use crate::parallel;
use crate::pca::{self, PcaAcceptance, PcaReplicationSummary};
use crate::risk::{self, BondRiskMetrics, KrdReport, SwapPv01Report, TwoBondHedgeReport};
use crate::rng::{derive_seed, Rng, SEED_DAY004A};
use crate::stats;

#[derive(Clone, Debug)]
pub struct PanelAResult {
    pub elapsed_seconds: f64,
    // Part 1
    pub max_integral_err_ns: f64,
    pub max_integral_err_nss: f64,
    pub max_bootstrap_df_err: f64,
    pub max_bootstrap_yield_err: f64,
    pub small_tau_rows: Vec<ns::SmallTauRow>,
    pub noise_exp: NoiseExperimentSummary,
    pub bond_risk: BondRiskMetrics,
    pub krd: KrdReport,
    pub swap: SwapPv01Report,
    pub hedge: TwoBondHedgeReport,
    // Part 2
    pub afns_riccati_max_a_err: f64,
    pub afns_riccati_max_adj_err: f64,
    pub afns_adj_10y_bps: f64,
    pub afns_adj_30y_bps: f64,
    pub afns_phi_max_err: f64,
    pub afns_q_vl_max_err: f64,
    pub afns_q_gl_max_err: f64,
    pub afns_mc_transition: crate::afns::TransitionMcReport,
    pub afns_tp_10y_closed_bps: f64,
    pub afns_tp_10y_mc_bps: f64,
    pub afns_tp_10y_mc_se_bps: f64,
    pub afns_tp_10y_z: f64,
    // Part 3
    pub anchor_t5: JointGaussianAnchor,
    pub single_mle: MleResult,
    pub dns_misspec_mle: MleResult,
    pub dns_ll_drop: f64,
    pub dns_lambda_bias: f64,
    pub mle_reps: Vec<MleResult>,
    pub mle_param_summaries: [MleParamSummary; 11],
    pub mle_ll_ge_true_count: usize,
    pub mle_converged_count: usize,
    pub mle_bound_hit_count: usize,
    // Part 4
    pub pca_summary: PcaReplicationSummary,
    pub pca_acceptance: PcaAcceptance,
}

#[derive(Clone, Debug)]
pub struct MleParamSummary {
    pub name: &'static str,
    pub true_val: f64,
    pub mean_est: f64,
    pub bias: f64,
    pub bias_mcse: f64,
    pub sd: f64,
    pub rmse: f64,
    pub mean_hess_se: f64,
    pub median_hess_se: f64,
    pub coverage_95: f64,
    pub covered_count: usize,
    pub total_reps: usize,
    pub bound_hit_count: usize,
}

pub fn run_panel_a(n_workers: usize, quick: bool) -> PanelAResult {
    let t0 = std::time::Instant::now();
    let ns_curve = NsParams::default_ns();
    let nss_curve = NsParams::default_nss();

    // A1: NS/NSS integral, small tau, bootstrap, interp, risk
    let mut max_int_ns = 0.0f64;
    let mut max_int_nss = 0.0f64;
    for &tau in &[0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] {
        max_int_ns =
            max_int_ns.max((ns_curve.yield_at(tau) - ns_curve.forward_integral_avg(tau, 32)).abs());
        max_int_nss = max_int_nss
            .max((nss_curve.yield_at(tau) - nss_curve.forward_integral_avg(tau, 32)).abs());
    }
    let small_tau_rows =
        ns::compare_small_tau_stability(&ns_curve, &[1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-15]);

    let mats = bootstrap::semi_annual_grid(30.0);
    let (pars, true_df, true_y) = bootstrap::par_rates_from_curve(&ns_curve, &mats, 0.5);
    let (boot_df, boot_y) = bootstrap::bootstrap_from_par_rates(&mats, &pars, 0.5).unwrap();
    let max_boot_df_err = boot_df
        .iter()
        .zip(true_df.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    let max_boot_y_err = boot_y
        .iter()
        .zip(true_y.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);

    let interp_reps = if quick { 20 } else { 200 };
    let noise_exp =
        interp::run_noise_experiment(&ns_curve, interp_reps, 0.00005, n_workers, SEED_DAY004A);

    let bond_risk = risk::analyze_bond_risk(&ns_curve, 10.0, 0.04, 2.0, 1e-6, 1e-4);
    let (b_times, b_cfs) = risk::bond_cashflows(10.0, 0.04, 2.0);
    let yfn = |t: f64| ns_curve.yield_at(t);
    let krd = risk::compute_krd_generic(&yfn, &b_times, &b_cfs, &risk::STANDARD_KRD_NODES, 1e-6);
    let swap = risk::analyze_swap_pv01(&ns_curve, 10.0, 2.0);
    let hedge = risk::hedge_two_krds(&ns_curve, &risk::STANDARD_KRD_NODES, 3, 5, 1e-6).unwrap();

    // A2: AFNS, Riccati, transition, term premium
    let afns_p = AfnsParams::default_params();
    let mut max_a_err = 0.0f64;
    let mut max_adj_err = 0.0f64;
    for &tau in &[0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] {
        let (adj_c, a_c) = afns_p.yield_adjustment(tau);
        let (_, _, _, a_r) = afns_p.solve_riccati_rk4(tau, 400);
        max_a_err = max_a_err.max((a_c - a_r).abs());
        max_adj_err = max_adj_err.max((adj_c - (-a_r / tau)).abs());
    }
    let afns_adj_10y_bps = afns_p.yield_adjustment(10.0).0 * 10_000.0;
    let afns_adj_30y_bps = afns_p.yield_adjustment(30.0).0 * 10_000.0;

    let dt_1y = 1.0;
    let phi_closed = afns_p.transition_matrix_closed(dt_1y);
    let mut neg_kp_dt = crate::linalg::zeros(3, 3);
    for i in 0..3 {
        for j in 0..3 {
            neg_kp_dt[i][j] = -afns_p.kp[i][j] * dt_1y;
        }
    }
    let phi_series = crate::linalg::mat_exp_series(&neg_kp_dt, 35);
    let mut afns_phi_max_err = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            afns_phi_max_err = afns_phi_max_err.max((phi_closed[i][j] - phi_series[i][j]).abs());
        }
    }
    let q_c = afns_p.transition_cov_closed(dt_1y);
    let q_vl = afns_p.transition_cov_van_loan(dt_1y);
    let q_gl = afns_p.transition_cov_gl(dt_1y, 32);
    let mut afns_q_vl_max_err = 0.0f64;
    let mut afns_q_gl_max_err = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            afns_q_vl_max_err = afns_q_vl_max_err.max((q_c[i][j] - q_vl[i][j]).abs());
            afns_q_gl_max_err = afns_q_gl_max_err.max((q_c[i][j] - q_gl[i][j]).abs());
        }
    }
    let mc_n_trans = if quick { 10_000 } else { 200_000 };
    let afns_mc_transition = verify_transition_mc(
        &afns_p,
        [0.05, -0.02, 0.01],
        dt_1y,
        mc_n_trans,
        n_workers.max(4),
        n_workers,
        SEED_DAY004A,
    );
    let mc_n_tp = if quick { 10_000 } else { 100_000 };
    let (tp_c, tp_m, tp_se, tp_z) = verify_long_run_tp_mc(
        &afns_p,
        10.0,
        mc_n_tp,
        n_workers.max(4),
        n_workers,
        SEED_DAY004A,
    );

    // A3: Kalman anchor T=5, single MLE, DNS misspec, R=50 MLE
    let true_spec = ModelSpec::default_spec();
    let anchor_taus = [0.5, 2.0, 10.0];
    let mut rg_anc = Rng::new(derive_seed(SEED_DAY004A, 0x3101, 0));
    let (_, y_anc) = simulate_stationary_path(
        &true_spec.afns,
        &anchor_taus,
        5,
        WEEKLY_DT,
        DEFAULT_SIGMA_EPS,
        &mut rg_anc,
    );
    let anchor_t5 =
        verify_joint_gaussian_anchor(&true_spec, &anchor_taus, WEEKLY_DT, &y_anc).unwrap();

    let t_mle = if quick { 250 } else { 1000 };
    let mut rg_single = Rng::new(derive_seed(SEED_DAY004A, 0x3201, 0));
    let (_, y_single) = simulate_stationary_path(
        &true_spec.afns,
        &STANDARD_TENORS,
        t_mle,
        WEEKLY_DT,
        DEFAULT_SIGMA_EPS,
        &mut rg_single,
    );
    let init_single = perturb_initial_guess(true_spec, &mut rg_single);
    let single_mle = estimate_mle(
        init_single,
        true_spec,
        &STANDARD_TENORS,
        WEEKLY_DT,
        &y_single,
        true,
    );
    let _ = filter_and_smooth(&true_spec, &STANDARD_TENORS, WEEKLY_DT, &y_single);
    let _ = fast_log_likelihood(&true_spec, &STANDARD_TENORS, WEEKLY_DT, &y_single);

    // DNS misspecification: fit with use_afns_adj = false on true AFNS data
    let mut rg_dns = Rng::new(derive_seed(SEED_DAY004A, 0x3301, 0));
    let (_, y_dns) = simulate_stationary_path_opt(
        &true_spec.afns,
        &STANDARD_TENORS,
        t_mle,
        WEEKLY_DT,
        DEFAULT_SIGMA_EPS,
        true,
        &mut rg_dns,
    );
    let mut init_dns = perturb_initial_guess(true_spec, &mut rg_dns);
    init_dns.use_afns_adj = false;
    let mut true_dns = true_spec;
    true_dns.use_afns_adj = false;
    let dns_misspec_mle = estimate_mle(
        init_dns,
        true_dns,
        &STANDARD_TENORS,
        WEEKLY_DT,
        &y_dns,
        true,
    );
    let dns_ll_drop = single_mle.est_log_lik - dns_misspec_mle.est_log_lik;
    let dns_lambda_bias = dns_misspec_mle.est_params[0] - true_spec.afns.lambda;

    let n_mle_reps = if quick { 4 } else { 50 };
    let mle_reps = parallel::run_chunks(n_mle_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004A, 0x3401, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t_mle,
            WEEKLY_DT,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        estimate_mle(init_sp, true_spec, &STANDARD_TENORS, WEEKLY_DT, &y, true)
    });

    let true_v = true_spec.to_vec();
    let mle_ll_ge_true_count = mle_reps
        .iter()
        .filter(|r| r.finite && r.est_log_lik >= r.true_log_lik - 1e-8)
        .count();
    let mle_converged_count = mle_reps.iter().filter(|r| r.converged).count();
    let mle_bound_hit_count = mle_reps.iter().filter(|r| r.hit_bound).count();

    let mle_param_summaries: [MleParamSummary; 11] = std::array::from_fn(|i| {
        let ests: Vec<f64> = mle_reps.iter().map(|r| r.est_params[i]).collect();
        let ses: Vec<f64> = mle_reps.iter().map(|r| r.hessian_se[i]).collect();
        let cov_cnt = mle_reps
            .iter()
            .filter(|r| {
                r.hessian_se[i] > 0.0
                    && (r.est_params[i] - true_v[i]).abs() <= stats::Z95 * r.hessian_se[i]
            })
            .count();
        let hit_cnt = mle_reps.iter().filter(|r| r.hit_bound_by_param[i]).count();
        let mean_est = stats::mean(&ests);
        MleParamSummary {
            name: PARAM_NAMES[i],
            true_val: true_v[i],
            mean_est,
            bias: mean_est - true_v[i],
            bias_mcse: stats::se(&ests),
            sd: stats::sd(&ests),
            rmse: stats::rmse(&ests, true_v[i]),
            mean_hess_se: stats::mean(&ses),
            median_hess_se: stats::quantile(&ses, 0.5),
            coverage_95: cov_cnt as f64 / n_mle_reps.max(1) as f64,
            covered_count: cov_cnt,
            total_reps: n_mle_reps,
            bound_hit_count: hit_cnt,
        }
    });

    // A4: PCA
    let pca_t = if quick { 500 } else { 5000 };
    let pca_r = if quick { 10 } else { 100 };
    let pca_summary = pca::run_pca_study(
        &true_spec,
        &STANDARD_TENORS,
        WEEKLY_DT,
        pca_t,
        pca_r,
        n_workers,
        SEED_DAY004A,
    );
    let pca_acceptance = pca::evaluate_acceptance(&pca_summary);

    PanelAResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        max_integral_err_ns: max_int_ns,
        max_integral_err_nss: max_int_nss,
        max_bootstrap_df_err: max_boot_df_err,
        max_bootstrap_yield_err: max_boot_y_err,
        small_tau_rows,
        noise_exp,
        bond_risk,
        krd,
        swap,
        hedge,
        afns_riccati_max_a_err: max_a_err,
        afns_riccati_max_adj_err: max_adj_err,
        afns_adj_10y_bps,
        afns_adj_30y_bps,
        afns_phi_max_err,
        afns_q_vl_max_err,
        afns_q_gl_max_err,
        afns_mc_transition,
        afns_tp_10y_closed_bps: tp_c * 10_000.0,
        afns_tp_10y_mc_bps: tp_m * 10_000.0,
        afns_tp_10y_mc_se_bps: tp_se * 10_000.0,
        afns_tp_10y_z: tp_z,
        anchor_t5,
        single_mle,
        dns_misspec_mle,
        dns_ll_drop,
        dns_lambda_bias,
        mle_reps,
        mle_param_summaries,
        mle_ll_ge_true_count,
        mle_converged_count,
        mle_bound_hit_count,
        pca_summary,
        pca_acceptance,
    }
}
