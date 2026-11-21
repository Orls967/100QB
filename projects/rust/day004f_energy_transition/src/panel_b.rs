//! Panel B: complete port of day004b (baseline 1000w, long-span 10/20/40/80y weekly with shared-prefix paths,
//! infill 20y monthly/weekly/daily, OLS/WLS 1/H bias trends, and optimizer diagnostics).

use crate::afns::simulate_stationary_path;
use crate::kalman::{
    estimate_mle, perturb_initial_guess, MleResult, ModelSpec, DEFAULT_SIGMA_EPS, PARAM_NAMES,
    STANDARD_TENORS, WEEKLY_DT,
};
use crate::optim;
use crate::parallel;
use crate::rng::{derive_seed, Rng, SEED_DAY004B};
use crate::stats;

#[derive(Clone, Debug)]
pub struct ConfigSpec {
    pub name: String,
    pub group: String,
    pub frequency_label: String,
    pub horizon_years: f64,
    pub dt: f64,
    pub t_steps: usize,
    pub num_reps: usize,
    pub salt: u64,
    pub reused_from: String,
}

#[derive(Clone, Debug)]
pub struct ParamSummaryB {
    pub name: &'static str,
    pub true_value: f64,
    pub mean_estimate: f64,
    pub bias: f64,
    pub bias_mcse: f64,
    pub bias_in_mcse: f64,
    pub empirical_sd: f64,
    pub rmse: f64,
    pub mean_hessian_se: f64,
    pub median_hessian_se: f64,
    pub mean_hess_to_emp_sd: f64,
    pub median_hess_to_emp_sd: f64,
    pub coverage_95: f64,
    pub coverage_wilson_lower: f64,
    pub coverage_wilson_upper: f64,
    pub covered_count: usize,
    pub total_reps: usize,
    pub bound_hit_count: usize,
    pub bound_hit_pct: f64,
}

#[derive(Clone, Debug)]
pub struct StudyConfigResult {
    pub spec: ConfigSpec,
    pub elapsed_seconds: f64,
    pub converged_count: usize,
    pub grad_tol_count: usize,
    pub boundary_stationary_count: usize,
    pub line_search_failed_count: usize,
    pub max_iter_count: usize,
    pub invalid_obj_count: usize,
    pub best_stage1_count: usize,
    pub best_stage2_count: usize,
    pub scaled_grad_le_1e6_count: usize,
    pub ll_ge_true_count: usize,
    pub bound_hit_count: usize,
    pub hessian_reg_count: usize,
    pub params: [ParamSummaryB; 11],
    pub replications: Vec<MleResult>,
}

#[derive(Clone, Debug)]
pub struct BiasInvHTrend {
    pub param_name: &'static str,
    pub true_value: f64,
    pub horizons: [f64; 4],
    pub biases: [f64; 4],
    pub bias_times_h: [f64; 4],
    pub ols_slope: f64,
    pub ols_intercept: f64,
    pub ols_slope_mcse: f64,
    pub wls_slope: f64,
    pub wls_intercept: f64,
    pub wls_slope_mcse: f64,
}

#[derive(Clone, Debug)]
pub struct PanelBResult {
    pub elapsed_seconds: f64,
    pub baseline: StudyConfigResult,
    pub longspan: [StudyConfigResult; 4],
    pub infill: [StudyConfigResult; 3],
    pub kp_trends: [BiasInvHTrend; 3],
}

fn summarize_replications(
    true_spec: ModelSpec,
    cfg: ConfigSpec,
    elapsed: f64,
    reps: Vec<MleResult>,
) -> StudyConfigResult {
    let true_v = true_spec.to_vec();
    let n_rep = reps.len();
    let mut conv_cnt = 0;
    let mut grad_tol = 0;
    let mut bnd_stat = 0;
    let mut ls_fail = 0;
    let mut max_iter = 0;
    let mut inv_obj = 0;
    let mut st1 = 0;
    let mut st2 = 0;
    let mut sgrad = 0;
    let mut ll_ge = 0;
    let mut hit = 0;
    let mut reg = 0;

    for r in &reps {
        if r.converged {
            conv_cnt += 1;
        }
        match r.termination_reason.as_str() {
            optim::TERM_GRADIENT_TOL => grad_tol += 1,
            optim::TERM_BOUNDARY_STATIONARY => bnd_stat += 1,
            optim::TERM_LINE_SEARCH_FAILED => ls_fail += 1,
            optim::TERM_MAX_ITERATIONS => max_iter += 1,
            optim::TERM_INVALID_OBJECTIVE => inv_obj += 1,
            _ => {}
        }
        if r.best_stage == 1 {
            st1 += 1;
        } else {
            st2 += 1;
        }
        if r.scaled_grad_norm_opt <= 1e-6 {
            sgrad += 1;
        }
        if r.finite && r.est_log_lik >= r.true_log_lik - 1e-8 {
            ll_ge += 1;
        }
        if r.hit_bound {
            hit += 1;
        }
        if r.hessian_regularized {
            reg += 1;
        }
    }

    let params: [ParamSummaryB; 11] = std::array::from_fn(|i| {
        let ests: Vec<f64> = reps.iter().map(|r| r.est_params[i]).collect();
        let ses: Vec<f64> = reps.iter().map(|r| r.hessian_se[i]).collect();
        let cov_cnt = reps
            .iter()
            .filter(|r| {
                r.hessian_se[i] > 0.0
                    && (r.est_params[i] - true_v[i]).abs() <= stats::Z95 * r.hessian_se[i]
            })
            .count();
        let hit_cnt = reps.iter().filter(|r| r.hit_bound_by_param[i]).count();
        let mean_est = stats::mean(&ests);
        let bias = mean_est - true_v[i];
        let emp_sd = stats::sd(&ests);
        let mcse = stats::se(&ests);
        let mean_hse = stats::mean(&ses);
        let med_hse = stats::quantile(&ses, 0.5);
        let (cov, w_lo, w_hi) = stats::wilson_ci(cov_cnt, n_rep, stats::Z95);
        ParamSummaryB {
            name: PARAM_NAMES[i],
            true_value: true_v[i],
            mean_estimate: mean_est,
            bias,
            bias_mcse: mcse,
            bias_in_mcse: if mcse > 0.0 { bias.abs() / mcse } else { 0.0 },
            empirical_sd: emp_sd,
            rmse: stats::rmse(&ests, true_v[i]),
            mean_hessian_se: mean_hse,
            median_hessian_se: med_hse,
            mean_hess_to_emp_sd: if emp_sd > 0.0 { mean_hse / emp_sd } else { 0.0 },
            median_hess_to_emp_sd: if emp_sd > 0.0 { med_hse / emp_sd } else { 0.0 },
            coverage_95: cov,
            coverage_wilson_lower: w_lo,
            coverage_wilson_upper: w_hi,
            covered_count: cov_cnt,
            total_reps: n_rep,
            bound_hit_count: hit_cnt,
            bound_hit_pct: 100.0 * hit_cnt as f64 / n_rep.max(1) as f64,
        }
    });

    StudyConfigResult {
        spec: cfg,
        elapsed_seconds: elapsed,
        converged_count: conv_cnt,
        grad_tol_count: grad_tol,
        boundary_stationary_count: bnd_stat,
        line_search_failed_count: ls_fail,
        max_iter_count: max_iter,
        invalid_obj_count: inv_obj,
        best_stage1_count: st1,
        best_stage2_count: st2,
        scaled_grad_le_1e6_count: sgrad,
        ll_ge_true_count: ll_ge,
        bound_hit_count: hit,
        hessian_reg_count: reg,
        params,
        replications: reps,
    }
}

pub fn run_panel_b(n_workers: usize, quick: bool) -> PanelBResult {
    let t0 = std::time::Instant::now();
    let true_spec = ModelSpec::default_spec();
    let n_reps = if quick { 4 } else { 50 };

    // 1. Baseline 1000w
    let tb = std::time::Instant::now();
    let t_base = if quick { 250 } else { 1000 };
    let base_reps = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004B, 0x1001, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t_base,
            WEEKLY_DT,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        estimate_mle(init_sp, true_spec, &STANDARD_TENORS, WEEKLY_DT, &y, true)
    });
    let baseline = summarize_replications(
        true_spec,
        ConfigSpec {
            name: "baseline_1000w".into(),
            group: "baseline".into(),
            frequency_label: "weekly".into(),
            horizon_years: t_base as f64 / 52.0,
            dt: WEEKLY_DT,
            t_steps: t_base,
            num_reps: n_reps,
            salt: 0x1001,
            reused_from: String::new(),
        },
        tb.elapsed().as_secs_f64(),
        base_reps,
    );

    // 2. Long-span weekly {10, 20, 40, 80}y with shared-prefix paths
    let horizons = [10.0, 20.0, 40.0, 80.0];
    let t_steps_ls = if quick {
        [130, 260, 520, 1040]
    } else {
        [520, 1040, 2080, 4160]
    };
    let max_ls_t = t_steps_ls[3];
    let shared_paths = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004B, 0x2001, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            max_ls_t,
            WEEKLY_DT,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        (y, init_sp)
    });

    let ls_results_vec: Vec<StudyConfigResult> = (0..4)
        .map(|idx| {
            let th = std::time::Instant::now();
            let t_len = t_steps_ls[idx];
            let reps = parallel::run_chunks(n_reps, n_workers, |r| {
                let (ref y_full, init_sp) = shared_paths[r];
                estimate_mle(
                    init_sp,
                    true_spec,
                    &STANDARD_TENORS,
                    WEEKLY_DT,
                    &y_full[..t_len],
                    true,
                )
            });
            summarize_replications(
                true_spec,
                ConfigSpec {
                    name: format!("longspan_{}y_weekly", horizons[idx] as usize),
                    group: "longspan".into(),
                    frequency_label: "weekly".into(),
                    horizon_years: horizons[idx],
                    dt: WEEKLY_DT,
                    t_steps: t_len,
                    num_reps: n_reps,
                    salt: 0x2001,
                    reused_from: String::new(),
                },
                th.elapsed().as_secs_f64(),
                reps,
            )
        })
        .collect();
    let longspan: [StudyConfigResult; 4] = ls_results_vec.try_into().unwrap();

    // 3. Infill 20y monthly (dt=1/12, T=240), weekly (reused from longspan_20y), daily (dt=1/252, T=5040)
    let tm = std::time::Instant::now();
    let t_mon = if quick { 60 } else { 240 };
    let mon_reps = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004B, 0x3001, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t_mon,
            1.0 / 12.0,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        estimate_mle(init_sp, true_spec, &STANDARD_TENORS, 1.0 / 12.0, &y, true)
    });
    let infill_monthly = summarize_replications(
        true_spec,
        ConfigSpec {
            name: "infill_20y_monthly".into(),
            group: "infill".into(),
            frequency_label: "monthly".into(),
            horizon_years: 20.0,
            dt: 1.0 / 12.0,
            t_steps: t_mon,
            num_reps: n_reps,
            salt: 0x3001,
            reused_from: String::new(),
        },
        tm.elapsed().as_secs_f64(),
        mon_reps,
    );

    let mut infill_weekly = longspan[1].clone();
    infill_weekly.spec.name = "infill_20y_weekly".into();
    infill_weekly.spec.group = "infill".into();
    infill_weekly.spec.reused_from = "longspan_20y_weekly".into();

    let td = std::time::Instant::now();
    let t_day = if quick { 504 } else { 5040 };
    let day_reps = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004B, 0x3003, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t_day,
            1.0 / 252.0,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        estimate_mle(init_sp, true_spec, &STANDARD_TENORS, 1.0 / 252.0, &y, true)
    });
    let infill_daily = summarize_replications(
        true_spec,
        ConfigSpec {
            name: "infill_20y_daily".into(),
            group: "infill".into(),
            frequency_label: "daily".into(),
            horizon_years: 20.0,
            dt: 1.0 / 252.0,
            t_steps: t_day,
            num_reps: n_reps,
            salt: 0x3003,
            reused_from: String::new(),
        },
        td.elapsed().as_secs_f64(),
        day_reps,
    );

    // 4. KPii bias vs 1/H OLS/WLS trends with prefix MCSE
    let inv_h: Vec<f64> = horizons.iter().map(|&h| 1.0 / h).collect();
    let kp_trends: [BiasInvHTrend; 3] = std::array::from_fn(|k| {
        let p_idx = k + 1;
        let biases = [
            longspan[0].params[p_idx].bias,
            longspan[1].params[p_idx].bias,
            longspan[2].params[p_idx].bias,
            longspan[3].params[p_idx].bias,
        ];
        let bias_times_h = [
            biases[0] * horizons[0],
            biases[1] * horizons[1],
            biases[2] * horizons[2],
            biases[3] * horizons[3],
        ];
        let weights: Vec<f64> = (0..4)
            .map(|h_i| {
                let s = longspan[h_i].params[p_idx].bias_mcse;
                if s > 0.0 {
                    1.0 / (s * s)
                } else {
                    1.0
                }
            })
            .collect();
        let (ols_s, ols_i) = stats::ols_slope_intercept(&inv_h, &biases);
        let (wls_s, wls_i) = stats::wls_slope_intercept(&inv_h, &biases, &weights);

        let mut rep_ols = Vec::with_capacity(n_reps);
        let mut rep_wls = Vec::with_capacity(n_reps);
        let tv = longspan[0].params[p_idx].true_value;
        for r in 0..n_reps {
            let rb: Vec<f64> = (0..4)
                .map(|h_i| longspan[h_i].replications[r].est_params[p_idx] - tv)
                .collect();
            rep_ols.push(stats::ols_slope_intercept(&inv_h, &rb).0);
            rep_wls.push(stats::wls_slope_intercept(&inv_h, &rb, &weights).0);
        }

        BiasInvHTrend {
            param_name: PARAM_NAMES[p_idx],
            true_value: tv,
            horizons,
            biases,
            bias_times_h,
            ols_slope: ols_s,
            ols_intercept: ols_i,
            ols_slope_mcse: stats::se(&rep_ols),
            wls_slope: wls_s,
            wls_intercept: wls_i,
            wls_slope_mcse: stats::se(&rep_wls),
        }
    });

    PanelBResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        baseline,
        longspan,
        infill: [infill_monthly, infill_weekly, infill_daily],
        kp_trends,
    }
}
