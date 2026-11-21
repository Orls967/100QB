//! Panel D: complete port of day004d (H=20y weekly T=1040, R=100 outer x B=100 inner = 10,000 refits
//! with Hessian SE*, comparing 5 interval methods + robust width summaries + kp11 stratification).

use crate::afns::simulate_stationary_path;
use crate::interval::{self, IntervalResult};
use crate::kalman::{
    estimate_mle, perturb_initial_guess, MleResult, ModelSpec, DEFAULT_SIGMA_EPS,
    IS_POSITIVE_PARAM, PARAM_NAMES, STANDARD_TENORS, WEEKLY_DT,
};
use crate::optim;
use crate::panel_c::perturb_bootstrap_guess;
use crate::parallel;
use crate::rng::{derive_boot_seed, derive_seed, Rng, SEED_DAY004D};
use crate::stats;

#[derive(Clone, Debug)]
pub struct OptAggregateDiag {
    pub total_fits: usize,
    pub finite_count: usize,
    pub converged_count: usize,
    pub scaled_grad_le_1e6_count: usize,
    pub bound_hit_count: usize,
    pub hessian_regularized_count: usize,
    pub grad_tol_count: usize,
    pub boundary_stationary_count: usize,
    pub line_search_failed_count: usize,
    pub max_iter_count: usize,
    pub invalid_obj_count: usize,
    pub ls_obj_bit_identical_count: usize,
    pub ls_armijo_not_met_count: usize,
    pub ls_step_underflow_count: usize,
}

fn aggregate_opt_diag(fits: &[&MleResult]) -> OptAggregateDiag {
    let mut d = OptAggregateDiag {
        total_fits: fits.len(),
        finite_count: 0,
        converged_count: 0,
        scaled_grad_le_1e6_count: 0,
        bound_hit_count: 0,
        hessian_regularized_count: 0,
        grad_tol_count: 0,
        boundary_stationary_count: 0,
        line_search_failed_count: 0,
        max_iter_count: 0,
        invalid_obj_count: 0,
        ls_obj_bit_identical_count: 0,
        ls_armijo_not_met_count: 0,
        ls_step_underflow_count: 0,
    };
    for &f in fits {
        if f.finite {
            d.finite_count += 1;
        }
        if f.converged {
            d.converged_count += 1;
        }
        if f.scaled_grad_norm_opt <= 1e-6 {
            d.scaled_grad_le_1e6_count += 1;
        }
        if f.hit_bound {
            d.bound_hit_count += 1;
        }
        if f.hessian_regularized {
            d.hessian_regularized_count += 1;
        }
        match f.termination_reason.as_str() {
            optim::TERM_GRADIENT_TOL => d.grad_tol_count += 1,
            optim::TERM_BOUNDARY_STATIONARY => d.boundary_stationary_count += 1,
            optim::TERM_LINE_SEARCH_FAILED => {
                d.line_search_failed_count += 1;
                if f.obj_bit_identical {
                    d.ls_obj_bit_identical_count += 1;
                }
                match f.line_search_failure_cause.as_str() {
                    optim::LS_CAUSE_ARMIJO_NOT_MET => d.ls_armijo_not_met_count += 1,
                    optim::LS_CAUSE_STEP_UNDERFLOW => d.ls_step_underflow_count += 1,
                    _ => {}
                }
            }
            optim::TERM_MAX_ITERATIONS => d.max_iter_count += 1,
            optim::TERM_INVALID_OBJECTIVE => d.invalid_obj_count += 1,
            _ => {}
        }
    }
    d
}

#[derive(Clone, Debug)]
pub struct MethodSummary {
    pub method: &'static str,
    pub applicable: bool,
    pub total_outer_reps: usize,
    pub valid_count: usize,
    pub invalid_count: usize,
    pub covered_count: usize,
    pub physically_valid_covered_count: usize,
    pub coverage_95: f64,
    pub wilson_lower: f64,
    pub wilson_upper: f64,
    pub mean_width: f64,
    pub median_width: f64,
    pub q25_width: f64,
    pub q75_width: f64,
    pub iqr_width: f64,
    pub p95_width: f64,
    pub max_width: f64,
    pub non_physical_count: usize,
    pub out_of_bounds_count: usize,
}

fn summarize_method(
    method: &'static str,
    applicable: bool,
    results: &[&IntervalResult],
) -> MethodSummary {
    let n = results.len();
    let mut valid_cnt = 0;
    let mut invalid_cnt = 0;
    let mut cov_cnt = 0;
    let mut phys_cov_cnt = 0;
    let mut non_phys_cnt = 0;
    let mut oob_cnt = 0;
    let mut widths = Vec::with_capacity(n);

    for &r in results {
        if r.valid {
            valid_cnt += 1;
            widths.push(r.width);
            if r.covered {
                cov_cnt += 1;
                if !r.non_physical {
                    phys_cov_cnt += 1;
                }
            }
            if r.non_physical {
                non_phys_cnt += 1;
            }
            if r.out_of_bounds {
                oob_cnt += 1;
            }
        } else {
            invalid_cnt += 1;
        }
    }
    let (cov, w_lo, w_hi) = stats::wilson_ci(cov_cnt, n, stats::Z95);
    let q25 = stats::quantile(&widths, 0.25);
    let q75 = stats::quantile(&widths, 0.75);
    let max_w = widths.iter().copied().fold(0.0f64, f64::max);

    MethodSummary {
        method,
        applicable,
        total_outer_reps: n,
        valid_count: valid_cnt,
        invalid_count: invalid_cnt,
        covered_count: cov_cnt,
        physically_valid_covered_count: phys_cov_cnt,
        coverage_95: cov,
        wilson_lower: w_lo,
        wilson_upper: w_hi,
        mean_width: stats::mean(&widths),
        median_width: stats::quantile(&widths, 0.5),
        q25_width: q25,
        q75_width: q75,
        iqr_width: q75 - q25,
        p95_width: stats::quantile(&widths, 0.95),
        max_width: if widths.is_empty() { f64::NAN } else { max_w },
        non_physical_count: non_phys_cnt,
        out_of_bounds_count: oob_cnt,
    }
}

#[derive(Clone, Debug)]
pub struct ParamCalibrationSummary {
    pub param_idx: usize,
    pub name: &'static str,
    pub true_value: f64,
    pub raw_bias: f64,
    pub raw_rmse: f64,
    pub wald: MethodSummary,
    pub percentile: MethodSummary,
    pub basic: MethodSummary,
    pub studentized: MethodSummary,
    pub log_basic: MethodSummary,
}

#[derive(Clone, Debug)]
pub struct StratificationSummary {
    pub median_outer_kp11: f64,
    pub low_kp11_count: usize,
    pub high_kp11_count: usize,
    pub low_kp11_theta1_wald_cov: f64,
    pub high_kp11_theta1_wald_cov: f64,
    pub low_kp11_theta1_stud_cov: f64,
    pub high_kp11_theta1_stud_cov: f64,
    pub low_kp11_theta1_wald_med_w: f64,
    pub high_kp11_theta1_wald_med_w: f64,
    pub low_kp11_theta1_stud_med_w: f64,
    pub high_kp11_theta1_stud_med_w: f64,
    pub corr_kp11_with_wald_w: f64,
    pub corr_kp11_with_stud_w: f64,
}

#[derive(Clone, Debug)]
pub struct PanelDResult {
    pub elapsed_seconds: f64,
    pub horizon_years: f64,
    pub t_steps: usize,
    pub num_outer_reps: usize,
    pub num_inner_boot: usize,
    pub total_inner_refits: usize,
    pub outer_opt_diag: OptAggregateDiag,
    pub inner_opt_diag: OptAggregateDiag,
    pub params: [ParamCalibrationSummary; 11],
    pub stratification: StratificationSummary,
    pub outer_fits: Vec<MleResult>,
}

pub fn run_panel_d(n_workers: usize, quick: bool) -> PanelDResult {
    let t0 = std::time::Instant::now();
    let true_spec = ModelSpec::default_spec();
    let t_steps = if quick { 260 } else { 1040 };
    let n_outer = if quick { 4 } else { 100 };
    let n_boot = if quick { 6 } else { 100 };

    let outer_fits = parallel::run_chunks(n_outer, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004D, 0x1001, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t_steps,
            WEEKLY_DT,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        estimate_mle(init_sp, true_spec, &STANDARD_TENORS, WEEKLY_DT, &y, true)
    });

    let total_inner = n_outer * n_boot;
    let inner_fits = parallel::run_chunks(total_inner, n_workers, |job_idx| {
        let r = job_idx / n_boot;
        let b = job_idx % n_boot;
        let outer = &outer_fits[r];
        let outer_spec = ModelSpec::from_vec(outer.est_params, true);
        let mut rg = Rng::new(derive_boot_seed(SEED_DAY004D, 0x5001, r as u64, b as u64));
        let (_, y_boot) = simulate_stationary_path(
            &outer_spec.afns,
            &STANDARD_TENORS,
            t_steps,
            WEEKLY_DT,
            outer_spec.sigma_eps,
            &mut rg,
        );
        let init_boot = perturb_bootstrap_guess(outer_spec, &mut rg);
        estimate_mle(
            init_boot,
            outer_spec,
            &STANDARD_TENORS,
            WEEKLY_DT,
            &y_boot,
            true,
        )
    });

    let outer_refs: Vec<&MleResult> = outer_fits.iter().collect();
    let inner_refs: Vec<&MleResult> = inner_fits.iter().collect();
    let outer_opt_diag = aggregate_opt_diag(&outer_refs);
    let inner_opt_diag = aggregate_opt_diag(&inner_refs);

    let true_v = true_spec.to_vec();
    let mut wald_all: Vec<[IntervalResult; 11]> = Vec::with_capacity(n_outer);
    let mut pct_all: Vec<[IntervalResult; 11]> = Vec::with_capacity(n_outer);
    let mut bas_all: Vec<[IntervalResult; 11]> = Vec::with_capacity(n_outer);
    let mut stud_all: Vec<[IntervalResult; 11]> = Vec::with_capacity(n_outer);
    let mut log_all: Vec<[IntervalResult; 11]> = Vec::with_capacity(n_outer);

    for r in 0..n_outer {
        let outer = &outer_fits[r];
        let w_arr = std::array::from_fn(|i| {
            interval::eval_wald(
                i,
                outer.est_params[i],
                outer.hessian_se[i],
                true_v[i],
                outer.finite,
            )
        });
        let p_arr = std::array::from_fn(|i| {
            let ests: Vec<f64> = (0..n_boot)
                .map(|b| inner_fits[r * n_boot + b].est_params[i])
                .collect();
            interval::eval_percentile(i, true_v[i], outer.finite, &ests, n_boot)
        });
        let b_arr = std::array::from_fn(|i| {
            let ests: Vec<f64> = (0..n_boot)
                .map(|b| inner_fits[r * n_boot + b].est_params[i])
                .collect();
            interval::eval_basic(
                i,
                outer.est_params[i],
                true_v[i],
                outer.finite,
                &ests,
                n_boot,
            )
        });
        let s_arr = std::array::from_fn(|i| {
            let ests: Vec<f64> = (0..n_boot)
                .map(|b| inner_fits[r * n_boot + b].est_params[i])
                .collect();
            let ses: Vec<f64> = (0..n_boot)
                .map(|b| inner_fits[r * n_boot + b].hessian_se[i])
                .collect();
            interval::eval_studentized(
                i,
                outer.est_params[i],
                outer.hessian_se[i],
                true_v[i],
                outer.finite,
                &ests,
                &ses,
                n_boot,
            )
        });
        let l_arr = std::array::from_fn(|i| {
            if !IS_POSITIVE_PARAM[i] {
                interval::eval_wald(
                    i,
                    outer.est_params[i],
                    outer.hessian_se[i],
                    true_v[i],
                    false,
                )
            } else {
                let ests: Vec<f64> = (0..n_boot)
                    .map(|b| inner_fits[r * n_boot + b].est_params[i])
                    .collect();
                interval::eval_log_basic(
                    i,
                    outer.est_params[i],
                    true_v[i],
                    outer.finite,
                    &ests,
                    n_boot,
                )
            }
        });
        wald_all.push(w_arr);
        pct_all.push(p_arr);
        bas_all.push(b_arr);
        stud_all.push(s_arr);
        log_all.push(l_arr);
    }

    let params: [ParamCalibrationSummary; 11] = std::array::from_fn(|i| {
        let raw_vals: Vec<f64> = outer_fits.iter().map(|r| r.est_params[i]).collect();
        let w_refs: Vec<&IntervalResult> = wald_all.iter().map(|a| &a[i]).collect();
        let p_refs: Vec<&IntervalResult> = pct_all.iter().map(|a| &a[i]).collect();
        let b_refs: Vec<&IntervalResult> = bas_all.iter().map(|a| &a[i]).collect();
        let s_refs: Vec<&IntervalResult> = stud_all.iter().map(|a| &a[i]).collect();
        let l_refs: Vec<&IntervalResult> = log_all.iter().map(|a| &a[i]).collect();
        ParamCalibrationSummary {
            param_idx: i,
            name: PARAM_NAMES[i],
            true_value: true_v[i],
            raw_bias: stats::mean(&raw_vals) - true_v[i],
            raw_rmse: stats::rmse(&raw_vals, true_v[i]),
            wald: summarize_method(interval::METHOD_WALD, true, &w_refs),
            percentile: summarize_method(interval::METHOD_PERCENTILE, true, &p_refs),
            basic: summarize_method(interval::METHOD_BASIC, true, &b_refs),
            studentized: summarize_method(interval::METHOD_STUDENTIZED, true, &s_refs),
            log_basic: summarize_method(interval::METHOD_LOG_BASIC, IS_POSITIVE_PARAM[i], &l_refs),
        }
    });

    // Stratification by outer kp11_hat
    let kp11_hats: Vec<f64> = outer_fits.iter().map(|r| r.est_params[1]).collect();
    let med_kp11 = stats::quantile(&kp11_hats, 0.5);
    let mut low_w_cov = 0usize;
    let mut high_w_cov = 0usize;
    let mut low_s_cov = 0usize;
    let mut high_s_cov = 0usize;
    let mut low_cnt = 0usize;
    let mut high_cnt = 0usize;
    let mut low_ww = Vec::new();
    let mut high_ww = Vec::new();
    let mut low_sw = Vec::new();
    let mut high_sw = Vec::new();
    let mut all_ww = Vec::with_capacity(n_outer);
    let mut all_sw = Vec::with_capacity(n_outer);

    for r in 0..n_outer {
        let w_res = &wald_all[r][4];
        let s_res = &stud_all[r][4];
        all_ww.push(w_res.width);
        all_sw.push(s_res.width);
        if kp11_hats[r] <= med_kp11 {
            low_cnt += 1;
            if w_res.covered {
                low_w_cov += 1;
            }
            if s_res.covered {
                low_s_cov += 1;
            }
            low_ww.push(w_res.width);
            low_sw.push(s_res.width);
        } else {
            high_cnt += 1;
            if w_res.covered {
                high_w_cov += 1;
            }
            if s_res.covered {
                high_s_cov += 1;
            }
            high_ww.push(w_res.width);
            high_sw.push(s_res.width);
        }
    }

    let stratification = StratificationSummary {
        median_outer_kp11: med_kp11,
        low_kp11_count: low_cnt,
        high_kp11_count: high_cnt,
        low_kp11_theta1_wald_cov: low_w_cov as f64 / low_cnt.max(1) as f64,
        high_kp11_theta1_wald_cov: high_w_cov as f64 / high_cnt.max(1) as f64,
        low_kp11_theta1_stud_cov: low_s_cov as f64 / low_cnt.max(1) as f64,
        high_kp11_theta1_stud_cov: high_s_cov as f64 / high_cnt.max(1) as f64,
        low_kp11_theta1_wald_med_w: stats::quantile(&low_ww, 0.5),
        high_kp11_theta1_wald_med_w: stats::quantile(&high_ww, 0.5),
        low_kp11_theta1_stud_med_w: stats::quantile(&low_sw, 0.5),
        high_kp11_theta1_stud_med_w: stats::quantile(&high_sw, 0.5),
        corr_kp11_with_wald_w: stats::correlation(&kp11_hats, &all_ww),
        corr_kp11_with_stud_w: stats::correlation(&kp11_hats, &all_sw),
    };

    PanelDResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        horizon_years: 20.0,
        t_steps,
        num_outer_reps: n_outer,
        num_inner_boot: n_boot,
        total_inner_refits: total_inner,
        outer_opt_diag,
        inner_opt_diag,
        params,
        stratification,
        outer_fits,
    }
}
