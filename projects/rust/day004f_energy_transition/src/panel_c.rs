//! Panel C: complete port of day004c (temporal jackknife on H=20y & 80y weekly, R=50;
//! parametric bootstrap 20y weekly R=50 x B=100 = 5,000 refits).

use crate::afns::simulate_stationary_path;
use crate::interval::{self, MIN_VALID_FRACTION};
use crate::jackknife;
use crate::kalman::{
    estimate_mle, perturb_initial_guess, MleResult, ModelSpec, DEFAULT_SIGMA_EPS,
    IS_POSITIVE_PARAM, LOWER_BOUNDS, PARAM_NAMES, STANDARD_TENORS, UPPER_BOUNDS, WEEKLY_DT,
};
use crate::optim;
use crate::parallel;
use crate::rng::{derive_boot_seed, derive_seed, Rng, SEED_DAY004C};
use crate::stats;

#[derive(Clone, Debug)]
pub struct SegmentDiag {
    pub segment: String,
    pub total_reps: usize,
    pub finite_count: usize,
    pub converged_count: usize,
    pub scaled_grad_le_1e6_count: usize,
    pub ll_ge_true_count: usize,
    pub bound_hit_count: usize,
    pub hessian_reg_count: usize,
    pub grad_tol_count: usize,
    pub boundary_stationary_count: usize,
    pub line_search_failed_count: usize,
    pub max_iter_count: usize,
    pub invalid_obj_count: usize,
    pub ls_obj_bit_identical_count: usize,
    pub ls_armijo_not_met_count: usize,
    pub ls_step_underflow_count: usize,
}

fn summarize_segment(name: &str, reps: &[&MleResult]) -> SegmentDiag {
    let mut d = SegmentDiag {
        segment: name.into(),
        total_reps: reps.len(),
        finite_count: 0,
        converged_count: 0,
        scaled_grad_le_1e6_count: 0,
        ll_ge_true_count: 0,
        bound_hit_count: 0,
        hessian_reg_count: 0,
        grad_tol_count: 0,
        boundary_stationary_count: 0,
        line_search_failed_count: 0,
        max_iter_count: 0,
        invalid_obj_count: 0,
        ls_obj_bit_identical_count: 0,
        ls_armijo_not_met_count: 0,
        ls_step_underflow_count: 0,
    };
    for &r in reps {
        if r.finite {
            d.finite_count += 1;
        }
        if r.converged {
            d.converged_count += 1;
        }
        if r.scaled_grad_norm_opt <= 1e-6 {
            d.scaled_grad_le_1e6_count += 1;
        }
        if r.finite && r.est_log_lik >= r.true_log_lik - 1e-8 {
            d.ll_ge_true_count += 1;
        }
        if r.hit_bound {
            d.bound_hit_count += 1;
        }
        if r.hessian_regularized {
            d.hessian_reg_count += 1;
        }
        match r.termination_reason.as_str() {
            optim::TERM_GRADIENT_TOL => d.grad_tol_count += 1,
            optim::TERM_BOUNDARY_STATIONARY => d.boundary_stationary_count += 1,
            optim::TERM_LINE_SEARCH_FAILED => {
                d.line_search_failed_count += 1;
                if r.obj_bit_identical {
                    d.ls_obj_bit_identical_count += 1;
                }
                match r.line_search_failure_cause.as_str() {
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
pub struct JackknifeParamSummary {
    pub param_idx: usize,
    pub name: &'static str,
    pub true_value: f64,
    pub raw_mean: f64,
    pub raw_bias: f64,
    pub raw_bias_mcse: f64,
    pub raw_sd: f64,
    pub raw_rmse: f64,
    pub jk_mean: f64,
    pub jk_bias: f64,
    pub jk_bias_mcse: f64,
    pub jk_sd: f64,
    pub jk_rmse: f64,
    pub rmse_ratio: f64,
    pub oob_count: usize,
    pub non_physical_count: usize,
    pub paired_bias_diff: f64,
    pub paired_bias_diff_mcse: f64,
    pub wald_cov_95: f64,
    pub wald_wilson_lo: f64,
    pub wald_wilson_hi: f64,
}

#[derive(Clone, Debug)]
pub struct JackknifeHorizonResult {
    pub horizon_years: f64,
    pub t_steps: usize,
    pub num_reps: usize,
    pub elapsed_seconds: f64,
    pub full_diag: SegmentDiag,
    pub half1_diag: SegmentDiag,
    pub half2_diag: SegmentDiag,
    pub params: [JackknifeParamSummary; 11],
    pub full_fits: Vec<MleResult>,
}

#[derive(Clone, Debug)]
pub struct BootParamSummaryC {
    pub param_idx: usize,
    pub name: &'static str,
    pub true_value: f64,
    pub raw_bias: f64,
    pub raw_rmse: f64,
    pub boot_mean_cond_bias: f64,
    pub boot_corr_mean: f64,
    pub boot_corr_bias: f64,
    pub boot_corr_bias_mcse: f64,
    pub boot_corr_sd: f64,
    pub boot_corr_rmse: f64,
    pub rmse_ratio: f64,
    pub oob_count: usize,
    pub non_physical_count: usize,
    pub wald_cov_95: f64,
    pub wald_wilson_lo: f64,
    pub wald_wilson_hi: f64,
    pub pct_cov_95: f64,
    pub pct_wilson_lo: f64,
    pub pct_wilson_hi: f64,
    pub pct_valid_outer: usize,
    pub pct_mean_width: f64,
    pub pct_median_width: f64,
}

#[derive(Clone, Debug)]
pub struct ParametricBootStudyC {
    pub horizon_years: f64,
    pub t_steps: usize,
    pub num_outer_reps: usize,
    pub num_inner_boot: usize,
    pub total_refits: usize,
    pub elapsed_seconds: f64,
    pub finite_refits: usize,
    pub converged_refits: usize,
    pub scaled_grad_le_1e6_refits: usize,
    pub bound_hit_refits: usize,
    pub grad_tol_refits: usize,
    pub line_search_failed_refits: usize,
    pub max_iter_refits: usize,
    pub ls_obj_bit_identical: usize,
    pub ls_armijo_not_met: usize,
    pub ls_step_underflow: usize,
    pub valid_outer_intervals: usize,
    pub params: [BootParamSummaryC; 11],
}

#[derive(Clone, Debug)]
pub struct PanelCResult {
    pub elapsed_seconds: f64,
    pub algebraic_c_over_h_max_err: f64,
    pub jackknife_20y: JackknifeHorizonResult,
    pub jackknife_80y: JackknifeHorizonResult,
    pub bootstrap_20y: ParametricBootStudyC,
}

pub fn perturb_bootstrap_guess(outer_fit: ModelSpec, rg: &mut Rng) -> ModelSpec {
    let v = outer_fit.to_vec();
    let mut init_v = [0.0; 11];
    for i in 0..11 {
        let mut mult = rg.uniform_range(0.92, 1.08);
        if (mult - 1.0).abs() < 0.02 {
            mult = if mult < 1.0 { 0.95 } else { 1.05 };
        }
        let lo = LOWER_BOUNDS[i];
        let hi = UPPER_BOUNDS[i];
        let mut cand = if IS_POSITIVE_PARAM[i] {
            (v[i] * mult).clamp(lo * 1.05, hi * 0.95)
        } else if v[i].abs() < 0.005 {
            (v[i] + rg.uniform_range(-0.002, 0.002)).clamp(lo + 0.005, hi - 0.005)
        } else {
            (v[i] * mult).clamp(lo + 0.005, hi - 0.005)
        };
        if !cand.is_finite() {
            cand = 0.5 * (lo + hi);
        }
        init_v[i] = cand;
    }
    ModelSpec::from_vec(init_v, outer_fit.use_afns_adj)
}

fn run_jackknife_horizon(
    true_spec: ModelSpec,
    horizon_years: f64,
    t_steps: usize,
    shared_paths: &[(Vec<Vec<f64>>, ModelSpec)],
    n_workers: usize,
) -> JackknifeHorizonResult {
    let th = std::time::Instant::now();
    let n_reps = shared_paths.len();
    let half_t = t_steps / 2;

    let rep_fits: Vec<(
        MleResult,
        MleResult,
        MleResult,
        jackknife::JackknifeEvalResult,
    )> = parallel::run_chunks(n_reps, n_workers, |r| {
        let (ref y_full, init_sp) = shared_paths[r];
        let full = estimate_mle(
            init_sp,
            true_spec,
            &STANDARD_TENORS,
            WEEKLY_DT,
            &y_full[..t_steps],
            true,
        );
        let half1 = estimate_mle(
            init_sp,
            true_spec,
            &STANDARD_TENORS,
            WEEKLY_DT,
            &y_full[..half_t],
            false,
        );
        let half2 = estimate_mle(
            init_sp,
            true_spec,
            &STANDARD_TENORS,
            WEEKLY_DT,
            &y_full[half_t..t_steps],
            false,
        );
        let jk = jackknife::combine_vec(&full.est_params, &half1.est_params, &half2.est_params);
        (full, half1, half2, jk)
    });

    let full_refs: Vec<&MleResult> = rep_fits.iter().map(|x| &x.0).collect();
    let half1_refs: Vec<&MleResult> = rep_fits.iter().map(|x| &x.1).collect();
    let half2_refs: Vec<&MleResult> = rep_fits.iter().map(|x| &x.2).collect();
    let full_diag = summarize_segment("full", &full_refs);
    let half1_diag = summarize_segment("half1", &half1_refs);
    let half2_diag = summarize_segment("half2", &half2_refs);

    let true_v = true_spec.to_vec();
    let params: [JackknifeParamSummary; 11] = std::array::from_fn(|i| {
        let raw_vals: Vec<f64> = rep_fits.iter().map(|x| x.0.est_params[i]).collect();
        let jk_vals: Vec<f64> = rep_fits.iter().map(|x| x.3.est_params[i]).collect();
        let diffs: Vec<f64> = raw_vals
            .iter()
            .zip(jk_vals.iter())
            .map(|(r, j)| j - r)
            .collect();
        let oob = rep_fits.iter().filter(|x| x.3.out_of_bounds[i]).count();
        let non_phys = rep_fits.iter().filter(|x| x.3.non_physical[i]).count();
        let wald_cov_cnt = rep_fits
            .iter()
            .filter(|x| {
                x.0.hessian_se[i] > 0.0
                    && (x.0.est_params[i] - true_v[i]).abs() <= stats::Z95 * x.0.hessian_se[i]
            })
            .count();
        let (w_cov, w_lo, w_hi) = stats::wilson_ci(wald_cov_cnt, n_reps, stats::Z95);
        let raw_mean = stats::mean(&raw_vals);
        let jk_mean = stats::mean(&jk_vals);
        let raw_rmse = stats::rmse(&raw_vals, true_v[i]);
        let jk_rmse = stats::rmse(&jk_vals, true_v[i]);
        JackknifeParamSummary {
            param_idx: i,
            name: PARAM_NAMES[i],
            true_value: true_v[i],
            raw_mean,
            raw_bias: raw_mean - true_v[i],
            raw_bias_mcse: stats::se(&raw_vals),
            raw_sd: stats::sd(&raw_vals),
            raw_rmse,
            jk_mean,
            jk_bias: jk_mean - true_v[i],
            jk_bias_mcse: stats::se(&jk_vals),
            jk_sd: stats::sd(&jk_vals),
            jk_rmse,
            rmse_ratio: if raw_rmse > 0.0 {
                jk_rmse / raw_rmse
            } else {
                0.0
            },
            oob_count: oob,
            non_physical_count: non_phys,
            paired_bias_diff: stats::mean(&diffs),
            paired_bias_diff_mcse: stats::se(&diffs),
            wald_cov_95: w_cov,
            wald_wilson_lo: w_lo,
            wald_wilson_hi: w_hi,
        }
    });

    let full_fits = rep_fits.into_iter().map(|x| x.0).collect();
    JackknifeHorizonResult {
        horizon_years,
        t_steps,
        num_reps: n_reps,
        elapsed_seconds: th.elapsed().as_secs_f64(),
        full_diag,
        half1_diag,
        half2_diag,
        params,
        full_fits,
    }
}

pub fn run_panel_c(n_workers: usize, quick: bool) -> PanelCResult {
    let t0 = std::time::Instant::now();
    let true_spec = ModelSpec::default_spec();
    let n_reps = if quick { 3 } else { 50 };
    let n_boot = if quick { 6 } else { 100 };

    let t20 = if quick { 260 } else { 1040 };
    let t80 = if quick { 1040 } else { 4160 };

    // Algebraic check
    let alg_err = [20.0_f64, 80.0_f64]
        .iter()
        .map(|&h| {
            let full = 0.05 + 4.0 / h;
            let half = 0.05 + 4.0 / (0.5 * h);
            (jackknife::combine_split(full, half, half) - 0.05).abs()
        })
        .fold(0.0f64, f64::max);

    let shared_paths = parallel::run_chunks(n_reps, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(SEED_DAY004C, 0x1001, r as u64));
        let (_, y) = simulate_stationary_path(
            &true_spec.afns,
            &STANDARD_TENORS,
            t80,
            WEEKLY_DT,
            DEFAULT_SIGMA_EPS,
            &mut rg,
        );
        let init_sp = perturb_initial_guess(true_spec, &mut rg);
        (y, init_sp)
    });

    let jk20 = run_jackknife_horizon(true_spec, 20.0, t20, &shared_paths, n_workers);
    let jk80 = run_jackknife_horizon(true_spec, 80.0, t80, &shared_paths, n_workers);

    // Parametric bootstrap on 20y weekly
    let tb = std::time::Instant::now();
    let total_refits = n_reps * n_boot;
    let boot_draws = parallel::run_chunks(total_refits, n_workers, |job_idx| {
        let r = job_idx / n_boot;
        let b = job_idx % n_boot;
        let outer = &jk20.full_fits[r];
        let outer_spec = ModelSpec::from_vec(outer.est_params, true);
        let mut rg = Rng::new(derive_boot_seed(SEED_DAY004C, 0x5001, r as u64, b as u64));
        let (_, y_boot) = simulate_stationary_path(
            &outer_spec.afns,
            &STANDARD_TENORS,
            t20,
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
            false,
        )
    });

    let mut fin_cnt = 0;
    let mut conv_cnt = 0;
    let mut sgrad_cnt = 0;
    let mut hit_cnt = 0;
    let mut grad_tol_cnt = 0;
    let mut ls_fail_cnt = 0;
    let mut max_iter_cnt = 0;
    let mut ls_obj_bit = 0;
    let mut ls_armijo = 0;
    let mut ls_underflow = 0;

    for d in &boot_draws {
        if d.finite {
            fin_cnt += 1;
        }
        if d.converged {
            conv_cnt += 1;
        }
        if d.scaled_grad_norm_opt <= 1e-6 {
            sgrad_cnt += 1;
        }
        if d.hit_bound {
            hit_cnt += 1;
        }
        match d.termination_reason.as_str() {
            optim::TERM_GRADIENT_TOL => grad_tol_cnt += 1,
            optim::TERM_LINE_SEARCH_FAILED => {
                ls_fail_cnt += 1;
                if d.obj_bit_identical {
                    ls_obj_bit += 1;
                }
                match d.line_search_failure_cause.as_str() {
                    optim::LS_CAUSE_ARMIJO_NOT_MET => ls_armijo += 1,
                    optim::LS_CAUSE_STEP_UNDERFLOW => ls_underflow += 1,
                    _ => {}
                }
            }
            optim::TERM_MAX_ITERATIONS => max_iter_cnt += 1,
            _ => {}
        }
    }

    let min_req = interval::min_required_valid(n_boot);
    let mut valid_outer = 0usize;
    for r in 0..n_reps {
        let fin_r = (0..n_boot)
            .filter(|&b| boot_draws[r * n_boot + b].finite)
            .count();
        if fin_r >= min_req {
            valid_outer += 1;
        }
    }

    let true_v = true_spec.to_vec();
    let params: [BootParamSummaryC; 11] = std::array::from_fn(|i| {
        let mut raw_vals = Vec::with_capacity(n_reps);
        let mut corr_vals = Vec::with_capacity(n_reps);
        let mut cond_biases = Vec::with_capacity(n_reps);
        let mut oob = 0;
        let mut non_phys = 0;
        let mut wald_cov = 0;
        let mut pct_cov = 0;
        let mut pct_valid = 0;
        let mut pct_widths = Vec::new();

        for r in 0..n_reps {
            let outer = &jk20.full_fits[r];
            let psi_hat = outer.est_params[i];
            raw_vals.push(psi_hat);
            if outer.hessian_se[i] > 0.0
                && (psi_hat - true_v[i]).abs() <= stats::Z95 * outer.hessian_se[i]
            {
                wald_cov += 1;
            }
            let boot_ests: Vec<f64> = (0..n_boot)
                .map(|b| boot_draws[r * n_boot + b].est_params[i])
                .collect();
            let fin_ests: Vec<f64> = boot_ests
                .iter()
                .copied()
                .filter(|x| x.is_finite())
                .collect();
            let cb = if !fin_ests.is_empty() {
                stats::mean(&fin_ests) - psi_hat
            } else {
                0.0
            };
            cond_biases.push(cb);
            let corr = psi_hat - cb;
            corr_vals.push(corr);
            if corr < LOWER_BOUNDS[i] || corr > UPPER_BOUNDS[i] {
                oob += 1;
            }
            if IS_POSITIVE_PARAM[i] && corr <= 0.0 {
                non_phys += 1;
            }
            let pct_res = interval::eval_percentile(i, true_v[i], outer.finite, &boot_ests, n_boot);
            if pct_res.valid {
                pct_valid += 1;
                pct_widths.push(pct_res.width);
                if pct_res.covered {
                    pct_cov += 1;
                }
            }
        }

        let (w_cov, w_lo, w_hi) = stats::wilson_ci(wald_cov, n_reps, stats::Z95);
        let (p_cov, p_lo, p_hi) = stats::wilson_ci(pct_cov, n_reps, stats::Z95);
        let raw_rmse = stats::rmse(&raw_vals, true_v[i]);
        let corr_mean = stats::mean(&corr_vals);
        let corr_rmse = stats::rmse(&corr_vals, true_v[i]);

        BootParamSummaryC {
            param_idx: i,
            name: PARAM_NAMES[i],
            true_value: true_v[i],
            raw_bias: stats::mean(&raw_vals) - true_v[i],
            raw_rmse,
            boot_mean_cond_bias: stats::mean(&cond_biases),
            boot_corr_mean: corr_mean,
            boot_corr_bias: corr_mean - true_v[i],
            boot_corr_bias_mcse: stats::se(&corr_vals),
            boot_corr_sd: stats::sd(&corr_vals),
            boot_corr_rmse: corr_rmse,
            rmse_ratio: if raw_rmse > 0.0 {
                corr_rmse / raw_rmse
            } else {
                0.0
            },
            oob_count: oob,
            non_physical_count: non_phys,
            wald_cov_95: w_cov,
            wald_wilson_lo: w_lo,
            wald_wilson_hi: w_hi,
            pct_cov_95: p_cov,
            pct_wilson_lo: p_lo,
            pct_wilson_hi: p_hi,
            pct_valid_outer: pct_valid,
            pct_mean_width: stats::mean(&pct_widths),
            pct_median_width: stats::quantile(&pct_widths, 0.5),
        }
    });

    let bootstrap_20y = ParametricBootStudyC {
        horizon_years: 20.0,
        t_steps: t20,
        num_outer_reps: n_reps,
        num_inner_boot: n_boot,
        total_refits,
        elapsed_seconds: tb.elapsed().as_secs_f64(),
        finite_refits: fin_cnt,
        converged_refits: conv_cnt,
        scaled_grad_le_1e6_refits: sgrad_cnt,
        bound_hit_refits: hit_cnt,
        grad_tol_refits: grad_tol_cnt,
        line_search_failed_refits: ls_fail_cnt,
        max_iter_refits: max_iter_cnt,
        ls_obj_bit_identical: ls_obj_bit,
        ls_armijo_not_met: ls_armijo,
        ls_step_underflow: ls_underflow,
        valid_outer_intervals: valid_outer,
        params,
    };

    let _ = MIN_VALID_FRACTION;
    PanelCResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        algebraic_c_over_h_max_err: alg_err,
        jackknife_20y: jk20,
        jackknife_80y: jk80,
        bootstrap_20y,
    }
}
