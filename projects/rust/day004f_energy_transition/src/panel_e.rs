//! Panel E: complete port of day004e + expanded-tail [-1000, 5000] sensitivity audit.

use crate::battery::BatterySpec;
use crate::dp::{
    build_kernel, compute_regret_vs_pf, eval_no_battery, eval_perfect_foresight,
    eval_stochastic_dp_helper, eval_strategy, policy_action_agreement, solve_deterministic_dp,
    solve_stochastic_dp, threshold_action, DpPolicy, PriceGrid, StrategyResult, TransitionKernel,
};
use crate::power::{
    jump_theoretical_moments, sample_stationary_regime, simulate_eval_days,
    stationary_distribution, stationary_residual, DayPath, PowerModelSpec, HOURS_PER_DAY,
    NUM_REGIMES, REGIME_NAMES, REGIME_OVERSUPPLY, REGIME_SCARCITY,
};
use crate::rng::{derive_seed, Rng, SEED_DAY004E_EVAL, SEED_DAY004E_KERNEL, SEED_DAY004E_VERIFY};
use crate::stats::{self, PairedDiffSummary};

#[derive(Clone, Debug)]
pub struct PathDiagnostics {
    pub total_hours: usize,
    pub neg_price_hours: usize,
    pub neg_price_pct: f64,
    pub normal_hours: usize,
    pub scarcity_hours: usize,
    pub oversupply_hours: usize,
    pub total_nonzero_jumps: usize,
    pub min_price: f64,
    pub q01_price: f64,
    pub q05_price: f64,
    pub median_price: f64,
    pub mean_price: f64,
    pub q95_price: f64,
    pub q99_price: f64,
    pub max_price: f64,
}

#[derive(Clone, Debug)]
pub struct PanelEResult {
    pub elapsed_seconds: f64,
    pub stationary_pi: [f64; 3],
    pub stationary_residual: f64,
    pub jump_verify_ok: bool,
    pub step_verify_ok: bool,
    pub bf_abs_diff: f64,
    pub path_diag: PathDiagnostics,
    pub kernel_17: TransitionKernel,
    pub policy_17: DpPolicy,
    pub kernel_33: TransitionKernel,
    pub policy_33: DpPolicy,
    pub kernel_exp38: TransitionKernel,
    pub policy_exp38: DpPolicy,
    pub res_dp17: StrategyResult,
    pub res_dp33: StrategyResult,
    pub res_exp38: StrategyResult,
    pub res_threshold: StrategyResult,
    pub res_pf: StrategyResult,
    pub res_no_battery: StrategyResult,
    pub diff_33_vs_thresh: PairedDiffSummary,
    pub diff_17_vs_thresh: PairedDiffSummary,
    pub diff_33_vs_17: PairedDiffSummary,
    pub diff_exp38_vs_33: PairedDiffSummary,
    pub eval_agree_exp38_vs_33: f64,
    pub interior_agree_exp38_vs_33: f64,
    pub eval_days: Vec<DayPath>,
    pub example_day_idx: usize,
    pub sample_window_7d: Vec<DayPath>,
}

pub fn run_panel_e(n_workers: usize, quick: bool) -> PanelEResult {
    let t0 = std::time::Instant::now();
    let spec = PowerModelSpec::default_spec();
    let bat = BatterySpec::default_spec();
    let pi = stationary_distribution(&spec.transition_matrix).unwrap();
    let pi_res = stationary_residual(&spec.transition_matrix, pi);

    let mc_draws = if quick { 20_000 } else { 200_000 };
    let mut jump_ok = true;
    for r in 0..NUM_REGIMES {
        let mut rg = Rng::new(derive_seed(SEED_DAY004E_VERIFY, 0x2100 + r as u64, 0));
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..mc_draws {
            let n = rg.poisson(spec.jumps[r].lambda);
            let jv = if spec.jumps[r].sign == 0.0 || n == 0 {
                0.0
            } else {
                let mut s = 0.0;
                for _ in 0..n {
                    s += rg.exp(spec.jumps[r].mean_jump);
                }
                spec.jumps[r].sign * s
            };
            sum += jv;
            sum2 += jv * jv;
        }
        let emp_m = sum / mc_draws as f64;
        let (tm, tv) = jump_theoretical_moments(spec.jumps[r]);
        if r == 0 {
            if emp_m != 0.0 {
                jump_ok = false;
            }
        } else {
            let mcse_m = (tv / mc_draws as f64).sqrt();
            if (emp_m - tm).abs() > 4.0 * mcse_m {
                jump_ok = false;
            }
        }
        let _ = (sum2, REGIME_NAMES[r]);
    }

    let step_cases = [
        (0, 18, 55.0),
        (REGIME_SCARCITY, 23, 210.0),
        (REGIME_OVERSUPPLY, 4, -85.0),
    ];
    let mut step_ok = true;
    for (idx, &(r, h, p)) in step_cases.iter().enumerate() {
        let mut rg = Rng::new(derive_seed(SEED_DAY004E_VERIFY, 0x3100 + idx as u64, 0));
        let mut sum = 0.0;
        for _ in 0..mc_draws {
            sum += spec.step(r, h, p, &mut rg).1;
        }
        let emp_m = sum / mc_draws as f64;
        let (tm, tv) = spec.one_step_conditional_moments(r, h, p);
        let mcse_m = (tv / mc_draws as f64).sqrt();
        if (emp_m - tm).abs() > 4.0 * mcse_m {
            step_ok = false;
        }
    }

    let bf_prices = [-40.0, 35.0, 240.0, 10.0];
    let (_, dp_bf_gross, _, _) = solve_deterministic_dp(&bf_prices, bat);
    let bf_abs_diff = (dp_bf_gross - 626.20).abs();

    let n_draws = if quick { 500 } else { 10_000 };
    let n_days = if quick { 50 } else { 2_000 };

    let k17 = build_kernel(
        &spec,
        PriceGrid::coarse_17(),
        n_draws,
        n_workers,
        SEED_DAY004E_KERNEL,
        0x4E17,
        false,
    );
    let pol17 = solve_stochastic_dp(&k17, bat);

    let k33 = build_kernel(
        &spec,
        PriceGrid::fine_33(),
        n_draws,
        n_workers,
        SEED_DAY004E_KERNEL,
        0x4E33,
        false,
    );
    let pol33 = solve_stochastic_dp(&k33, bat);

    let k_exp = build_kernel(
        &spec,
        PriceGrid::expanded_tail_38(),
        n_draws,
        n_workers,
        SEED_DAY004E_KERNEL,
        0x4E38,
        true,
    );
    let pol_exp = solve_stochastic_dp(&k_exp, bat);

    let eval_days = simulate_eval_days(n_days, n_workers, &spec, pi, SEED_DAY004E_EVAL, 0x4E99);
    let path_diag = analyze_paths(&eval_days);

    let mut res17 = eval_stochastic_dp_helper("stochastic_dp_17bin", &eval_days, &pol17, bat);
    let mut res33 = eval_stochastic_dp_helper("stochastic_dp_33bin", &eval_days, &pol33, bat);
    let mut res_exp =
        eval_stochastic_dp_helper("stochastic_dp_expanded_38bin", &eval_days, &pol_exp, bat);
    let mut res_thresh = eval_strategy("threshold_policy", &eval_days, bat, |_, d, h, s| {
        threshold_action(d.prices[h], s, bat)
    });
    let res_pf = eval_perfect_foresight(&eval_days, bat);
    let mut res_nobat = eval_no_battery(&eval_days);

    compute_regret_vs_pf(&mut res17, &res_pf, 1e-9);
    compute_regret_vs_pf(&mut res33, &res_pf, 1e-9);
    compute_regret_vs_pf(&mut res_exp, &res_pf, 1e-9);
    compute_regret_vs_pf(&mut res_thresh, &res_pf, 1e-9);
    compute_regret_vs_pf(&mut res_nobat, &res_pf, 1e-9);

    let diff_33_vs_thresh = stats::paired_diff_ci(&res33.daily_profits, &res_thresh.daily_profits);
    let diff_17_vs_thresh = stats::paired_diff_ci(&res17.daily_profits, &res_thresh.daily_profits);
    let diff_33_vs_17 = stats::paired_diff_ci(&res33.daily_profits, &res17.daily_profits);
    let diff_exp38_vs_33 = stats::paired_diff_ci(&res_exp.daily_profits, &res33.daily_profits);
    let (eval_agree, _, _, int_agree, _, _) =
        policy_action_agreement(&eval_days, &pol33, &pol_exp, bat);

    let ex_idx = select_example_day(&eval_days);
    let sample_win = build_continuous_sample_window(&spec, pi);

    PanelEResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        stationary_pi: pi,
        stationary_residual: pi_res,
        jump_verify_ok: jump_ok,
        step_verify_ok: step_ok,
        bf_abs_diff,
        path_diag,
        kernel_17: k17,
        policy_17: pol17,
        kernel_33: k33,
        policy_33: pol33,
        kernel_exp38: k_exp,
        policy_exp38: pol_exp,
        res_dp17: res17,
        res_dp33: res33,
        res_exp38: res_exp,
        res_threshold: res_thresh,
        res_pf,
        res_no_battery: res_nobat,
        diff_33_vs_thresh,
        diff_17_vs_thresh,
        diff_33_vs_17,
        diff_exp38_vs_33,
        eval_agree_exp38_vs_33: eval_agree,
        interior_agree_exp38_vs_33: int_agree,
        eval_days,
        example_day_idx: ex_idx,
        sample_window_7d: sample_win,
    }
}

fn analyze_paths(days: &[DayPath]) -> PathDiagnostics {
    let mut all_p = Vec::with_capacity(days.len() * HOURS_PER_DAY);
    let mut neg_h = 0;
    let mut norm_h = 0;
    let mut scar_h = 0;
    let mut over_h = 0;
    let mut tot_j = 0;
    let mut min_p = f64::INFINITY;
    let mut max_p = f64::NEG_INFINITY;
    for d in days {
        for h in 0..HOURS_PER_DAY {
            let p = d.prices[h];
            all_p.push(p);
            min_p = min_p.min(p);
            max_p = max_p.max(p);
            if p < 0.0 {
                neg_h += 1;
            }
            match d.regimes[h] {
                0 => norm_h += 1,
                1 => {
                    scar_h += 1;
                    tot_j += d.num_jumps[h];
                }
                _ => {
                    over_h += 1;
                    tot_j += d.num_jumps[h];
                }
            }
        }
    }
    let tot_h = all_p.len();
    PathDiagnostics {
        total_hours: tot_h,
        neg_price_hours: neg_h,
        neg_price_pct: 100.0 * neg_h as f64 / tot_h.max(1) as f64,
        normal_hours: norm_h,
        scarcity_hours: scar_h,
        oversupply_hours: over_h,
        total_nonzero_jumps: tot_j,
        min_price: min_p,
        q01_price: stats::quantile(&all_p, 0.01),
        q05_price: stats::quantile(&all_p, 0.05),
        median_price: stats::quantile(&all_p, 0.50),
        mean_price: stats::mean(&all_p),
        q95_price: stats::quantile(&all_p, 0.95),
        q99_price: stats::quantile(&all_p, 0.99),
        max_price: max_p,
    }
}

fn select_example_day(days: &[DayPath]) -> usize {
    for (i, d) in days.iter().enumerate() {
        let has_neg = d.prices.iter().any(|&p| p < 0.0);
        let has_hi = d.prices.iter().any(|&p| p >= 200.0);
        if has_neg && has_hi {
            return i;
        }
    }
    0
}

fn build_continuous_sample_window(spec: &PowerModelSpec, pi: [f64; 3]) -> Vec<DayPath> {
    let mut rg = Rng::new(derive_seed(SEED_DAY004E_VERIFY, 0x5A00, 0));
    let mut r = sample_stationary_regime(pi, &mut rg);
    let mut p = spec.seasonal_mean(r, 0);
    for step in 0..(spec.burn_in_days * HOURS_PER_DAY) {
        (r, p, _, _) = spec.step(r, step, p, &mut rg);
    }
    let mut win = Vec::with_capacity(7);
    for d in 0..7 {
        let mut day = DayPath {
            day_idx: d,
            prices: [0.0; HOURS_PER_DAY],
            regimes: [0; HOURS_PER_DAY],
            num_jumps: [0; HOURS_PER_DAY],
            jump_vals: [0.0; HOURS_PER_DAY],
        };
        for h in 0..HOURS_PER_DAY {
            let (rn, pn, nj, jv) = spec.step(r, d * HOURS_PER_DAY + h, p, &mut rg);
            r = rn;
            p = pn;
            day.prices[h] = pn;
            day.regimes[h] = rn;
            day.num_jumps[h] = nj;
            day.jump_vals[h] = jv;
        }
        win.push(day);
    }
    win
}
