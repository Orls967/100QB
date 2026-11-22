#![allow(clippy::too_many_arguments)]

use crate::cpmm::{lvr_rate_analytic, pool_delta, pool_value};
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::stats;

/// Closed-form jump LVR rate per unit pool value for Merton log-normal jumps
/// `Y ~ N(mu_j, sigma_j^2)` with Poisson intensity `lambda_j`:
/// loss per jump = `(V/2)(e^{Y/2} - 1)^2`, so
/// `LVR_jump / V = (lambda_j / 2) * (exp(mu_j + sigma_j^2 / 2) - 2 * exp(mu_j / 2 + sigma_j^2 / 8) + 1)`.
#[inline]
pub fn merton_jump_lvr_rate(lambda_j: f64, mu_j: f64, sigma_j: f64) -> f64 {
    let e_exp_y = (mu_j + 0.5 * sigma_j * sigma_j).exp();
    let e_exp_half_y = (0.5 * mu_j + 0.125 * sigma_j * sigma_j).exp();
    0.5 * lambda_j * (e_exp_y - 2.0 * e_exp_half_y + 1.0)
}

/// Merton jump-diffusion (`lambda_j = 20/yr, mu_j = -0.02, sigma_j = 0.05`,
/// diffusion `sigma = 0.8`). Returns
/// `(lvr_realized_ann, lvr_diffusion_only_ann, lvr_total_pred_ann, se_realized)`.
pub fn misspec_merton_jumps(
    sigma: f64,
    lambda_j: f64,
    mu_j: f64,
    sigma_j: f64,
    t_yr: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let dt = t_yr / n_steps as f64;
    let kappa_j = (mu_j + 0.5 * sigma_j * sigma_j).exp() - 1.0;
    let drift = (-0.5 * sigma * sigma - lambda_j * kappa_j) * dt;
    let vol_dt = sigma * dt.sqrt();
    let l = 0.5_f64;
    let jump_rate = merton_jump_lvr_rate(lambda_j, mu_j, sigma_j);

    let rows = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x6001, r as u64));
        let mut p = 1.0_f64;
        let mut v_rebal = 1.0_f64;
        let mut lvr_diff_int = 0.0_f64;
        let mut lvr_jump_int = 0.0_f64;
        for _ in 0..n_steps {
            let v_t = pool_value(l, p);
            let d_t = pool_delta(l, p);
            lvr_diff_int += lvr_rate_analytic(v_t, sigma) * dt;
            lvr_jump_int += jump_rate * v_t * dt;

            let nj = rg.poisson(lambda_j * dt);
            let mut jump = 0.0_f64;
            for _ in 0..nj {
                jump += mu_j + sigma_j * rg.normal();
            }
            let p_next = p * (drift + vol_dt * rg.normal() + jump).exp();
            v_rebal += d_t * (p_next - p);
            p = p_next;
        }
        let lvr_real = (v_rebal - pool_value(l, p)) / t_yr;
        let lvr_diff = lvr_diff_int / t_yr;
        let lvr_pred = (lvr_diff_int + lvr_jump_int) / t_yr;
        (lvr_real, lvr_diff, lvr_pred)
    });

    let reals: Vec<f64> = rows.iter().map(|x| x.0).collect();
    let diffs: Vec<f64> = rows.iter().map(|x| x.1).collect();
    let preds: Vec<f64> = rows.iter().map(|x| x.2).collect();
    (
        stats::mean(&reals),
        stats::mean(&diffs),
        stats::mean(&preds),
        stats::se(&reals),
    )
}

/// Stochastic volatility (CIR/Heston with `E[v] = theta = sigma^2 = 0.64, kappa = 3, xi = 1`).
/// Returns `(lvr_realized_ann, expected_v_v_over_8_ann, se_realized, corr_lvr_v)`.
pub fn misspec_stoch_vol(
    theta: f64,
    kappa: f64,
    xi: f64,
    t_yr: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let dt = t_yr / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let l = 0.5_f64;

    let rows = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x6002, r as u64));
        let mut p = 1.0_f64;
        let mut v_var = theta;
        let mut v_rebal = 1.0_f64;
        let mut lvr_pred = 0.0_f64;
        let mut v_avg = 0.0_f64;
        for _ in 0..n_steps {
            let v_pos = v_var.max(0.0);
            let v_pool = pool_value(l, p);
            let d_t = pool_delta(l, p);
            lvr_pred += (v_pos * v_pool / 8.0) * dt;
            v_avg += v_pos * dt;

            let z_p = rg.normal();
            let z_v = rg.normal();
            let p_next = p * (-0.5 * v_pos * dt + v_pos.sqrt() * sqrt_dt * z_p).exp();
            v_var =
                (v_pos + kappa * (theta - v_pos) * dt + xi * v_pos.sqrt() * sqrt_dt * z_v).max(0.0);
            v_rebal += d_t * (p_next - p);
            p = p_next;
        }
        (
            (v_rebal - pool_value(l, p)) / t_yr,
            lvr_pred / t_yr,
            v_avg / t_yr,
        )
    });

    let reals: Vec<f64> = rows.iter().map(|x| x.0).collect();
    let preds: Vec<f64> = rows.iter().map(|x| x.1).collect();
    let vs: Vec<f64> = rows.iter().map(|x| x.2).collect();
    (
        stats::mean(&reals),
        stats::mean(&preds),
        stats::se(&reals),
        stats::corr(&reals, &vs),
    )
}

/// Non-zero drift `mu in {-0.5, 0.5}`: decomposes `V_pool - V_hodl` (IL) and `V_rebal - V_pool` (LVR),
/// and verifies that `LVR_realized / integral(V_t dt) = sigma^2 / 8` is invariant to `mu`.
pub fn misspec_drift(
    sigma: f64,
    mu: f64,
    t_yr: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let dt = t_yr / n_steps as f64;
    let drift = (mu - 0.5 * sigma * sigma) * dt;
    let vol_dt = sigma * dt.sqrt();
    let l = 0.5_f64;

    let rows = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x6003 ^ (mu.to_bits()), r as u64));
        let mut p = 1.0_f64;
        let mut v_rebal = 1.0_f64;
        let mut int_v = 0.0_f64;
        for _ in 0..n_steps {
            let v_t = pool_value(l, p);
            let d_t = pool_delta(l, p);
            int_v += v_t * dt;
            let p_next = p * (drift + vol_dt * rg.normal()).exp();
            v_rebal += d_t * (p_next - p);
            p = p_next;
        }
        let v_pool = pool_value(l, p);
        let v_hodl = 0.5 * p + 0.5;
        let lvr = v_rebal - v_pool;
        let il = v_pool - v_hodl;
        (v_pool - 1.0, il, lvr, int_v)
    });

    let mean_lp_ret = stats::mean(&rows.iter().map(|x| x.0).collect::<Vec<_>>());
    let mean_il = stats::mean(&rows.iter().map(|x| x.1).collect::<Vec<_>>());
    let mean_lvr = stats::mean(&rows.iter().map(|x| x.2).collect::<Vec<_>>());
    let mean_int_v = stats::mean(&rows.iter().map(|x| x.3).collect::<Vec<_>>());
    let lvr_per_v_yr = mean_lvr / mean_int_v;
    (mean_lp_ret, mean_il, mean_lvr / t_yr, lvr_per_v_yr)
}
