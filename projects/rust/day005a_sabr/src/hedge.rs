#![allow(clippy::too_many_arguments)]

use crate::black::black_delta;
use crate::cev::cev_call_exact;
use crate::hagan::{hagan_call, hagan_implied_vol};
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[inline]
pub fn sabr_naive_delta(f: f64, k: f64, tau: f64, alpha: f64, beta: f64, rho: f64, nu: f64) -> f64 {
    if tau <= 1e-8 || f <= 1e-8 {
        return if f > k { 1.0 } else { 0.0 };
    }
    let h = (1e-4 * f).max(1e-7);
    let cp = hagan_call(f + h, k, tau, alpha, beta, rho, nu);
    let cm = hagan_call((f - h).max(1e-8), k, tau, alpha, beta, rho, nu);
    (cp - cm) / (f + h - (f - h).max(1e-8))
}

#[inline]
pub fn sabr_bartlett_delta(
    f: f64,
    k: f64,
    tau: f64,
    alpha: f64,
    beta: f64,
    rho: f64,
    nu: f64,
) -> f64 {
    if tau <= 1e-8 || f <= 1e-8 {
        return if f > k { 1.0 } else { 0.0 };
    }
    let d_f = sabr_naive_delta(f, k, tau, alpha, beta, rho, nu);
    let ha = (1e-4 * alpha).max(1e-7);
    let cap = hagan_call(f, k, tau, alpha + ha, beta, rho, nu);
    let cam = hagan_call(f, k, tau, (alpha - ha).max(1e-7), beta, rho, nu);
    let d_alpha = (cap - cam) / (alpha + ha - (alpha - ha).max(1e-7));
    d_f + d_alpha * (rho * nu / f.powf(beta))
}

#[inline]
pub fn cev_delta_fd(f: f64, k: f64, tau: f64, alpha: f64, beta: f64) -> f64 {
    if tau <= 1e-8 || f <= 1e-8 {
        return if f > k { 1.0 } else { 0.0 };
    }
    let h = (1e-4 * f).max(1e-7);
    let cp = cev_call_exact(f + h, k, alpha, beta, tau);
    let cm = cev_call_exact((f - h).max(1e-8), k, alpha, beta, tau);
    (cp - cm) / (f + h - (f - h).max(1e-8))
}

#[derive(Clone, Debug)]
pub struct HedgeComparison {
    pub initial_price_bps: f64,
    pub mean_pnl_bps: [f64; 4], // Bartlett, Naive, BlackATM, CEV
    pub sd_pnl_bps: [f64; 4],
    pub paired_sd_diff_bartlett_minus_naive_bps: f64,
    pub paired_se_diff_bps: f64,
    pub paired_var_diff_bps2: f64,
    pub paired_var_diff_se_bps2: f64,
    pub pnl_bartlett_bps: Vec<f64>,
    pub pnl_naive_bps: Vec<f64>,
}

/// Simulates daily delta-hedging (`n_steps = 252`) of a short ATM call (`K = F0`, `T = 1`)
/// under the true SABR DGP `(alpha_true, beta_true, rho_true, nu_true)`, hedging with
/// calibrated model `(alpha_cal, beta_cal, rho_cal, nu_cal)`.
/// All PnL metrics are returned in basis points of notional (`* 1e4`).
pub fn run_sabr_hedging_crn(
    f0: f64,
    k: f64,
    t: f64,
    alpha_true: f64,
    beta_true: f64,
    rho_true: f64,
    nu_true: f64,
    alpha_cal: f64,
    beta_cal: f64,
    rho_cal: f64,
    nu_cal: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> HedgeComparison {
    let dt = t / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let drift_a = -0.5 * nu_true * nu_true * dt;
    let vol_a = nu_true * sqrt_dt;
    let one_m_r2 = (1.0 - rho_true * rho_true).max(0.0).sqrt();

    let c0 = hagan_call(f0, k, t, alpha_cal, beta_cal, rho_cal, nu_cal);
    let atm_iv0 = hagan_implied_vol(f0, f0, t, alpha_cal, beta_cal, rho_cal, nu_cal);
    // For CEV benchmark with calibrated beta_cal, match initial instantaneous vol:
    // alpha_cev * F0^beta_cal = alpha_cal * F0^beta_cal => alpha_cev = alpha_cal
    let alpha_cev = alpha_cal;

    let rows = run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x5001, r as u64));
        let mut f = f0;
        let mut a_true = alpha_true;
        // Scale calibrated alpha proportionally to true alpha change to reflect vol state
        let a_scale = alpha_cal / alpha_true;

        let mut g_bart = 0.0_f64;
        let mut g_naiv = 0.0_f64;
        let mut g_blk = 0.0_f64;
        let mut g_cev = 0.0_f64;

        for step in 0..n_steps {
            let tau = t - step as f64 * dt;
            let a_cal_t = a_true * a_scale;
            let db = if f > 0.0 {
                sabr_bartlett_delta(f, k, tau, a_cal_t, beta_cal, rho_cal, nu_cal)
            } else {
                0.0
            };
            let dn = if f > 0.0 {
                sabr_naive_delta(f, k, tau, a_cal_t, beta_cal, rho_cal, nu_cal)
            } else {
                0.0
            };
            let dk = if f > 0.0 {
                black_delta(f, k, atm_iv0, tau)
            } else {
                0.0
            };
            let dc = if f > 0.0 {
                cev_delta_fd(f, k, tau, alpha_cev, beta_cal)
            } else {
                0.0
            };

            let (z1, z2_ind) = rg.normal_pair();
            let dw2 = rho_true * z1 + one_m_r2 * z2_ind;
            let f_next = if f > 0.0 {
                (f + a_true * f.powf(beta_true) * sqrt_dt * z1).max(0.0)
            } else {
                0.0
            };
            let df = f_next - f;
            g_bart += db * df;
            g_naiv += dn * df;
            g_blk += dk * df;
            g_cev += dc * df;

            f = f_next;
            a_true *= (drift_a + vol_a * dw2).exp();
        }
        let pay = (f - k).max(0.0);
        let to_bps = 1e4;
        [
            (c0 + g_bart - pay) * to_bps,
            (c0 + g_naiv - pay) * to_bps,
            (c0 + g_blk - pay) * to_bps,
            (c0 + g_cev - pay) * to_bps,
        ]
    });

    let p0: Vec<f64> = rows.iter().map(|x| x[0]).collect();
    let p1: Vec<f64> = rows.iter().map(|x| x[1]).collect();
    let p2: Vec<f64> = rows.iter().map(|x| x[2]).collect();
    let p3: Vec<f64> = rows.iter().map(|x| x[3]).collect();

    let diff_sq: Vec<f64> = p0
        .iter()
        .zip(p1.iter())
        .map(|(&a, &b)| a * a - b * b)
        .collect();
    let diff_pnl: Vec<f64> = p0.iter().zip(p1.iter()).map(|(&a, &b)| a - b).collect();

    HedgeComparison {
        initial_price_bps: c0 * 1e4,
        mean_pnl_bps: [
            stats::mean(&p0),
            stats::mean(&p1),
            stats::mean(&p2),
            stats::mean(&p3),
        ],
        sd_pnl_bps: [
            stats::sd(&p0),
            stats::sd(&p1),
            stats::sd(&p2),
            stats::sd(&p3),
        ],
        paired_sd_diff_bartlett_minus_naive_bps: stats::sd(&p0) - stats::sd(&p1),
        paired_se_diff_bps: stats::se(&diff_pnl),
        paired_var_diff_bps2: stats::mean(&diff_sq),
        paired_var_diff_se_bps2: stats::se(&diff_sq),
        pnl_bartlett_bps: p0,
        pnl_naive_bps: p1,
    }
}
