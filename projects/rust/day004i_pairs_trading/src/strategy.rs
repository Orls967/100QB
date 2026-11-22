//! Bertram optimal symmetric OU trading threshold `a*` (for `m = -a`, `a < 0`):
//! `R(a) = (-2a - c) / (2 * E[tau(a -> -a)])`, Sharpe criterion, root-finding of `R'(a) = 0`,
//! 200-year MC verification with batch means, and sensitivity `a*(theta, c)`.

use std::f64::consts::PI;

use crate::fpt::{expected_fpt_exact, solve_fpt_bvp};
use crate::ou::{simulate_ou_path, OuParams};
use crate::rng::Rng;
use crate::special::{erfi, erfi_deriv};
use crate::stats;

/// Expected return rate per year `R(a) = (-2a - c) / (2 * E[tau(a -> -a)])` for `a < -c/2`.
#[inline]
pub fn expected_return_rate(a: f64, theta: f64, sigma: f64, c: f64) -> f64 {
    if -2.0 * a <= c {
        return 0.0;
    }
    let e_tau = expected_fpt_exact(a, -a, theta, sigma);
    (-2.0 * a - c) / (2.0 * e_tau)
}

/// Finds optimal `a* < -c/2` maximizing `R(a)` via golden-section search.
pub fn optimal_threshold_direct(theta: f64, sigma: f64, c: f64) -> (f64, f64) {
    let stat_sd = sigma / (2.0 * theta).sqrt();
    let mut lo = -4.0 * stat_sd;
    let mut hi = -0.5 * c - 1e-7;
    if hi <= lo {
        return (hi, 0.0);
    }
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut f1 = expected_return_rate(x1, theta, sigma, c);
    let mut f2 = expected_return_rate(x2, theta, sigma, c);
    for _ in 0..80 {
        if (hi - lo).abs() < 1e-13 {
            break;
        }
        if f1 < f2 {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = expected_return_rate(x2, theta, sigma, c);
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = expected_return_rate(x1, theta, sigma, c);
        }
    }
    let a_star = 0.5 * (lo + hi);
    (a_star, expected_return_rate(a_star, theta, sigma, c))
}

/// Finds optimal `a*` by solving the stationary condition `R'(a) = 0`:
/// `2 * erfi(z) - (2*z - c*sqrt(theta)/sigma) * erfi'(z) = 0`, where `z = -a * sqrt(theta)/sigma > 0`.
pub fn optimal_threshold_foc(theta: f64, sigma: f64, c: f64) -> f64 {
    let scale = theta.sqrt() / sigma;
    let c_scaled = c * scale;
    let f = |z: f64| -> f64 { 2.0 * erfi(z) - (2.0 * z - c_scaled) * erfi_deriv(z) };
    let mut lo = 0.5 * c_scaled + 1e-8;
    let mut hi = 4.0_f64;
    while f(hi) > 0.0 && hi < 8.0 {
        hi += 1.0;
    }
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() < 1e-14 {
            break;
        }
    }
    let z_star = 0.5 * (lo + hi);
    -z_star / scale
}

/// Analytical/BVP Sharpe rate `(R(a) * sqrt(2*E[tau])) / sqrt(Var(tau_cycle))` using `solve_fpt_bvp`.
pub fn theoretical_sharpe_at_a(a: f64, theta: f64, sigma: f64, c: f64) -> f64 {
    let (e_tau, var_tau) = solve_fpt_bvp(a, -a, theta, sigma, 1200);
    let profit_per_cycle = -2.0 * a - c;
    // Per unit time (year): mean rate = profit / (2*e_tau), var rate = profit^2 * (2*var_tau) / (2*e_tau)^3
    // So Sharpe = (mean_rate) / sqrt(var_rate) = sqrt(2*e_tau) / sqrt(2*var_tau) * ... wait:
    // Renewal reward variance rate = E[R_cycle]^2 * Var(T_cycle) / E[T_cycle]^3 => SD/yr = |profit| * sqrt(2*var_tau) / (2*e_tau)^{3/2}
    // Thus Sharpe = (profit / (2*e_tau)) / (profit * sqrt(2*var_tau) / (2*e_tau)^{3/2}) = sqrt(e_tau / var_tau)!
    let _ = profit_per_cycle;
    (e_tau / var_tau.max(1e-12)).sqrt()
}

#[derive(Clone, Debug)]
pub struct StrategyMcCheck {
    pub total_years: usize,
    pub trades_per_year_mc: f64,
    pub se_trades_per_year: f64,
    pub trades_per_year_exact: f64,
    pub profit_per_year_mc: f64,
    pub se_profit_per_year: f64,
    pub profit_per_year_exact: f64,
    pub z_trades: f64,
    pub z_profit: f64,
}

/// Simulates 200 years (200 x 1-year batches) of OU trading with threshold `a` (`m = -a`, cost `c` per round-trip)
/// on fine sub-steps (`sub = 10` per day -> `dt_fine = 1/2520`) to verify `R(a)` and `1/(2*E[tau])` within 3 SE.
pub fn verify_strategy_mc_200yr(
    theta: f64,
    mu: f64,
    sigma: f64,
    a: f64,
    c: f64,
    seed: u64,
) -> StrategyMcCheck {
    let total_years = 200_usize;
    let sub = 10_usize;
    let steps_per_yr = 252 * sub;
    let p_fine = OuParams {
        theta,
        mu,
        sigma,
        dt: 1.0 / (252 * sub) as f64,
    };
    let mut rg = Rng::new(seed);
    let path = simulate_ou_path(p_fine, total_years * steps_per_yr, true, &mut rg);

    let mut state_wait_sell = false;
    let mut yr_profits = vec![0.0_f64; total_years];
    let mut yr_trades = vec![0.0_f64; total_years];
    let cycle_profit = -2.0 * a - c;
    let var_dt = sigma * sigma * p_fine.dt;

    for yr in 0..total_years {
        let s_idx = yr * steps_per_yr;
        let e_idx = (yr + 1) * steps_per_yr;
        for idx in s_idx..e_idx {
            let x_prev = path[idx];
            let x = path[idx + 1];
            if !state_wait_sell {
                let b = mu + a;
                let hit = x <= b
                    || (x_prev > b
                        && rg.uniform() < (-2.0 * (x_prev - b) * (x - b) / var_dt).exp());
                if hit {
                    state_wait_sell = true;
                }
            } else {
                let b = mu - a;
                let hit = x >= b
                    || (x_prev < b
                        && rg.uniform() < (-2.0 * (b - x_prev) * (b - x) / var_dt).exp());
                if hit {
                    state_wait_sell = false;
                    yr_trades[yr] += 1.0;
                    yr_profits[yr] += cycle_profit;
                }
            }
        }
    }

    let e_tau = (PI / theta) * erfi(-a * theta.sqrt() / sigma);
    let exact_tr = 1.0 / (2.0 * e_tau);
    let exact_pr = cycle_profit / (2.0 * e_tau);
    let m_tr = stats::mean(&yr_trades);
    let se_tr = stats::se(&yr_trades);
    let m_pr = stats::mean(&yr_profits);
    let se_pr = stats::se(&yr_profits);

    StrategyMcCheck {
        total_years,
        trades_per_year_mc: m_tr,
        se_trades_per_year: se_tr,
        trades_per_year_exact: exact_tr,
        profit_per_year_mc: m_pr,
        se_profit_per_year: se_pr,
        profit_per_year_exact: exact_pr,
        z_trades: (m_tr - exact_tr).abs() / se_tr.max(1e-12),
        z_profit: (m_pr - exact_pr).abs() / se_pr.max(1e-12),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optimal_threshold_direct_matches_foc() {
        for &theta in &[8.0, 16.0, 32.0] {
            for &c in &[0.001, 0.002, 0.005] {
                let (a_dir, _) = optimal_threshold_direct(theta, 0.20, c);
                let a_foc = optimal_threshold_foc(theta, 0.20, c);
                assert!(
                    (a_dir - a_foc).abs() < 1e-7,
                    "theta={theta} c={c}: {a_dir} vs {a_foc}"
                );
            }
        }
    }
}
