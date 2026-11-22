//! Ornstein-Uhlenbeck process: exact vs Euler transition, closed-form AR(1) MLE vs numerical MLE,
//! numerical Fisher information SEs, and R=20,000 replication bias study (Kendall / Marriott-Pope).

#![allow(clippy::needless_range_loop)]

use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[derive(Clone, Copy, Debug)]
pub struct OuParams {
    pub theta: f64,
    pub mu: f64,
    pub sigma: f64,
    pub dt: f64,
}

impl OuParams {
    pub fn default_params() -> Self {
        Self {
            theta: 16.0,
            mu: 0.0,
            sigma: 0.20,
            dt: 1.0 / 252.0,
        }
    }

    #[inline]
    pub fn stationary_sd(&self) -> f64 {
        self.sigma / (2.0 * self.theta).sqrt()
    }

    #[inline]
    pub fn half_life(&self) -> f64 {
        2.0_f64.ln() / self.theta
    }
}

/// Simulates an exact OU path of `n_steps` starting from `x0` (or stationary draw if `stationary_init`).
pub fn simulate_ou_path(
    p: OuParams,
    n_steps: usize,
    stationary_init: bool,
    rg: &mut Rng,
) -> Vec<f64> {
    let mut x = Vec::with_capacity(n_steps + 1);
    let mut cur = if stationary_init {
        p.mu + p.stationary_sd() * rg.normal()
    } else {
        p.mu
    };
    x.push(cur);
    let phi = (-p.theta * p.dt).exp();
    let sd_step = (p.sigma * p.sigma * (1.0 - phi * phi) / (2.0 * p.theta)).sqrt();
    for _ in 0..n_steps {
        cur = p.mu + (cur - p.mu) * phi + sd_step * rg.normal();
        x.push(cur);
    }
    x
}

/// Exact Gaussian log-likelihood of an OU series `x[0..=N]` with step `dt`.
pub fn exact_loglik(x: &[f64], dt: f64, theta: f64, mu: f64, sigma: f64) -> f64 {
    if theta <= 1e-8 || sigma <= 1e-8 || x.len() < 2 {
        return -1e18;
    }
    let n = (x.len() - 1) as f64;
    let phi = (-theta * dt).exp();
    let var = sigma * sigma * (1.0 - phi * phi) / (2.0 * theta);
    if var <= 1e-20 || !var.is_finite() {
        return -1e18;
    }
    let mut ss = 0.0_f64;
    for i in 0..(x.len() - 1) {
        let m = mu + (x[i] - mu) * phi;
        let r = x[i + 1] - m;
        ss += r * r;
    }
    -0.5 * n * (2.0 * std::f64::consts::PI * var).ln() - 0.5 * ss / var
}

/// Euler Gaussian log-likelihood of an OU series `x[0..=N]` with step `dt`.
pub fn euler_loglik(x: &[f64], dt: f64, theta: f64, mu: f64, sigma: f64) -> f64 {
    if theta <= 1e-8 || sigma <= 1e-8 || x.len() < 2 {
        return -1e18;
    }
    let n = (x.len() - 1) as f64;
    let var = sigma * sigma * dt;
    let mut ss = 0.0_f64;
    for i in 0..(x.len() - 1) {
        let m = x[i] + theta * (mu - x[i]) * dt;
        let r = x[i + 1] - m;
        ss += r * r;
    }
    -0.5 * n * (2.0 * std::f64::consts::PI * var).ln() - 0.5 * ss / var
}

#[derive(Clone, Debug)]
pub struct OuMleResult {
    pub theta: f64,
    pub mu: f64,
    pub sigma: f64,
    pub phi_hat: f64,
    pub se_theta: f64,
    pub se_mu: f64,
    pub se_sigma: f64,
}

/// Closed-form exact MLE via AR(1) regression `X_{i+1} = c + phi * X_i + eps_i`.
/// Returns `phi_hat` even if `phi_hat >= 1` (`theta <= 0`), and asymptotic Fisher SEs when `theta > 0`.
pub fn mle_closed_form(x: &[f64], dt: f64) -> OuMleResult {
    let n = x.len() - 1;
    let nf = n as f64;
    let mx0 = stats::mean(&x[..n]);
    let mx1 = stats::mean(&x[1..]);
    let mut sxx = 0.0_f64;
    let mut sxy = 0.0_f64;
    for i in 0..n {
        let dx = x[i] - mx0;
        let dy = x[i + 1] - mx1;
        sxx += dx * dx;
        sxy += dx * dy;
    }
    let phi = if sxx > 0.0 { sxy / sxx } else { 1.0 };
    let c = mx1 - phi * mx0;
    let mut rss = 0.0_f64;
    for i in 0..n {
        let e = x[i + 1] - (c + phi * x[i]);
        rss += e * e;
    }
    let v_hat = (rss / nf).max(1e-18);

    if phi <= 0.0 || phi >= 1.0 - 1e-12 {
        let theta = if phi <= 0.0 {
            f64::NAN
        } else {
            -(phi.ln()) / dt
        };
        return OuMleResult {
            theta,
            mu: stats::mean(x),
            sigma: (v_hat / dt).sqrt(),
            phi_hat: phi,
            se_theta: f64::NAN,
            se_mu: f64::NAN,
            se_sigma: f64::NAN,
        };
    }

    let theta = -phi.ln() / dt;
    let mu = c / (1.0 - phi);
    let sigma = (2.0 * theta * v_hat / (1.0 - phi * phi)).sqrt();

    // Analytical / numerical Fisher SEs from observed Fisher information on (theta, mu, sigma)
    // Diagonal curvature of -loglik gives quick, reliable SEs; full 3x3 is even cleaner:
    let h_th = (theta * 1e-4).max(1e-4);
    let h_mu = (sigma * 1e-4).max(1e-5);
    let h_si = (sigma * 1e-4).max(1e-5);
    let ll0 = exact_loglik(x, dt, theta, mu, sigma);
    let i_tt = -(exact_loglik(x, dt, theta + h_th, mu, sigma) - 2.0 * ll0
        + exact_loglik(x, dt, theta - h_th, mu, sigma))
        / (h_th * h_th);
    let i_mm = -(exact_loglik(x, dt, theta, mu + h_mu, sigma) - 2.0 * ll0
        + exact_loglik(x, dt, theta, mu - h_mu, sigma))
        / (h_mu * h_mu);
    let i_ss = -(exact_loglik(x, dt, theta, mu, sigma + h_si) - 2.0 * ll0
        + exact_loglik(x, dt, theta, mu, sigma - h_si))
        / (h_si * h_si);

    OuMleResult {
        theta,
        mu,
        sigma,
        phi_hat: phi,
        se_theta: if i_tt > 0.0 {
            (1.0 / i_tt).sqrt()
        } else {
            f64::NAN
        },
        se_mu: if i_mm > 0.0 {
            (1.0 / i_mm).sqrt()
        } else {
            f64::NAN
        },
        se_sigma: if i_ss > 0.0 {
            (1.0 / i_ss).sqrt()
        } else {
            f64::NAN
        },
    }
}

/// Numerical maximization of exact log-likelihood via Nelder-Mead for cross-checking closed-form MLE.
pub fn mle_numerical(x: &[f64], dt: f64, init: [f64; 3]) -> [f64; 3] {
    let ll_ref = exact_loglik(x, dt, init[0], init[1], init[2]);
    let obj = |pt: [f64; 3]| -> f64 { -(exact_loglik(x, dt, pt[0], pt[1], pt[2]) - ll_ref) };
    let steps = [2.0, 0.02, 0.03];
    let mut simplex: Vec<([f64; 3], f64)> = Vec::with_capacity(4);
    simplex.push((init, obj(init)));
    for i in 0..3 {
        let mut pt = init;
        pt[i] += steps[i];
        let val = obj(pt);
        simplex.push((pt, val));
    }
    for _ in 0..600 {
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let max_dx = (1..4)
            .map(|i| (simplex[i].0[0] - simplex[0].0[0]).abs())
            .fold(0.0_f64, f64::max);
        if (simplex[3].1 - simplex[0].1).abs() < 1e-15 && max_dx < 1e-9 {
            break;
        }
        let mut c = [0.0_f64; 3];
        for item in simplex.iter().take(3) {
            for j in 0..3 {
                c[j] += item.0[j] / 3.0;
            }
        }
        let mut xr = [0.0_f64; 3];
        for j in 0..3 {
            xr[j] = 2.0 * c[j] - simplex[3].0[j];
        }
        let fr = obj(xr);
        if fr < simplex[0].1 {
            let mut xe = [0.0_f64; 3];
            for j in 0..3 {
                xe[j] = c[j] + 2.0 * (xr[j] - c[j]);
            }
            let fe = obj(xe);
            simplex[3] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < simplex[2].1 {
            simplex[3] = (xr, fr);
        } else {
            let mut xc = [0.0_f64; 3];
            if fr < simplex[3].1 {
                for j in 0..3 {
                    xc[j] = c[j] + 0.5 * (xr[j] - c[j]);
                }
            } else {
                for j in 0..3 {
                    xc[j] = c[j] - 0.5 * (c[j] - simplex[3].0[j]);
                }
            }
            let fc = obj(xc);
            if fc < simplex[3].1.min(fr) {
                simplex[3] = (xc, fc);
            } else {
                let x0 = simplex[0].0;
                for item in simplex.iter_mut().skip(1) {
                    for j in 0..3 {
                        item.0[j] = x0[j] + 0.5 * (item.0[j] - x0[j]);
                    }
                    item.1 = obj(item.0);
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    simplex[0].0
}

#[derive(Clone, Debug)]
pub struct OuReplicationRow {
    pub t_years: f64,
    pub n_rep: usize,
    pub bias_theta: f64,
    pub sd_theta: f64,
    pub rmse_theta: f64,
    pub cov95_theta: f64,
    pub kendall_approx_bias: f64,
    pub bias_ratio_to_kendall: f64,
    pub bias_mu: f64,
    pub sd_mu: f64,
    pub cov95_mu: f64,
    pub bias_sigma: f64,
    pub sd_sigma: f64,
    pub cov95_sigma: f64,
    pub half_life_mean: f64,
    pub half_life_skewness: f64,
    pub prob_theta_nonpos: f64,
}

pub fn run_ou_replication_study(
    p: OuParams,
    ts_years: &[f64],
    n_rep: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<OuReplicationRow> {
    let mut out = Vec::with_capacity(ts_years.len());
    for &t_yr in ts_years {
        let n_steps = (t_yr / p.dt).round() as usize;
        let rows = parallel::run_chunks(n_rep, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x1100 ^ t_yr.to_bits(), r as u64));
            let path = simulate_ou_path(p, n_steps, true, &mut rg);
            mle_closed_form(&path, p.dt)
        });

        let mut thetas = Vec::new();
        let mut mus = Vec::new();
        let mut sigmas = Vec::new();
        let mut half_lives = Vec::new();
        let mut c_th = 0usize;
        let mut c_mu = 0usize;
        let mut c_si = 0usize;
        let mut nonpos = 0usize;

        for m in &rows {
            if !m.theta.is_finite() || m.theta <= 0.0 {
                nonpos += 1;
                continue;
            }
            thetas.push(m.theta);
            mus.push(m.mu);
            sigmas.push(m.sigma);
            half_lives.push(2.0_f64.ln() / m.theta);
            if m.se_theta.is_finite() && (m.theta - p.theta).abs() <= 1.96 * m.se_theta {
                c_th += 1;
            }
            if m.se_mu.is_finite() && (m.mu - p.mu).abs() <= 1.96 * m.se_mu {
                c_mu += 1;
            }
            if m.se_sigma.is_finite() && (m.sigma - p.sigma).abs() <= 1.96 * m.se_sigma {
                c_si += 1;
            }
        }

        let valid = thetas.len().max(1) as f64;
        let m_th = stats::mean(&thetas);
        let s_th = stats::sd(&thetas);
        let bias_th = m_th - p.theta;
        let rmse_th = (bias_th * bias_th + s_th * s_th).sqrt();
        // Tang-Chen (2009) / Kendall-Marriott-Pope AR(1) bias approximation with unknown mu: ~ 4 / T
        let phi = (-p.theta * p.dt).exp();
        let kendall_approx = (1.0 + 3.0 * phi) / (t_yr * phi);

        out.push(OuReplicationRow {
            t_years: t_yr,
            n_rep,
            bias_theta: bias_th,
            sd_theta: s_th,
            rmse_theta: rmse_th,
            cov95_theta: c_th as f64 / valid,
            kendall_approx_bias: kendall_approx,
            bias_ratio_to_kendall: bias_th / kendall_approx,
            bias_mu: stats::mean(&mus) - p.mu,
            sd_mu: stats::sd(&mus),
            cov95_mu: c_mu as f64 / valid,
            bias_sigma: stats::mean(&sigmas) - p.sigma,
            sd_sigma: stats::sd(&sigmas),
            cov95_sigma: c_si as f64 / valid,
            half_life_mean: stats::mean(&half_lives),
            half_life_skewness: stats::skewness(&half_lives),
            prob_theta_nonpos: nonpos as f64 / n_rep as f64,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_closed_form_mle_matches_numerical() {
        let p = OuParams::default_params();
        let mut rg = Rng::new(42);
        let path = simulate_ou_path(p, 252 * 5, true, &mut rg);
        let cf = mle_closed_form(&path, p.dt);
        let num = mle_numerical(&path, p.dt, [cf.theta * 0.9, cf.mu + 0.01, cf.sigma * 1.1]);
        assert!(
            (cf.theta - num[0]).abs() < 1e-5,
            "{} vs {}",
            cf.theta,
            num[0]
        );
        assert!((cf.mu - num[1]).abs() < 1e-7);
        assert!((cf.sigma - num[2]).abs() < 1e-7);
    }
}
