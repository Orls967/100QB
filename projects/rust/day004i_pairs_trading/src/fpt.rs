//! First-passage time for zero-mean OU process:
//! (1) Exact closed-form (derived from Dynkin BVP: symmetric (pi/theta)*erfi(m*sqrt(theta)/sigma),
//!     plus general asymmetric erf-integral correction),
//! (2) Finite-difference BVP solver for E[tau] and E[tau^2] (hence Var(tau)),
//! (3) Fine-step Monte Carlo first-passage simulator with Brownian-bridge boundary correction.

use std::f64::consts::{FRAC_2_SQRT_PI, PI};

use crate::linalg::solve_tridiagonal;
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::special::erfi;
use crate::stats;

/// Exact expected first-passage time `E[tau(a -> m)]` (`a < m`) for zero-mean OU `dX = -theta*X dt + sigma dW`.
/// Derived from BVP `0.5*sigma^2 u'' - theta*x u' = -1`, `u(m)=0, u'(-inf)=0`:
/// `E[tau(a->m)] = (pi/(2*theta))*[erfi(zm) - erfi(za)] + (sqrt(pi)/theta) * int_{za}^{zm} exp(z^2)*erf(z) dz`.
/// For symmetric `a = -m`, the `erf` integral vanishes and `E[tau(-m->m)] = (pi/theta)*erfi(m*sqrt(theta)/sigma)`.
pub fn expected_fpt_exact(a: f64, m: f64, theta: f64, sigma: f64) -> f64 {
    assert!(a < m, "a must be < m");
    let scale = theta.sqrt() / sigma;
    let za = a * scale;
    let zm = m * scale;
    let term1 = (0.5 * PI / theta) * (erfi(zm) - erfi(za));
    if (a + m).abs() < 1e-13 * m.abs().max(1.0) {
        return term1;
    }
    // Composite Simpson for int_{za}^{zm} exp(z^2) * erf(z) dz
    let n = 400;
    let h = (zm - za) / n as f64;
    let erf_approx = |z: f64| -> f64 {
        // Taylor for |z|<=2.5, 1 - erfc otherwise (using series)
        let az = z.abs();
        if az <= 2.8 {
            let z2 = az * az;
            let mut t = az;
            let mut s = az;
            for k in 1..45 {
                t *= -z2 / k as f64;
                s += t / (2 * k + 1) as f64;
                if t.abs() < 1e-16 * s.abs() {
                    break;
                }
            }
            z.signum() * FRAC_2_SQRT_PI * s
        } else {
            z.signum() * 1.0
        }
    };
    let integrand = |z: f64| (z * z).exp() * erf_approx(z);
    let mut acc = integrand(za) + integrand(zm);
    for i in 1..n {
        let z = za + i as f64 * h;
        acc += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(z);
    }
    term1 + (PI.sqrt() / theta) * (h / 3.0) * acc
}

/// Solves BVP `0.5*sigma^2 u'' - theta*x u' = -1` on `[x_min, m]` with `u'(x_min) = 0, u(m) = 0`,
/// and second-moment BVP `0.5*sigma^2 v'' - theta*x v' = -2*u(x)` to return `(E[tau(a->m)], Var(tau(a->m)))`.
pub fn solve_fpt_bvp(a: f64, m: f64, theta: f64, sigma: f64, n_grid: usize) -> (f64, f64) {
    let stat_sd = sigma / (2.0 * theta).sqrt();
    let x_min = (a - 6.0 * stat_sd).min(-8.0 * stat_sd);
    let h = (m - x_min) / n_grid as f64;
    let diff = 0.5 * sigma * sigma;

    // Unknowns u[0..n_grid] at x_i = x_min + i*h; u[n_grid] = 0 (Dirichlet at m).
    // At i = 0: Neumann u'(x_min) = 0 => (u_1 - u_{-1})/(2h) = 0 => u_{-1} = u_1
    // Central difference: diff*(u_{i+1} - 2u_i + u_{i-1})/h^2 - theta*x_i*(u_{i+1} - u_{i-1})/(2h) = -1
    let mut sub = vec![0.0_f64; n_grid];
    let mut diag = vec![0.0_f64; n_grid];
    let mut sup = vec![0.0_f64; n_grid];
    let mut rhs = vec![-1.0_f64; n_grid];

    for i in 0..n_grid {
        let x = x_min + i as f64 * h;
        let c_m = diff / (h * h) + theta * x / (2.0 * h);
        let c_0 = -2.0 * diff / (h * h);
        let c_p = diff / (h * h) - theta * x / (2.0 * h);
        if i == 0 {
            diag[0] = c_0;
            sup[0] = c_m + c_p;
        } else {
            sub[i] = c_m;
            diag[i] = c_0;
            if i + 1 < n_grid {
                sup[i] = c_p;
            }
        }
    }

    let u = solve_tridiagonal(&sub, &diag, &sup, &rhs);
    for i in 0..n_grid {
        rhs[i] = -2.0 * u[i];
    }
    let v = solve_tridiagonal(&sub, &diag, &sup, &rhs);

    let pos = (a - x_min) / h;
    let i0 = (pos.floor() as usize).min(n_grid - 1);
    let i1 = (i0 + 1).min(n_grid);
    let w = pos - i0 as f64;
    let u1 = if i1 == n_grid { 0.0 } else { u[i1] };
    let v1 = if i1 == n_grid { 0.0 } else { v[i1] };
    let e_tau = u[i0] * (1.0 - w) + u1 * w;
    let e_tau2 = v[i0] * (1.0 - w) + v1 * w;
    (e_tau, (e_tau2 - e_tau * e_tau).max(0.0))
}

#[derive(Clone, Debug)]
pub struct FptMcResult {
    pub mean_tau: f64,
    pub se_tau: f64,
    pub sd_tau: f64,
}

/// Fine-step Monte Carlo first-passage time simulator with Brownian-bridge overshoot/exit correction.
#[allow(clippy::too_many_arguments)]
pub fn simulate_fpt_mc(
    a: f64,
    m: f64,
    theta: f64,
    sigma: f64,
    dt: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    salt: u64,
) -> FptMcResult {
    let phi = (-theta * dt).exp();
    let sd_step = (sigma * sigma * (1.0 - phi * phi) / (2.0 * theta)).sqrt();
    let taus = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, salt, r as u64));
        let mut x = a;
        let mut t = 0.0_f64;
        loop {
            let x_next = x * phi + sd_step * rg.normal();
            if x_next >= m {
                // Linear interpolation fraction within step
                let frac = ((m - x) / (x_next - x)).clamp(0.0, 1.0);
                t += frac * dt;
                break;
            }
            // Brownian-bridge intra-step hitting probability: exp(-2*(m - x)*(m - x_next)/(sigma^2 * dt))
            let p_hit = (-2.0 * (m - x) * (m - x_next) / (sigma * sigma * dt)).exp();
            if rg.uniform() < p_hit {
                t += 0.5 * dt;
                break;
            }
            x = x_next;
            t += dt;
        }
        t
    });
    FptMcResult {
        mean_tau: stats::mean(&taus),
        se_tau: stats::se(&taus),
        sd_tau: stats::sd(&taus),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fpt_closed_form_matches_bvp_symmetric_and_asymmetric() {
        let theta = 16.0_f64;
        let sigma = 0.20_f64;
        let stat_sd = sigma / (2.0 * theta).sqrt();
        for &(a, m) in &[
            (-stat_sd, stat_sd),
            (-1.5 * stat_sd, 1.5 * stat_sd),
            (0.0, stat_sd),
        ] {
            let ex = expected_fpt_exact(a, m, theta, sigma);
            let (bvp, _) = solve_fpt_bvp(a, m, theta, sigma, 2000);
            let rel = (ex - bvp).abs() / ex;
            assert!(rel < 1e-3, "a={a} m={m} ex={ex} bvp={bvp} rel={rel}");
        }
    }
}
