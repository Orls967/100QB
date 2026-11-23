use crate::linalg::solve_tridiagonal;
use crate::merton::{merton_call_series, MertonParams};
use crate::special::norm_cdf;

#[derive(Clone, Debug)]
pub struct PideCompareRow {
    pub k: f64,
    pub exact_call: f64,
    pub pide_call: f64,
    pub abs_err: f64,
    pub rel_err: f64,
}

#[derive(Clone, Debug)]
pub struct PideConvergenceRow {
    pub n_x: usize,
    pub n_t: usize,
    pub max_abs_err: f64,
    pub atm_abs_err: f64,
}

/// Solves the Merton PIDE for a European call option on a uniform log-spot grid `x = ln(S)`
/// using Crank-Nicolson with 4 Rannacher half-steps at `tau = 0` and predictor-corrector
/// jump integral evaluation with exact far-field asymptotic continuation.
pub fn price_call_pide_cn(p: &MertonParams, k: f64, t: f64, n_x: usize, n_t: usize) -> f64 {
    if t <= 0.0 {
        return (p.s0 - k).max(0.0);
    }
    let x0 = p.s0.ln();
    let x_half = 2.5_f64;
    let x_min = x0 - x_half;
    let x_max = x0 + x_half;
    let dx = (x_max - x_min) / (n_x - 1) as f64;

    let x_grid: Vec<f64> = (0..n_x).map(|i| x_min + i as f64 * dx).collect();
    let mut v: Vec<f64> = x_grid.iter().map(|&x| (x.exp() - k).max(0.0)).collect();

    let kappa = p.kappa();
    let sig2 = p.sigma * p.sigma;
    let drift = p.r - p.q - 0.5 * sig2 - p.lambda * kappa;
    let reaction = p.r + p.lambda;

    // Precompute discrete jump-size probabilities w[m] for y = m * dx, m in [-m_max..=m_max]
    let m_max = (((p.mu_j.abs() + 6.5 * p.delta_j) / dx).ceil() as usize).max(4);
    let mut jump_w = Vec::with_capacity(2 * m_max + 1);
    if p.lambda > 0.0 && p.delta_j > 0.0 {
        for m in -(m_max as isize)..=(m_max as isize) {
            let y_lo = (m as f64 - 0.5) * dx;
            let y_hi = (m as f64 + 0.5) * dx;
            let w = norm_cdf((y_hi - p.mu_j) / p.delta_j) - norm_cdf((y_lo - p.mu_j) / p.delta_j);
            jump_w.push((m, w));
        }
    }

    let eval_jump = |vals: &[f64], tau: f64| -> Vec<f64> {
        let mut out = vec![0.0_f64; n_x];
        if p.lambda <= 0.0 {
            return out;
        }
        let disc_q = (-p.q * tau).exp();
        let disc_r = (-p.r * tau).exp();
        for (i, slot) in out.iter_mut().enumerate().take(n_x - 1).skip(1) {
            let mut integral = 0.0_f64;
            for &(m, w) in &jump_w {
                if w == 0.0 {
                    continue;
                }
                let j = i as isize + m;
                let v_j = if j < 0 {
                    0.0
                } else if j >= n_x as isize {
                    let x_far = x_min + j as f64 * dx;
                    (x_far.exp() * disc_q - k * disc_r).max(0.0)
                } else {
                    vals[j as usize]
                };
                integral += w * v_j;
            }
            *slot = p.lambda * integral;
        }
        out
    };

    // Spatial operator coefficients: L V = alpha_l V_{i-1} + beta_l V_i + gamma_l V_{i+1}
    let alpha_l = 0.5 * sig2 / (dx * dx) - 0.5 * drift / dx;
    let beta_l = -sig2 / (dx * dx) - reaction;
    let gamma_l = 0.5 * sig2 / (dx * dx) + 0.5 * drift / dx;

    let n_int = n_x - 2;
    let dt_full = t / n_t as f64;

    // 4 Rannacher half-steps (theta = 1.0) for the first 2 full steps, then CN (theta = 0.5)
    let mut step_schedule = Vec::with_capacity(n_t + 2);
    for _ in 0..4 {
        step_schedule.push((0.5 * dt_full, 1.0_f64));
    }
    for _ in 2..n_t {
        step_schedule.push((dt_full, 0.5_f64));
    }

    let mut tau = 0.0_f64;
    for (dt, theta) in step_schedule {
        let tau_next = tau + dt;
        let bc_left = 0.0_f64;
        let bc_right =
            (x_max.exp() * (-p.q * tau_next).exp() - k * (-p.r * tau_next).exp()).max(0.0);

        let a_tri = vec![-theta * dt * alpha_l; n_int];
        let b_tri = vec![1.0 - theta * dt * beta_l; n_int];
        let c_tri = vec![-theta * dt * gamma_l; n_int];

        let j_curr = eval_jump(&v, tau);
        let mut rhs_base = vec![0.0_f64; n_int];
        for (idx, rb) in rhs_base.iter_mut().enumerate() {
            let i = idx + 1;
            let explicit_lv = alpha_l * v[i - 1] + beta_l * v[i] + gamma_l * v[i + 1];
            *rb = v[i] + (1.0 - theta) * dt * explicit_lv;
        }
        rhs_base[0] += theta * dt * alpha_l * bc_left;
        rhs_base[n_int - 1] += theta * dt * gamma_l * bc_right;

        // Predictor step with explicit jump at tau
        let mut rhs_pred = rhs_base.clone();
        for idx in 0..n_int {
            rhs_pred[idx] += dt * j_curr[idx + 1];
        }
        let u_pred = solve_tridiagonal(&a_tri, &b_tri, &c_tri, &rhs_pred);
        let mut v_pred = v.clone();
        v_pred[0] = bc_left;
        v_pred[1..(n_x - 1)].copy_from_slice(&u_pred);
        v_pred[n_x - 1] = bc_right;

        // Corrector step averaging jump at tau and tau_next (for CN theta = 0.5)
        let u_final = if p.lambda > 0.0 && theta < 0.99 {
            let j_next = eval_jump(&v_pred, tau_next);
            let mut rhs_corr = rhs_base;
            for idx in 0..n_int {
                rhs_corr[idx] += 0.5 * dt * (j_curr[idx + 1] + j_next[idx + 1]);
            }
            solve_tridiagonal(&a_tri, &b_tri, &c_tri, &rhs_corr)
        } else {
            u_pred
        };

        v[0] = bc_left;
        v[1..(n_x - 1)].copy_from_slice(&u_final);
        v[n_x - 1] = bc_right;
        tau = tau_next;
    }

    // x0 is exactly at the center grid point when n_x is odd, or interpolate linearly
    let pos = (x0 - x_min) / dx;
    let i0 = (pos.floor() as usize).min(n_x - 2);
    let frac = pos - i0 as f64;
    v[i0] * (1.0 - frac) + v[i0 + 1] * frac
}

pub fn compare_pide_vs_exact(
    p: &MertonParams,
    ks: &[f64],
    t: f64,
    n_x: usize,
    n_t: usize,
) -> Vec<PideCompareRow> {
    ks.iter()
        .map(|&k| {
            let exact_call = merton_call_series(p, k, t);
            let pide_call = price_call_pide_cn(p, k, t, n_x, n_t);
            let abs_err = (pide_call - exact_call).abs();
            let rel_err = abs_err / exact_call.max(1e-12);
            PideCompareRow {
                k,
                exact_call,
                pide_call,
                abs_err,
                rel_err,
            }
        })
        .collect()
}

pub fn pide_grid_convergence(
    p: &MertonParams,
    ks: &[f64],
    t: f64,
    grids: &[(usize, usize)],
) -> Vec<PideConvergenceRow> {
    grids
        .iter()
        .map(|&(n_x, n_t)| {
            let rows = compare_pide_vs_exact(p, ks, t, n_x, n_t);
            let max_abs_err = rows.iter().map(|r| r.abs_err).fold(0.0_f64, f64::max);
            let atm_abs_err = rows
                .iter()
                .find(|r| (r.k - p.s0).abs() < 1e-6)
                .map(|r| r.abs_err)
                .unwrap_or(max_abs_err);
            PideConvergenceRow {
                n_x,
                n_t,
                max_abs_err,
                atm_abs_err,
            }
        })
        .collect()
}
