//! Direct nonlinear theta-form HJB backward finite-difference solver (without v-transform)
//! for independent verification and dt convergence-order testing.

use crate::glt::GltParams;

fn golden_max<F>(mut f: F, mut a: f64, mut b: f64) -> (f64, f64)
where
    F: FnMut(f64) -> f64,
{
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut c = b - phi * (b - a);
    let mut d = a + phi * (b - a);
    let mut fc = f(c);
    let mut fd = f(d);
    for _ in 0..55 {
        if (b - a).abs() < 1e-13 {
            break;
        }
        if fc < fd {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = f(d);
        } else {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = f(c);
        }
    }
    let x = 0.5 * (a + b);
    (x, f(x))
}

/// Solves nonlinear theta-form HJB directly backward from `t=T` to `t=0` with step `dt`.
/// Optional `eps_adverse` shifts execution mark-to-market by `-(q+1)*eps` (bid) and `+(q-1)*eps` (ask).
pub fn solve_nonlinear_hjb_theta(
    p: &GltParams,
    dt: f64,
    eps_adverse: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = p.dim();
    let n_steps = (p.t_horizon / dt).round() as usize;
    let actual_dt = p.t_horizon / n_steps as f64;
    let mut theta = vec![0.0; n];
    let mut next_theta = vec![0.0; n];
    let mut opt_bid = vec![f64::INFINITY; n];
    let mut opt_ask = vec![f64::INFINITY; n];

    for _step in 0..n_steps {
        for i in 0..n {
            let q = i as isize - p.q_max;
            let drift = -0.5 * p.gamma * p.sigma * p.sigma * (q * q) as f64;

            let mut fb_val = 0.0;
            if i + 1 < n {
                let d_th_b = theta[i + 1] - theta[i] - (q + 1) as f64 * eps_adverse;
                let (db, val) = golden_max(
                    |d| {
                        (p.a_intens / p.gamma)
                            * (-p.k * d).exp()
                            * (1.0 - (-p.gamma * (d + d_th_b)).exp())
                    },
                    -0.5,
                    4.0,
                );
                opt_bid[i] = db;
                fb_val = val;
            }

            let mut fa_val = 0.0;
            if i > 0 {
                let d_th_a = theta[i - 1] - theta[i] + (q - 1) as f64 * eps_adverse;
                let (da, val) = golden_max(
                    |d| {
                        (p.a_intens / p.gamma)
                            * (-p.k * d).exp()
                            * (1.0 - (-p.gamma * (d + d_th_a)).exp())
                    },
                    -0.5,
                    4.0,
                );
                opt_ask[i] = da;
                fa_val = val;
            }

            next_theta[i] = theta[i] + actual_dt * (drift + fb_val + fa_val);
        }
        theta.copy_from_slice(&next_theta);
    }
    (theta, opt_bid, opt_ask)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glt::{eigendecompose_m, solve_v_jacobi, theta_from_v};

    #[test]
    fn test_nonlinear_hjb_matches_v_system() {
        let p = GltParams::default_params();
        let eig = eigendecompose_m(&p);
        let v_exact = solve_v_jacobi(&p, &eig, p.t_horizon);
        let th_exact = theta_from_v(&p, &v_exact);
        let (th_fd, _, _) = solve_nonlinear_hjb_theta(&p, 5e-4, 0.0);
        let max_err = th_exact
            .iter()
            .zip(th_fd.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(max_err < 1e-2, "nonlinear HJB FD max_err={max_err}");
    }
}
