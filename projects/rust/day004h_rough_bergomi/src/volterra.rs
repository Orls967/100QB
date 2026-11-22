//! Fractional Volterra / Riemann-Liouville process `W_tilde^H(t) = sqrt(2H) int_0^t (t-s)^{H-1/2} dB_s`:
//! exact covariance (closed-form 2F1 + singularity-free Gauss-Legendre),
//! 2n x 2n Cholesky exact reference simulator, and Bennedsen-Lunde-Pakkanen hybrid scheme (kappa=1).

#![allow(clippy::needless_range_loop)]

use crate::linalg;
use crate::parallel;
use crate::quad::GaussLegendre;
use crate::rng::{derive_seed, Rng};
use crate::stats;

/// Exact variance `Var(W_tilde^H_t) = t^{2H}`.
#[inline]
pub fn exact_var(t: f64, h: f64) -> f64 {
    if t <= 0.0 {
        0.0
    } else {
        t.powf(2.0 * h)
    }
}

/// Evaluates Gauss hypergeometric function `2F1(a, b; c; z)` for `|z| < 1` via power series.
fn hyp2f1_series(a: f64, b: f64, c: f64, z: f64) -> f64 {
    let mut sum = 1.0_f64;
    let mut term = 1.0_f64;
    for n in 0..400 {
        let nf = n as f64;
        term *= (a + nf) * (b + nf) / ((c + nf) * (nf + 1.0)) * z;
        sum += term;
        if term.abs() < 1e-15 * sum.abs() {
            break;
        }
    }
    sum
}

/// Lanczos approximation for `ln(Gamma(x))` for `x > 0`.
fn lgamma(x: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        let s = (std::f64::consts::PI * x).sin();
        return std::f64::consts::PI.ln() - s.abs().ln() - lgamma(1.0 - x);
    }
    let z = x - 1.0;
    let mut a = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        a += c / (z + i as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}

fn gamma_fn(x: f64) -> f64 {
    if x < 0.0 {
        // Gamma(x) = pi / (sin(pi*x) * Gamma(1-x))
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * lgamma(1.0 - x).exp())
    } else {
        lgamma(x).exp()
    }
}

/// Closed-form covariance `Cov(W_tilde^H_s, W_tilde^H_t)` using `2F1(1/2 - H, 1; H + 3/2; s/t)`
/// with DLMF 15.8.4 `1-z` analytic continuation near `s/t -> 1`.
pub fn exact_cov_closed_form(s_in: f64, t_in: f64, h: f64) -> f64 {
    let (s, t) = if s_in <= t_in {
        (s_in, t_in)
    } else {
        (t_in, s_in)
    };
    if s <= 0.0 || t <= 0.0 {
        return 0.0;
    }
    if (s - t).abs() <= 1e-15 * t {
        return exact_var(t, h);
    }
    if (h - 0.5).abs() < 1e-14 {
        return s;
    }
    let z = s / t;
    let a = 0.5 - h;
    let b = 1.0;
    let c = h + 1.5;
    let pre = (2.0 * h / (h + 0.5)) * s.powf(h + 0.5) * t.powf(h - 0.5);
    let f_val = if z <= 0.5 {
        hyp2f1_series(a, b, c, z)
    } else {
        // DLMF 15.8.4 with b = 1:
        // Term 1: (Gamma(c) Gamma(c-a-1) / (Gamma(c-a) Gamma(c-1))) * 2F1(a, 1; 2+a-c; 1-z)
        //       = ((h+0.5)/(2h)) * 2F1(0.5-h, 1; 1-2h; 1-z)
        // Term 2: (Gamma(c) Gamma(1+a-c) / Gamma(a)) * (1-z)^{2h} * z^{-(h+0.5)}
        let omz = 1.0 - z;
        let t1 = ((h + 0.5) / (2.0 * h)) * hyp2f1_series(a, 1.0, 1.0 - 2.0 * h, omz);
        let coef2 = gamma_fn(c) * gamma_fn(-2.0 * h) / gamma_fn(a);
        let t2 = coef2 * omz.powf(2.0 * h) * z.powf(-(h + 0.5));
        t1 + t2
    };
    pre * f_val
}

/// Singularity-free Gauss-Legendre quadrature for `Cov(W_tilde^H_s, W_tilde^H_t)`
/// via substitution `s - u = s * y^{1/(H + 1/2)}` which cancels `(s-u)^{H-1/2}`.
pub fn exact_cov_quadrature(s_in: f64, t_in: f64, h: f64, gl: &GaussLegendre) -> f64 {
    let (s, t) = if s_in <= t_in {
        (s_in, t_in)
    } else {
        (t_in, s_in)
    };
    if s <= 0.0 || t <= 0.0 {
        return 0.0;
    }
    if (s - t).abs() <= 1e-15 * t {
        return exact_var(t, h);
    }
    let hp = h + 0.5;
    let inv_hp = 1.0 / hp;
    let alpha = h - 0.5;
    let pre = (2.0 * h * inv_hp) * s.powf(hp);
    pre * gl.integrate(0.0, 1.0, |y| {
        let v = s * y.powf(inv_hp);
        (t - s + v).powf(alpha)
    })
}

/// Exact cross-covariance `Cov(W_tilde^H_t, B_{t_j} - B_{t_{j-1}})`.
#[inline]
pub fn exact_cross_cov(t: f64, t_prev: f64, t_curr: f64, h: f64) -> f64 {
    let hp = h + 0.5;
    let a = (t - t_prev).max(0.0).powf(hp);
    let b = (t - t_curr).max(0.0).powf(hp);
    ((2.0 * h).sqrt() / hp) * (a - b)
}

/// Builds the `2n x 2n` exact covariance matrix for `(W_tilde^H_{t_1..t_n}, dB_1..dB_n)`
/// where `t_i = i * dt` and `dt = t_horizon / n`.
pub fn build_joint_cov_matrix(n: usize, t_horizon: f64, h: f64) -> Vec<Vec<f64>> {
    let dt = t_horizon / n as f64;
    let dim = 2 * n;
    let mut cov = vec![vec![0.0_f64; dim]; dim];
    for i in 0..n {
        let ti = (i + 1) as f64 * dt;
        for j in i..n {
            let tj = (j + 1) as f64 * dt;
            let c = exact_cov_closed_form(ti, tj, h);
            cov[i][j] = c;
            cov[j][i] = c;
        }
        for j in 0..n {
            let t_prev = j as f64 * dt;
            let t_curr = (j + 1) as f64 * dt;
            let cc = exact_cross_cov(ti, t_prev, t_curr, h);
            cov[i][n + j] = cc;
            cov[n + j][i] = cc;
        }
        cov[n + i][n + i] = dt;
    }
    cov
}

#[derive(Clone, Debug)]
pub struct CholeskyRefSummary {
    pub n: usize,
    pub n_paths: usize,
    pub max_chol_recon_err: f64,
    pub max_cov_z_score: f64,
    pub var_w_terminal_emp: f64,
    pub var_w_terminal_exact: f64,
}

/// Simulates `n_paths` joint `(W_tilde^H, dB)` via `2n x 2n` Cholesky and checks empirical vs exact covariance.
pub fn run_cholesky_reference_check(
    n: usize,
    t_horizon: f64,
    h: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    salt: u64,
) -> CholeskyRefSummary {
    let cov = build_joint_cov_matrix(n, t_horizon, h);
    let chol = linalg::cholesky(&cov).expect("joint covariance must be SPD");
    let dim = 2 * n;

    // Accumulate outer products for empirical second moments E[X_a X_b] (mean is 0)
    // Only check W-block and selected cross-entries + terminal variance for speed and clarity.
    let m_len = dim * dim;
    let sum_xx = parallel::run_chunk_reduce(
        n_paths,
        n_workers,
        |start, end| {
            let mut local = vec![0.0_f64; m_len];
            let mut z = vec![0.0_f64; dim];
            let mut x = vec![0.0_f64; dim];
            for r in start..end {
                let mut rg = Rng::new(derive_seed(seed, salt, r as u64));
                for zj in z.iter_mut() {
                    *zj = rg.normal();
                }
                for i in 0..dim {
                    let row = &chol[i];
                    let mut s = 0.0_f64;
                    for j in 0..=i {
                        s += row[j] * z[j];
                    }
                    x[i] = s;
                }
                for i in 0..dim {
                    let xi = x[i];
                    let base = i * dim;
                    for j in i..dim {
                        local[base + j] += xi * x[j];
                    }
                }
            }
            local
        },
        |parts| {
            let mut tot = vec![0.0_f64; m_len];
            for p in parts {
                for (a, b) in tot.iter_mut().zip(p.iter()) {
                    *a += *b;
                }
            }
            tot
        },
    );

    // Cholesky reconstruction error max |L L^T - Cov|
    let mut max_chol_err = 0.0_f64;
    for i in 0..dim {
        for j in i..dim {
            let mut s = 0.0_f64;
            for k in 0..=i {
                s += chol[i][k] * chol[j][k];
            }
            max_chol_err = max_chol_err.max((s - cov[i][j]).abs());
        }
    }

    let nf = n_paths as f64;
    let mut max_z = 0.0_f64;
    for i in 0..dim {
        for j in i..dim {
            let emp = sum_xx[i * dim + j] / nf;
            let ex = cov[i][j];
            let var_prod = cov[i][i] * cov[j][j] + ex * ex;
            let se = (var_prod / nf).sqrt().max(1e-12);
            let z = (emp - ex).abs() / se;
            if z > max_z {
                max_z = z;
            }
        }
    }

    let var_w_emp = sum_xx[(n - 1) * dim + (n - 1)] / nf;
    let var_w_ex = cov[n - 1][n - 1];

    CholeskyRefSummary {
        n,
        n_paths,
        max_chol_recon_err: max_chol_err,
        max_cov_z_score: max_z,
        var_w_terminal_emp: var_w_emp,
        var_w_terminal_exact: var_w_ex,
    }
}

/// Precomputed Bennedsen-Lunde-Pakkanen hybrid scheme (`kappa = 1`) weights on grid `n`, `dt = T/n`.
/// Cell 1 (`[t_{i-1}, t_i]`): exact stochastic integral correlated with `dB_i` (`rho_h`).
/// Cells `k = 2..=n`: right-endpoint kernel weight `sqrt(2H) * ((k-1)*dt)^{H - 1/2}`.
#[derive(Clone, Debug)]
pub struct HybridWeights {
    pub n: usize,
    pub dt: f64,
    pub sqrt_dt: f64,
    pub h: f64,
    pub dt_pow_h: f64,
    pub rho_h: f64,
    pub rho_perp: f64,
    pub conv_weights: Vec<f64>,
}

impl HybridWeights {
    pub fn new(n: usize, t_horizon: f64, h: f64) -> Self {
        let dt = t_horizon / n as f64;
        let sqrt_dt = dt.sqrt();
        let dt_pow_h = dt.powf(h);
        let rho_h = ((2.0 * h).sqrt() / (h + 0.5)).clamp(-1.0, 1.0);
        let rho_perp = (1.0 - rho_h * rho_h).max(0.0).sqrt();
        let alpha = h - 0.5;
        let pre = (2.0 * h).sqrt();
        let mut conv_weights = vec![0.0_f64; n + 1];
        for (k, w) in conv_weights.iter_mut().enumerate().skip(2) {
            *w = pre * (((k - 1) as f64) * dt).powf(alpha);
        }
        Self {
            n,
            dt,
            sqrt_dt,
            h,
            dt_pow_h,
            rho_h,
            rho_perp,
            conv_weights,
        }
    }

    /// Generates `w_tilde[0..n]` (at `t_1..t_n`) and `db[0..n]` from `2n` independent standard normals
    /// `(u_0..u_{n-1}, v_0..v_{n-1})`.
    pub fn generate_path(&self, u: &[f64], v: &[f64], w_tilde: &mut [f64], db: &mut [f64]) {
        let n = self.n;
        for i in 0..n {
            let db_i = self.sqrt_dt * u[i];
            db[i] = db_i;
            let z_i1 = self.dt_pow_h * (self.rho_h * u[i] + self.rho_perp * v[i]);
            let mut acc = z_i1;
            // sum_{k=2}^{i+1} conv_weights[k] * db[i + 1 - k]
            // Let j = i + 1 - k (from 0..i), so k = i + 1 - j >= 2
            for j in 0..i {
                acc += self.conv_weights[i + 1 - j] * db[j];
            }
            w_tilde[i] = acc;
        }
    }
}

#[derive(Clone, Debug)]
pub struct HybridBiasRow {
    pub h: f64,
    pub n: usize,
    pub var_w_exact: f64,
    pub var_w_hybrid_theoretical: f64,
    pub var_w_hybrid_mc: f64,
    pub var_bias_theoretical: f64,
    pub max_cov_abs_err_theoretical: f64,
}

/// Computes exact theoretical variance & covariance of the `kappa=1` right-endpoint hybrid scheme
/// and compares against MC across `n in {50, 100, 250, 500}`.
pub fn hybrid_bias_study(
    h: f64,
    t_horizon: f64,
    ns: &[usize],
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<HybridBiasRow> {
    let mut rows = Vec::with_capacity(ns.len());
    let var_ex = exact_var(t_horizon, h);
    for &n in ns {
        let hw = HybridWeights::new(n, t_horizon, h);
        // Theoretical hybrid variance at T = t_n:
        // Var(W_n) = dt^{2H} + sum_{k=2}^n conv_weights[k]^2 * dt
        let mut var_hyb_th = hw.dt_pow_h * hw.dt_pow_h;
        for k in 2..=n {
            var_hyb_th += hw.conv_weights[k] * hw.conv_weights[k] * hw.dt;
        }
        // Max covariance error against exact on 10 evenly spaced indices
        let stride = (n / 10).max(1);
        let mut max_cov_err = 0.0_f64;
        for i in (stride - 1..n).step_by(stride) {
            let ti = (i + 1) as f64 * hw.dt;
            for j in (i..n).step_by(stride) {
                let tj = (j + 1) as f64 * hw.dt;
                let ex_c = exact_cov_closed_form(ti, tj, h);
                // Theoretical hybrid Cov(W_i, W_j) for i <= j:
                // For i == j: dt^{2H} + sum_{m=0}^{i-1} w[i+1-m]^2 * dt
                // For i < j: cov(Z_{i,1}, dB_i)*w[j+1-i] + sum_{m=0}^{i-1} w[i+1-m]*w[j+1-m]*dt
                let hyb_c = if i == j {
                    let mut s = hw.dt_pow_h * hw.dt_pow_h;
                    for m in 0..i {
                        s += hw.conv_weights[i + 1 - m] * hw.conv_weights[i + 1 - m] * hw.dt;
                    }
                    s
                } else {
                    let mut s = hw.dt_pow_h * hw.rho_h * hw.sqrt_dt * hw.conv_weights[j + 1 - i];
                    for m in 0..i {
                        s += hw.conv_weights[i + 1 - m] * hw.conv_weights[j + 1 - m] * hw.dt;
                    }
                    s
                };
                max_cov_err = max_cov_err.max((hyb_c - ex_c).abs());
            }
        }

        let w_end_sq = parallel::run_chunk_reduce(
            n_paths,
            n_workers,
            |start, end| {
                let mut sum_sq = 0.0_f64;
                let mut u = vec![0.0_f64; n];
                let mut v = vec![0.0_f64; n];
                let mut wt = vec![0.0_f64; n];
                let mut db = vec![0.0_f64; n];
                for r in start..end {
                    let mut rg = Rng::new(derive_seed(seed, 0x1300 ^ (n as u64), r as u64));
                    for idx in 0..n {
                        u[idx] = rg.normal();
                        v[idx] = rg.normal();
                    }
                    hw.generate_path(&u, &v, &mut wt, &mut db);
                    let w_last = wt[n - 1];
                    sum_sq += w_last * w_last;
                }
                sum_sq
            },
            |parts| parts.into_iter().sum::<f64>(),
        );
        let var_mc = w_end_sq / n_paths as f64;
        rows.push(HybridBiasRow {
            h,
            n,
            var_w_exact: var_ex,
            var_w_hybrid_theoretical: var_hyb_th,
            var_w_hybrid_mc: var_mc,
            var_bias_theoretical: var_hyb_th - var_ex,
            max_cov_abs_err_theoretical: max_cov_err,
        });
    }
    let _ = stats::mean(&[]);
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_cov_closed_form_matches_quadrature() {
        let gl = GaussLegendre::new(64);
        let h = 0.1;
        let times = [0.02, 0.1, 0.25, 0.5, 0.77, 1.0];
        for &s in &times {
            for &t in &times {
                let c_cf = exact_cov_closed_form(s, t, h);
                let c_q = exact_cov_quadrature(s, t, h, &gl);
                let rel = (c_cf - c_q).abs() / c_q.abs().max(1e-12);
                assert!(rel < 1e-10, "s={s} t={t} cf={c_cf} q={c_q} rel={rel}");
            }
        }
    }

    #[test]
    fn test_brownian_limit_h_half() {
        let h = 0.5;
        let gl = GaussLegendre::new(32);
        for &s in &[0.1_f64, 0.4, 0.9] {
            for &t in &[0.1_f64, 0.4, 0.9] {
                let expected = s.min(t);
                assert!((exact_cov_closed_form(s, t, h) - expected).abs() < 1e-12);
                assert!((exact_cov_quadrature(s, t, h, &gl) - expected).abs() < 1e-12);
            }
        }
        let hw = HybridWeights::new(20, 1.0, 0.5);
        let mut u = vec![0.1; 20];
        u[3] = -0.7;
        let v = vec![0.9; 20];
        let mut wt = vec![0.0; 20];
        let mut db = vec![0.0; 20];
        hw.generate_path(&u, &v, &mut wt, &mut db);
        let mut cum = 0.0;
        for i in 0..20 {
            cum += db[i];
            assert!((wt[i] - cum).abs() < 1e-13);
        }
    }
}
