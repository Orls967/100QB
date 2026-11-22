//! Guéant–Lehalle–Tapia (GLT) exact v-system solvers (Jacobi eigen, scaling-and-squaring, RK4+Richardson),
//! optimal quotes, stationary ground-state limit, closed-form approximation, and Q sensitivity.

use crate::linalg;
use crate::ode;

#[derive(Clone, Copy, Debug)]
pub struct GltParams {
    pub t_horizon: f64,
    pub sigma: f64,
    pub a_intens: f64,
    pub k: f64,
    pub gamma: f64,
    pub q_max: isize,
}

impl GltParams {
    pub fn default_params() -> Self {
        Self {
            t_horizon: 1.0,
            sigma: 0.3,
            a_intens: 140.0,
            k: 1.5,
            gamma: 0.1,
            q_max: 20,
        }
    }

    #[inline]
    pub fn alpha(&self) -> f64 {
        0.5 * self.k * self.gamma * self.sigma * self.sigma
    }

    #[inline]
    pub fn eta(&self) -> f64 {
        self.a_intens * (1.0 + self.gamma / self.k).powf(-(1.0 + self.k / self.gamma))
    }

    #[inline]
    pub fn half_spread_zero_inv(&self) -> f64 {
        (1.0 / self.gamma) * (1.0 + self.gamma / self.k).ln()
    }

    #[inline]
    pub fn dim(&self) -> usize {
        (2 * self.q_max + 1) as usize
    }

    #[inline]
    pub fn idx(&self, q: isize) -> usize {
        (q + self.q_max) as usize
    }
}

pub fn build_m_matrix(p: &GltParams) -> Vec<Vec<f64>> {
    let n = p.dim();
    let alpha = p.alpha();
    let eta = p.eta();
    let mut m = linalg::zeros(n, n);
    for i in 0..n {
        let q = i as isize - p.q_max;
        m[i][i] = alpha * (q * q) as f64;
        if i > 0 {
            m[i][i - 1] = -eta;
        }
        if i + 1 < n {
            m[i][i + 1] = -eta;
        }
    }
    m
}

pub struct EigenDecomposition {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<Vec<f64>>,
    pub proj_ones: Vec<f64>,
}

pub fn eigendecompose_m(p: &GltParams) -> EigenDecomposition {
    let m = build_m_matrix(p);
    let n = m.len();
    let (vals, mut vecs, _) = linalg::symmetric_eigen_jacobi_asc(&m, 120, 1e-14).unwrap();
    // Orient each eigenvector so sum(v) >= 0 (ground state is strictly positive)
    for vec in &mut vecs {
        let s: f64 = vec.iter().sum();
        if s < 0.0 || (s.abs() < 1e-12 && vec[n / 2] < 0.0) {
            for x in vec.iter_mut() {
                *x = -*x;
            }
        }
    }
    let proj_ones: Vec<f64> = vecs.iter().map(|v| v.iter().sum()).collect();
    EigenDecomposition {
        eigenvalues: vals,
        eigenvectors: vecs,
        proj_ones,
    }
}

pub fn solve_v_jacobi(p: &GltParams, eig: &EigenDecomposition, tau: f64) -> Vec<f64> {
    let n = p.dim();
    let mut v = vec![0.0; n];
    for (m_idx, &lam) in eig.eigenvalues.iter().enumerate() {
        let coef = (-lam * tau).exp() * eig.proj_ones[m_idx];
        let vec_m = &eig.eigenvectors[m_idx];
        for i in 0..n {
            v[i] += coef * vec_m[i];
        }
    }
    v
}

pub fn solve_v_mat_exp(p: &GltParams, tau: f64) -> Vec<f64> {
    let m = build_m_matrix(p);
    let n = m.len();
    let mut neg_tau_m = linalg::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            neg_tau_m[i][j] = -tau * m[i][j];
        }
    }
    let exp_m = linalg::mat_exp_scaling_squaring(&neg_tau_m);
    let ones = vec![1.0; n];
    linalg::mat_vec(&exp_m, &ones)
}

pub fn solve_v_rk4_richardson(p: &GltParams, tau: f64, n_coarse: usize) -> Vec<f64> {
    let n = p.dim();
    let alpha = p.alpha();
    let eta = p.eta();
    let y0 = vec![1.0; n];
    let deriv = |_t: f64, y: &[f64], dydt: &mut [f64]| {
        for i in 0..n {
            let q = i as isize - p.q_max;
            let mut val = -alpha * (q * q) as f64 * y[i];
            if i > 0 {
                val += eta * y[i - 1];
            }
            if i + 1 < n {
                val += eta * y[i + 1];
            }
            dydt[i] = val;
        }
    };
    ode::rk4_richardson(0.0, tau, &y0, n_coarse, deriv)
}

pub fn theta_from_v(p: &GltParams, v: &[f64]) -> Vec<f64> {
    v.iter().map(|&vi| vi.ln() / p.k).collect()
}

/// Returns `(delta_bid, delta_ask)` for each inventory index `0..dim`.
/// At `q = Q` (`idx = dim-1`), `delta_bid = +inf` (quote turned off);
/// at `q = -Q` (`idx = 0`), `delta_ask = +inf` (quote turned off).
pub fn quotes_from_v(p: &GltParams, v: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = p.dim();
    let base = p.half_spread_zero_inv();
    let inv_k = 1.0 / p.k;
    let mut db = vec![f64::INFINITY; n];
    let mut da = vec![f64::INFINITY; n];
    for i in 0..n {
        if i + 1 < n {
            db[i] = inv_k * (v[i] / v[i + 1]).ln() + base;
        }
        if i > 0 {
            da[i] = inv_k * (v[i] / v[i - 1]).ln() + base;
        }
    }
    (db, da)
}

/// Closed-form asymptotic approximation (small gamma, large Q) for bid quote at inventory `q`.
pub fn approx_bid_quote(p: &GltParams, q: isize) -> f64 {
    let base = p.half_spread_zero_inv();
    let slope = (p.sigma * p.sigma * p.gamma / (2.0 * p.k * p.eta())).sqrt();
    base + 0.5 * (2 * q + 1) as f64 * slope
}

pub fn approx_ask_quote(p: &GltParams, q: isize) -> f64 {
    approx_bid_quote(p, -q)
}

#[derive(Clone, Debug)]
pub struct QuoteTable {
    pub params: GltParams,
    pub n_time: usize,
    pub dt: f64,
    pub bid: Vec<Vec<f64>>,
    pub ask: Vec<Vec<f64>>,
}

impl QuoteTable {
    pub fn build(p: GltParams, n_time: usize) -> Self {
        let eig = eigendecompose_m(&p);
        let dt = p.t_horizon / n_time as f64;
        let mut bid = Vec::with_capacity(n_time + 1);
        let mut ask = Vec::with_capacity(n_time + 1);
        for step in 0..=n_time {
            let t = step as f64 * dt;
            let tau = (p.t_horizon - t).max(0.0);
            let v = solve_v_jacobi(&p, &eig, tau);
            let (db, da) = quotes_from_v(&p, &v);
            bid.push(db);
            ask.push(da);
        }
        Self {
            params: p,
            n_time,
            dt,
            bid,
            ask,
        }
    }

    #[inline]
    pub fn interp_quotes(&self, t: f64, q: isize) -> (f64, f64) {
        let idx_q = self.params.idx(q);
        let pos = (t / self.dt).clamp(0.0, self.n_time as f64);
        let k0 = (pos.floor() as usize).min(self.n_time);
        let k1 = (k0 + 1).min(self.n_time);
        let w = pos - k0 as f64;
        let db = if q >= self.params.q_max {
            f64::INFINITY
        } else {
            self.bid[k0][idx_q] * (1.0 - w) + self.bid[k1][idx_q] * w
        };
        let da = if q <= -self.params.q_max {
            f64::INFINITY
        } else {
            self.ask[k0][idx_q] * (1.0 - w) + self.ask[k1][idx_q] * w
        };
        (db, da)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_three_solvers_agree_and_quote_symmetry() {
        let p = GltParams::default_params();
        let eig = eigendecompose_m(&p);
        for &tau in &[0.1, 0.5, 1.0] {
            let vj = solve_v_jacobi(&p, &eig, tau);
            let ve = solve_v_mat_exp(&p, tau);
            let vr = solve_v_rk4_richardson(&p, tau, 6000);
            for i in 0..p.dim() {
                let rel_je = (vj[i] - ve[i]).abs() / vj[i].abs().max(1e-15);
                let rel_jr = (vj[i] - vr[i]).abs() / vj[i].abs().max(1e-15);
                assert!(
                    rel_je < 1e-10,
                    "Jacobi vs MatExp rel err {rel_je} at tau={tau} i={i}"
                );
                assert!(
                    rel_jr < 1e-10,
                    "Jacobi vs RK4 rel err {rel_jr} at tau={tau} i={i}"
                );
            }
            let (db, da) = quotes_from_v(&p, &vj);
            for q in -p.q_max..p.q_max {
                let i_pos = p.idx(q);
                let i_neg = p.idx(-q);
                assert!((db[i_pos] - da[i_neg]).abs() < 1e-12);
                if q + 1 < p.q_max {
                    assert!(db[p.idx(q + 1)] > db[i_pos]);
                }
            }
        }
    }
}
