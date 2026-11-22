use crate::hedge::PathData;
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};

#[derive(Clone, Copy, Debug)]
pub struct HestonParams {
    pub s0: f64,
    pub k: f64,
    pub r: f64,
    pub kappa: f64,
    pub theta: f64,
    pub xi: f64,
    pub rho: f64,
    pub v0: f64,
    pub t_mat: f64,
    pub n_steps: usize,
}

/// Simulates Heston paths using Andersen (2008) Quadratic Exponential (QE) scheme
/// with exact martingale-preserving drift correction `K_0^*`.
pub fn generate_heston_qe_paths(
    p: &HestonParams,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    tag: u64,
) -> Vec<PathData> {
    let dt = p.t_mat / p.n_steps as f64;
    let exp_kdt = (-p.kappa * dt).exp();
    let one_m_exp = 1.0 - exp_kdt;
    let psi_c = 1.5_f64;

    let gamma1 = 0.5_f64;
    let gamma2 = 0.5_f64;
    let k1 = gamma1 * dt * (p.kappa * p.rho / p.xi - 0.5) - p.rho / p.xi;
    let k2 = gamma2 * dt * (p.kappa * p.rho / p.xi - 0.5) + p.rho / p.xi;
    let k3 = gamma1 * dt * (1.0 - p.rho * p.rho);
    let k4 = gamma2 * dt * (1.0 - p.rho * p.rho);
    let a_mart = k2 + 0.5 * k4;

    let c1_s2 = p.xi * p.xi * exp_kdt * one_m_exp / p.kappa;
    let c2_s2 = p.theta * p.xi * p.xi * one_m_exp * one_m_exp / (2.0 * p.kappa);

    run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, tag, r as u64));
        let mut s_vec = Vec::with_capacity(p.n_steps + 1);
        let mut v_vec = Vec::with_capacity(p.n_steps + 1);
        let mut log_s = p.s0.ln();
        let mut v = p.v0;
        s_vec.push(p.s0);
        v_vec.push(v);

        for _ in 0..p.n_steps {
            let m = p.theta + (v - p.theta) * exp_kdt;
            let s2 = v * c1_s2 + c2_s2;
            let psi = s2 / (m * m).max(1e-30);

            let (v_next, k0_star) = if psi <= psi_c {
                let inv_psi = 2.0 / psi;
                let b2 = inv_psi - 1.0 + (inv_psi * (inv_psi - 1.0)).sqrt().max(0.0);
                let b = b2.sqrt();
                let a = m / (1.0 + b2);
                let zv = rg.normal();
                let vn = a * (b + zv) * (b + zv);
                let denom = (1.0 - 2.0 * a_mart * a).max(1e-12);
                let k0 = -a_mart * b2 * a / denom + 0.5 * denom.ln() - (k1 + 0.5 * k3) * v;
                (vn, k0)
            } else {
                let p_exp = (psi - 1.0) / (psi + 1.0);
                let beta = (1.0 - p_exp) / m.max(1e-15);
                let u = rg.uniform();
                let vn = if u <= p_exp {
                    0.0
                } else {
                    ((1.0 - p_exp) / (1.0 - u)).ln() / beta
                };
                let k0 = -(p_exp + beta * (1.0 - p_exp) / (beta - a_mart).max(1e-12)).ln()
                    - (k1 + 0.5 * k3) * v;
                (vn, k0)
            };

            let zs = rg.normal();
            let var_s = (k3 * v + k4 * v_next).max(0.0);
            log_s += p.r * dt + k0_star + k1 * v + k2 * v_next + var_s.sqrt() * zs;
            v = v_next;
            s_vec.push(log_s.exp());
            v_vec.push(v);
        }

        PathData {
            s: s_vec,
            v: Some(v_vec),
        }
    })
}
