use crate::parallel::reduce_tuple4_fixed_blocks;
use crate::rng::{derive_seed, Rng};

#[inline]
pub fn reserves(l: f64, p: f64) -> (f64, f64) {
    let sp = p.sqrt();
    (l / sp, l * sp)
}

#[inline]
pub fn pool_value(l: f64, p: f64) -> f64 {
    2.0 * l * p.sqrt()
}

#[inline]
pub fn pool_delta(l: f64, p: f64) -> f64 {
    l / p.sqrt()
}

#[inline]
pub fn pool_gamma(l: f64, p: f64) -> f64 {
    -0.5 * l * p.powf(-1.5)
}

#[inline]
pub fn impermanent_loss(r: f64) -> f64 {
    2.0 * r.sqrt() / (1.0 + r) - 1.0
}

#[inline]
pub fn expected_pool_value(v0: f64, sigma: f64, t: f64) -> f64 {
    v0 * (-sigma * sigma * t / 8.0).exp()
}

#[inline]
pub fn variance_pool_value(v0: f64, sigma: f64, t: f64) -> f64 {
    v0 * v0 * (1.0 - (-sigma * sigma * t / 4.0).exp())
}

#[inline]
pub fn lvr_rate_analytic(v: f64, sigma: f64) -> f64 {
    sigma * sigma * v / 8.0
}

/// Simulates `n_paths` GBM paths under `mu = 0` and computes MC `E[V_T]`, `SE(E[V_T])`, `Var[V_T]`.
pub fn mc_pool_moments(
    v0: f64,
    sigma: f64,
    t: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64) {
    let drift = -0.5 * sigma * sigma * t;
    let vol_dt = sigma * t.sqrt();
    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x1001, r as u64));
        let pt = (drift + vol_dt * rg.normal()).exp();
        let vt = v0 * pt.sqrt();
        [vt, vt * vt, 0.0, 0.0]
    });
    let nf = n_paths as f64;
    let mean_v = sums[0] / nf;
    let var_v = (sums[1] / nf - mean_v * mean_v).max(0.0) * (nf / (nf - 1.0).max(1.0));
    let se_v = (var_v / nf).sqrt();
    (mean_v, se_v, var_v)
}

/// Verifies pathwise Ito LVR: discrete delta-rebalancing portfolio `V_rebal - V_pool`
/// vs integrated `(sigma^2 / 8) * int_0^T V_t dt`.
pub fn pathwise_lvr_ito(
    v0: f64,
    sigma: f64,
    t_mat: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64) {
    let l = 0.5 * v0;
    let dt = t_mat / n_steps as f64;
    let drift = -0.5 * sigma * sigma * dt;
    let vol_dt = sigma * dt.sqrt();
    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x1002, r as u64));
        let mut p = 1.0_f64;
        let mut v_rebal = v0;
        let mut lvr_int = 0.0_f64;
        for _ in 0..n_steps {
            let v_t = pool_value(l, p);
            let d_t = pool_delta(l, p);
            lvr_int += lvr_rate_analytic(v_t, sigma) * dt;
            let p_next = p * (drift + vol_dt * rg.normal()).exp();
            v_rebal += d_t * (p_next - p);
            p = p_next;
        }
        let lvr_rep = v_rebal - pool_value(l, p);
        let diff = lvr_rep - lvr_int;
        [lvr_rep, lvr_int, diff * diff, 0.0]
    });
    let nf = n_paths as f64;
    let mean_rep = sums[0] / nf;
    let mean_int = sums[1] / nf;
    let rmse_path = (sums[2] / nf).sqrt();
    (mean_rep, mean_int, rmse_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpmm_identities_and_lvr_ito() {
        let l = 0.5_f64;
        for &r in &[0.25_f64, 0.5, 1.0, 2.0, 4.0] {
            let (x, y) = reserves(l, r);
            assert!((x * y - l * l).abs() < 1e-14);
            let v_pool = pool_value(l, r);
            let v_hodl = 0.5 * r + 0.5;
            assert!((v_pool / v_hodl - 1.0 - impermanent_loss(r)).abs() < 1e-14);
        }
        let (rep, int_a, rmse) = pathwise_lvr_ito(1.0, 0.8, 1.0 / 365.0, 200, 10_000, 4, 20261122);
        assert!((rep - int_a).abs() < 1e-4, "rep={rep} int={int_a}");
        assert!(rmse < 1e-3, "rmse={rmse}");
    }
}
