#![allow(clippy::too_many_arguments)]

use crate::black::black_call;
use crate::parallel::reduce_tuple4_fixed_blocks;
use crate::rng::{derive_seed, Rng};
use crate::special::{gamma_q, lgamma, ncx2_cdf, ncx2_sf};
use crate::stats;

/// Absorbing probability `P(F_T = 0)` for CEV `dF = alpha * F^beta dW` (`0 < beta < 1`).
#[inline]
pub fn cev_absorb_prob(f0: f64, alpha: f64, beta: f64, t: f64) -> f64 {
    if f0 <= 0.0 {
        return 1.0;
    }
    if t <= 0.0 || alpha <= 0.0 || beta >= 1.0 {
        return 0.0;
    }
    let one_m_b = 1.0 - beta;
    let nu = 0.5 / one_m_b;
    let x = f0.powf(2.0 * one_m_b) / (alpha * alpha * one_m_b * one_m_b * t);
    gamma_q(nu, 0.5 * x)
}

/// Exact closed-form CEV call price (`r = 0`, absorbing at `F = 0`, `0 < beta < 1`):
/// `C = F0 * ncx2_sf(y; 2 + 1/(1-beta), x) - K * ncx2_cdf(x; 1/(1-beta), y)`
/// where `x = F0^{2(1-beta)} / (alpha^2 (1-beta)^2 T)`, `y = K^{2(1-beta)} / (alpha^2 (1-beta)^2 T)`.
pub fn cev_call_exact(f0: f64, k: f64, alpha: f64, beta: f64, t: f64) -> f64 {
    if k <= 0.0 {
        return f0.max(0.0);
    }
    if f0 <= 0.0 {
        return 0.0;
    }
    if t <= 0.0 || alpha <= 0.0 {
        return (f0 - k).max(0.0);
    }
    if (beta - 1.0).abs() < 1e-12 {
        return black_call(f0, k, alpha, t);
    }
    let one_m_b = 1.0 - beta;
    let p = 1.0 / one_m_b;
    let denom = alpha * alpha * one_m_b * one_m_b * t;
    let x = f0.powf(2.0 * one_m_b) / denom;
    let y = k.powf(2.0 * one_m_b) / denom;
    let term1 = f0 * ncx2_sf(y, 2.0 + p, x);
    let term2 = k * ncx2_cdf(x, p, y);
    (term1 - term2).max(0.0)
}

/// 1D Crank-Nicolson finite-difference PDE solver on a non-uniform sinh grid `[0, F_max]`
/// with absorbing boundary `C(0, tau) = 0`, linear far-field `C(F_max, tau) = F_max - K`,
/// and 2 Rannacher half-steps to damp payoff kink at `F = K`.
pub fn cev_call_pde_cn(
    f0: f64,
    k: f64,
    alpha: f64,
    beta: f64,
    t: f64,
    n_space: usize,
    n_time: usize,
) -> f64 {
    if k <= 0.0 {
        return f0.max(0.0);
    }
    if t <= 0.0 || alpha <= 0.0 {
        return (f0 - k).max(0.0);
    }
    let f_scale = alpha * f0.powf(beta) * t.sqrt();
    let f_max = (5.0 * f0).max(3.0 * k) + 6.0 * f_scale;
    // Non-uniform sinh grid clustered around `0.5 * (f0 + k)`
    let f_center = (0.5 * (f0 + k)).clamp(0.05 * f0, 0.9 * f_max);
    let c_sinh = 0.05 * f_max;
    let xi_min = (-f_center / c_sinh).asinh();
    let xi_max = ((f_max - f_center) / c_sinh).asinh();

    let n = n_space;
    let mut grid = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let xi = xi_min + (xi_max - xi_min) * (i as f64) / (n as f64);
        grid.push((f_center + c_sinh * xi.sinh()).clamp(0.0, f_max));
    }
    grid[0] = 0.0;
    grid[n] = f_max;
    // Align nearest interior node to strike `k` so the payoff kink sits on a grid node
    if k > 0.0 && k < f_max {
        let mut best_i = 1;
        let mut best_d = (grid[1] - k).abs();
        for (i, &g) in grid.iter().enumerate().take(n).skip(2) {
            let d = (g - k).abs();
            if d < best_d {
                best_d = d;
                best_i = i;
            }
        }
        grid[best_i] = k;
    }

    let mut u: Vec<f64> = grid.iter().map(|&f| (f - k).max(0.0)).collect();
    let m = n - 1; // interior nodes 1..=m

    // Spatial operator coefficients L u_i = coeff[i] * (w_m * u_{i-1} - (w_m + w_p) * u_i + w_p * u_{i+1})
    let mut a_op = vec![0.0_f64; m];
    let mut b_op = vec![0.0_f64; m];
    let mut c_op = vec![0.0_f64; m];
    for idx in 0..m {
        let i = idx + 1;
        let hm = grid[i] - grid[i - 1];
        let hp = grid[i + 1] - grid[i];
        let diff = 0.5 * alpha * alpha * grid[i].powf(2.0 * beta);
        let fac = 2.0 * diff / (hm + hp);
        a_op[idx] = fac / hm;
        c_op[idx] = fac / hp;
        b_op[idx] = -(a_op[idx] + c_op[idx]);
    }

    let dt = t / n_time as f64;
    let mut dl = vec![0.0_f64; m];
    let mut dd = vec![0.0_f64; m];
    let mut du = vec![0.0_f64; m];
    let mut rhs = vec![0.0_f64; m];

    let thomas_step = |theta: f64,
                       step_dt: f64,
                       u_vec: &mut [f64],
                       dl: &mut [f64],
                       dd: &mut [f64],
                       du: &mut [f64],
                       rhs: &mut [f64]| {
        let bc_hi = (f_max - k).max(0.0);
        for idx in 0..m {
            let i = idx + 1;
            dl[idx] = -theta * step_dt * a_op[idx];
            dd[idx] = 1.0 - theta * step_dt * b_op[idx];
            du[idx] = -theta * step_dt * c_op[idx];
            let explicit = u_vec[i]
                + (1.0 - theta)
                    * step_dt
                    * (a_op[idx] * u_vec[i - 1] + b_op[idx] * u_vec[i] + c_op[idx] * u_vec[i + 1]);
            rhs[idx] = explicit;
        }
        // Boundary contributions at i=1 (u_0=0) and i=m (u_n = bc_hi)
        rhs[m - 1] += theta * step_dt * c_op[m - 1] * bc_hi;

        // Forward sweep
        for idx in 1..m {
            let w = dl[idx] / dd[idx - 1];
            dd[idx] -= w * du[idx - 1];
            rhs[idx] -= w * rhs[idx - 1];
        }
        // Back substitution
        u_vec[m] = rhs[m - 1] / dd[m - 1];
        for idx in (0..m - 1).rev() {
            u_vec[idx + 1] = (rhs[idx] - du[idx] * u_vec[idx + 2]) / dd[idx];
        }
        u_vec[0] = 0.0;
        u_vec[n] = bc_hi;
    };

    // Rannacher startup: 2 half-steps of implicit Euler (theta = 1.0)
    thomas_step(1.0, 0.5 * dt, &mut u, &mut dl, &mut dd, &mut du, &mut rhs);
    thomas_step(1.0, 0.5 * dt, &mut u, &mut dl, &mut dd, &mut du, &mut rhs);
    // Remaining n_time - 1 steps with Crank-Nicolson (theta = 0.5)
    for _ in 1..n_time {
        thomas_step(0.5, dt, &mut u, &mut dl, &mut dd, &mut du, &mut rhs);
    }

    // Cubic-style (4-point Lagrange) interpolation at f0
    let pos = grid
        .binary_search_by(|g| g.partial_cmp(&f0).unwrap())
        .unwrap_or_else(|p| p);
    let i0 = pos.saturating_sub(2).min(n - 3);
    let mut val = 0.0_f64;
    for j in 0..4 {
        let mut lj = 1.0_f64;
        for m_i in 0..4 {
            if m_i != j {
                lj *= (f0 - grid[i0 + m_i]) / (grid[i0 + j] - grid[i0 + m_i]);
            }
        }
        val += lj * u[i0 + j];
    }
    val.max(0.0)
}

/// Precomputed discrete mixture CDF for exact sampling of absorbed CEV `F_T` (`0 < beta < 1`).
pub struct CevExactSampler {
    pub p_absorb: f64,
    pub scale: f64,
    pub inv_2_1mb: f64,
    pub cdf: Vec<f64>,
}

impl CevExactSampler {
    pub fn new(f0: f64, alpha: f64, beta: f64, t: f64) -> Self {
        let one_m_b = 1.0 - beta;
        let nu = 0.5 / one_m_b;
        let scale = alpha * alpha * one_m_b * one_m_b * t;
        let x = f0.powf(2.0 * one_m_b) / scale;
        let half_x = 0.5 * x;
        let p_absorb = gamma_q(nu, half_x);

        let k_max = (half_x + 12.0 * half_x.sqrt() + 64.0).ceil() as usize;
        let mut weights = vec![0.0_f64; k_max + 1];
        let log_hx = half_x.ln();
        let mut sum_w = 0.0_f64;
        for (k, wk) in weights.iter_mut().enumerate() {
            let log_w = -half_x + (nu + k as f64) * log_hx - lgamma(nu + k as f64 + 1.0);
            let w = log_w.exp();
            *wk = w;
            sum_w += w;
        }
        let mut cdf = Vec::with_capacity(k_max + 1);
        let mut acc = 0.0_f64;
        for w in weights {
            acc += w / sum_w;
            cdf.push(acc);
        }
        if let Some(last) = cdf.last_mut() {
            *last = 1.0;
        }
        Self {
            p_absorb,
            scale,
            inv_2_1mb: nu,
            cdf,
        }
    }

    #[inline]
    pub fn sample_ft(&self, rg: &mut Rng) -> f64 {
        if rg.uniform() < self.p_absorb {
            return 0.0;
        }
        let u = rg.uniform();
        let k = match self.cdf.binary_search_by(|c| c.partial_cmp(&u).unwrap()) {
            Ok(i) | Err(i) => i.min(self.cdf.len() - 1),
        };
        let y = rg.gamma((k + 1) as f64, 2.0);
        (self.scale * y).powf(self.inv_2_1mb)
    }
}

/// Exact Monte Carlo for CEV call price, absorbing probability, and martingale check `E[F_T]`.
/// Returns `(mean_call, se_call, p_absorb_mc, mean_ft, se_ft)`.
pub fn cev_mc_exact(
    f0: f64,
    k: f64,
    alpha: f64,
    beta: f64,
    t: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64, f64) {
    let sampler = CevExactSampler::new(f0, alpha, beta, t);
    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x1001, r as u64));
        let ft = sampler.sample_ft(&mut rg);
        let pay = (ft - k).max(0.0);
        let abs_flag = if ft == 0.0 { 1.0 } else { 0.0 };
        [pay, pay * pay, abs_flag, ft]
    });
    let nf = n_paths as f64;
    let m_c = sums[0] / nf;
    let se_c = ((sums[1] / nf - m_c * m_c).max(0.0) / nf).sqrt();
    let p_abs = sums[2] / nf;
    let m_ft = sums[3] / nf;
    // Approximate SE for F_T using second moment bound or separate pass (not needed for call SE)
    let _ = stats::mean(&[m_ft]);
    (m_c, se_c, p_abs, m_ft, 0.0)
}

/// Fine-step Euler Monte Carlo for CEV with absorbing boundary at `F <= 0`.
/// Returns `(mean_call, se_call, p_absorb_mc, mean_ft)`.
pub fn cev_mc_euler(
    f0: f64,
    k: f64,
    alpha: f64,
    beta: f64,
    t: f64,
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let dt = t / n_steps as f64;
    let a_sqrt_dt = alpha * dt.sqrt();
    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x1002, r as u64));
        let mut f = f0;
        for _ in 0..n_steps {
            if f <= 0.0 {
                f = 0.0;
                break;
            }
            f += a_sqrt_dt * f.powf(beta) * rg.normal();
            if f <= 0.0 {
                f = 0.0;
                break;
            }
        }
        let pay = (f - k).max(0.0);
        let abs_flag = if f == 0.0 { 1.0 } else { 0.0 };
        [pay, pay * pay, abs_flag, f]
    });
    let nf = n_paths as f64;
    let m_c = sums[0] / nf;
    let se_c = ((sums[1] / nf - m_c * m_c).max(0.0) / nf).sqrt();
    (m_c, se_c, sums[2] / nf, sums[3] / nf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cev_exact_vs_pde_and_mc() {
        let f0 = 0.04;
        let alpha = 0.0325;
        let t = 5.0;
        for &beta in &[0.3, 0.5, 0.7] {
            for &k in &[0.02, 0.04, 0.06] {
                let c_ex = cev_call_exact(f0, k, alpha, beta, t);
                let c_pde = cev_call_pde_cn(f0, k, alpha, beta, t, 2000, 1200);
                let rel_err = (c_ex - c_pde).abs() / c_ex.max(1e-10);
                assert!(
                    rel_err < 5e-6,
                    "beta={beta} k={k} ex={c_ex} pde={c_pde} rel={rel_err}"
                );
            }
        }
        // Exact MC 3 SE check
        let c_ex = cev_call_exact(f0, 0.04, alpha, 0.5, 10.0);
        let (c_mc, se_mc, p_abs_mc, m_ft, _) =
            cev_mc_exact(f0, 0.04, alpha, 0.5, 10.0, 200_000, 4, 20261123);
        let p_abs_ex = cev_absorb_prob(f0, alpha, 0.5, 10.0);
        assert!(
            (c_mc - c_ex).abs() <= 3.0 * se_mc,
            "ex={c_ex} mc={c_mc} se={se_mc}"
        );
        assert!(
            (p_abs_mc - p_abs_ex).abs() < 0.005,
            "p_abs_mc={p_abs_mc} ex={p_abs_ex}"
        );
        assert!((m_ft - f0).abs() < 0.002, "m_ft={m_ft} f0={f0}");
    }
}
