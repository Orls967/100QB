use crate::bs::bs_call;
use crate::parallel::reduce_vec_fixed_blocks;
use crate::rng::{derive_seed, Rng};

#[derive(Clone, Copy, Debug)]
pub struct MertonParams {
    pub s0: f64,
    pub r: f64,
    pub q: f64,
    pub sigma: f64,
    pub lambda: f64,
    pub mu_j: f64,
    pub delta_j: f64,
}

impl Default for MertonParams {
    fn default() -> Self {
        Self {
            s0: 100.0,
            r: 0.05,
            q: 0.01,
            sigma: 0.15,
            lambda: 0.80,
            mu_j: -0.10,
            delta_j: 0.25,
        }
    }
}

impl MertonParams {
    /// Expected relative jump size `kappa = E[e^J - 1] = exp(mu_j + 0.5 * delta_j^2) - 1`.
    #[inline]
    pub fn kappa(&self) -> f64 {
        (self.mu_j + 0.5 * self.delta_j * self.delta_j).exp() - 1.0
    }
}

/// Merton (1976) exact series for European call option, truncated at `n_max = 60`
/// (Poisson tail < 1e-16 for moderate lambda*T).
pub fn merton_call_series(p: &MertonParams, k: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return (p.s0 - k).max(0.0);
    }
    if p.lambda <= 0.0 {
        return bs_call(p.s0, k, p.r, p.q, p.sigma, t);
    }
    let kappa = p.kappa();
    let lam_t = p.lambda * t;
    let mut weight = (-lam_t).exp();
    let mut total = 0.0_f64;
    let sig2 = p.sigma * p.sigma;
    let dj2 = p.delta_j * p.delta_j;
    let mean_log_j = p.mu_j + 0.5 * dj2;

    for n in 0..=80 {
        let nf = n as f64;
        let s_n = p.s0 * (nf * mean_log_j - lam_t * kappa).exp();
        let sigma_n = (sig2 + nf * dj2 / t).sqrt();
        total += weight * bs_call(s_n, k, p.r, p.q, sigma_n, t);
        if n > 5 && weight < 1e-17 {
            break;
        }
        weight *= lam_t / ((n + 1) as f64);
    }
    total
}

/// Merton European put price via exact put-call parity.
pub fn merton_put_series(p: &MertonParams, k: f64, t: f64) -> f64 {
    let c = merton_call_series(p, k, t);
    c - p.s0 * (-p.q * t).exp() + k * (-p.r * t).exp()
}

#[derive(Clone, Debug)]
pub struct MertonMcRow {
    pub k: f64,
    pub call_exact: f64,
    pub call_mc: f64,
    pub call_se: f64,
    pub call_z: f64,
    pub put_exact: f64,
    pub put_mc: f64,
    pub put_se: f64,
    pub put_z: f64,
}

#[derive(Clone, Debug)]
pub struct MertonMcSummary {
    pub t: f64,
    pub r_paths: usize,
    pub martingale_mean: f64,
    pub martingale_se: f64,
    pub martingale_z: f64,
    pub rows: Vec<MertonMcRow>,
}

/// Exact terminal simulation of Merton jump-diffusion (`N_T ~ Poisson(lambda*T)`,
/// `sum_{i=1}^{N_T} J_i ~ N(N_T * mu_j, N_T * delta_j^2)`).
pub fn simulate_merton_mc(
    p: &MertonParams,
    ks: &[f64],
    t: f64,
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> MertonMcSummary {
    let kappa = p.kappa();
    let lam_t = p.lambda * t;
    let drift = (p.r - p.q - 0.5 * p.sigma * p.sigma - p.lambda * kappa) * t;
    let vol_sqrt_t = p.sigma * t.sqrt();
    let disc_r = (-p.r * t).exp();
    let disc_rq = (-(p.r - p.q) * t).exp();
    let n_k = ks.len();

    // 2 accumulators for martingale S_T * e^{-(r-q)T} + 4 per strike (call sum/sq, put sum/sq)
    let sums = reduce_vec_fixed_blocks(r_paths, n_workers, 2 + n_k * 4, |idx, acc| {
        let mut rng = Rng::new(derive_seed(main_seed, 100, idx as u64));
        let z_w = rng.normal();
        let n_jumps = rng.poisson(lam_t);
        let jump_sum = if n_jumps == 0 {
            0.0
        } else {
            let nf = n_jumps as f64;
            nf * p.mu_j + (nf.sqrt() * p.delta_j) * rng.normal()
        };
        let s_t = p.s0 * (drift + vol_sqrt_t * z_w + jump_sum).exp();
        let disc_st = disc_rq * s_t;
        acc[0] += disc_st;
        acc[1] += disc_st * disc_st;

        for (ki, &k) in ks.iter().enumerate() {
            let c = disc_r * (s_t - k).max(0.0);
            let pt = disc_r * (k - s_t).max(0.0);
            let off = 2 + ki * 4;
            acc[off] += c;
            acc[off + 1] += c * c;
            acc[off + 2] += pt;
            acc[off + 3] += pt * pt;
        }
    });

    let rf = r_paths as f64;
    let stat = |sum: f64, sum_sq: f64| -> (f64, f64) {
        let m = sum / rf;
        let se = ((sum_sq / rf - m * m).max(0.0) / rf).sqrt();
        (m, se)
    };
    let (martingale_mean, martingale_se) = stat(sums[0], sums[1]);
    let martingale_z = (martingale_mean - p.s0) / martingale_se.max(1e-15);

    let rows = ks
        .iter()
        .enumerate()
        .map(|(ki, &k)| {
            let off = 2 + ki * 4;
            let (call_mc, call_se) = stat(sums[off], sums[off + 1]);
            let (put_mc, put_se) = stat(sums[off + 2], sums[off + 3]);
            let call_exact = merton_call_series(p, k, t);
            let put_exact = merton_put_series(p, k, t);
            MertonMcRow {
                k,
                call_exact,
                call_mc,
                call_se,
                call_z: (call_mc - call_exact) / call_se.max(1e-15),
                put_exact,
                put_mc,
                put_se,
                put_z: (put_mc - put_exact) / put_se.max(1e-15),
            }
        })
        .collect();

    MertonMcSummary {
        t,
        r_paths,
        martingale_mean,
        martingale_se,
        martingale_z,
        rows,
    }
}
