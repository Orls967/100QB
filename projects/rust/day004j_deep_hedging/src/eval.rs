use crate::dp::DpResult;
use crate::fourier::heston_delta_fourier;
use crate::hedge::{FeatureSpec, PathData};
use crate::nn::{loss_and_d_pnl, Activation, LossType, MlpPolicy, StepCache};
use crate::optim::Adam;
use crate::parallel::{reduce_vec_fixed_blocks, run_chunks};
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[derive(Clone, Copy, Debug)]
pub struct PairedDiff {
    pub mean: f64,
    pub se: f64,
    pub z: f64,
}

pub fn paired_diff(a: &[f64], b: &[f64]) -> PairedDiff {
    let diffs: Vec<f64> = a.iter().zip(b.iter()).map(|(&x, &y)| x - y).collect();
    let m = stats::mean(&diffs);
    let se = stats::se(&diffs);
    let z = if se > 1e-15 { m / se } else { 0.0 };
    PairedDiff { mean: m, se, z }
}

/// Ordinary least squares slope `b1` of `y ~ b0 + b1 * x`.
pub fn ols_slope(x: &[f64], y: &[f64]) -> f64 {
    let mx = stats::mean(x);
    let my = stats::mean(y);
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mx;
        num += dx * (yi - my);
        den += dx * dx;
    }
    if den > 1e-18 {
        num / den
    } else {
        0.0
    }
}

/// Measures the NN no-trade band half-width at fixed `(t, S)` by sweeping `delta_{t-1}`
/// on a fine `>= 401`-point grid (`step = 0.00125`) and finding the half-width where `|delta_t - delta_{t-1}| <= 0.005`.
pub fn measure_nn_notrade_half_width(
    policy: &MlpPolicy,
    spec: FeatureSpec,
    t: usize,
    n_steps: usize,
    s: f64,
    d_center: f64,
) -> (f64, Vec<f64>, Vec<f64>) {
    let m = 401_usize;
    let d_lo = (d_center - 0.25).max(-0.1);
    let d_hi = (d_center + 0.25).min(1.1);
    let tau = spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
    let mut d_prevs = Vec::with_capacity(m);
    let mut d_outs = Vec::with_capacity(m);
    let mut cache = StepCache::new(policy.d_in, policy.hidden);
    let mut feat = vec![0.0_f64; spec.d_in()];
    feat[0] = if spec.scale_inputs {
        (s / spec.k).ln() / 0.1
    } else {
        (s / spec.k).ln()
    };
    feat[1] = if spec.scale_inputs {
        tau / spec.t_mat
    } else {
        tau
    };

    let mut in_band_min = f64::INFINITY;
    let mut in_band_max = f64::NEG_INFINITY;
    for i in 0..m {
        let dp = d_lo + (d_hi - d_lo) * i as f64 / (m - 1) as f64;
        if let Some(pidx) = spec.prev_delta_idx() {
            feat[pidx] = dp;
        }
        let d_out = policy.forward_step(&feat, &mut cache);
        d_prevs.push(dp);
        d_outs.push(d_out);
        if (d_out - dp).abs() < 0.005 {
            in_band_min = in_band_min.min(dp);
            in_band_max = in_band_max.max(dp);
        }
    }
    let hw = if in_band_max > in_band_min {
        0.5 * (in_band_max - in_band_min)
    } else {
        0.0
    };
    (hw, d_prevs, d_outs)
}

/// Evaluates DP vs NN policy agreement weighted by empirical state visitation frequency along simulated paths.
pub fn visitation_weighted_dp_agreement(
    paths: &[PathData],
    dp: &DpResult,
    policy: &MlpPolicy,
    spec: FeatureSpec,
    half_grid: f64,
    n_workers: usize,
) -> (f64, f64) {
    let n_paths = paths.len();
    let n_steps = dp.n_steps;
    let per_path: Vec<(usize, f64)> = run_chunks(n_paths, n_workers, |r_idx| {
        let path = &paths[r_idx];
        let mut cache = StepCache::new(policy.d_in, policy.hidden);
        let mut feat = vec![0.0_f64; spec.d_in()];
        let mut prev_d = 0.0_f64;
        let mut agree = 0_usize;
        let mut sum_err = 0.0_f64;
        for t in 0..n_steps {
            let d_dp = dp.act(t, path.s[t], prev_d);
            spec.fill_features(path, t, n_steps, prev_d, &mut feat);
            let d_nn = policy.forward_step(&feat, &mut cache);
            let err = (d_dp - d_nn).abs();
            if err <= half_grid {
                agree += 1;
            }
            sum_err += err;
            prev_d = d_dp;
        }
        (agree, sum_err)
    });
    let total_steps = (n_paths * n_steps) as f64;
    let total_agree: usize = per_path.iter().map(|x| x.0).sum();
    let total_err: f64 = per_path.iter().map(|x| x.1).sum();
    (total_agree as f64 / total_steps, total_err / total_steps)
}

/// Verifies 20-iteration training bit-identity between `1` and `8` workers,
/// returning `(n_params, bit_diff_count)`.
pub fn measure_worker_invariance_1_vs_8() -> (usize, usize) {
    let run_with = |n_workers: usize| -> Vec<f64> {
        let n_paths = 256_usize;
        let n_steps = 10_usize;
        let d_in = 3_usize;
        let hidden = 16_usize;
        let seed = 20261122_u64;
        let paths: Vec<Vec<f64>> = run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x1001, r as u64));
            let mut s = 1.0_f64;
            let mut out = Vec::with_capacity(n_steps + 1);
            out.push(s);
            let dt = 0.25 / n_steps as f64;
            for _ in 0..n_steps {
                s *= (-0.5 * 0.04 * dt + 0.2 * dt.sqrt() * rg.normal()).exp();
                out.push(s);
            }
            out
        });
        let mut pol = MlpPolicy::new(d_in, hidden, Activation::Tanh, seed ^ 0x99);
        let mut opt = Adam::new(pol.n_params(), 1e-3);
        let p0 = 0.04_f64;
        let cost = 0.001_f64;
        for _ in 0..20 {
            let pnl = run_chunks(n_paths, n_workers, |r| {
                let path = &paths[r];
                let mut cache = StepCache::new(d_in, hidden);
                let mut prev_d = 0.0_f64;
                let mut wealth = p0;
                for (t, &s_t) in path.iter().enumerate().take(n_steps) {
                    let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
                    let d = pol.forward_step(&[s_t.ln(), tau, prev_d], &mut cache);
                    wealth += d * (path[t + 1] - s_t) - cost * (d - prev_d).abs() * s_t;
                    prev_d = d;
                }
                wealth - (path[n_steps] - 1.0).max(0.0)
            });
            let (_loss, d_pnl) = loss_and_d_pnl(&pnl, LossType::Entropic { lambda: 10.0 });
            let (_s, grad) =
                reduce_vec_fixed_blocks(n_paths, pol.n_params(), n_workers, |r, local_g| {
                    let path = &paths[r];
                    let mut caches = Vec::with_capacity(n_steps);
                    let mut prev_d = 0.0_f64;
                    for (t, &s_t) in path.iter().enumerate().take(n_steps) {
                        let tau = 0.25 * (n_steps - t) as f64 / n_steps as f64;
                        let mut c = StepCache::new(d_in, hidden);
                        let d = pol.forward_step(&[s_t.ln(), tau, prev_d], &mut c);
                        prev_d = d;
                        caches.push(c);
                    }
                    let dp = d_pnl[r];
                    let mut d_in_buf = vec![0.0_f64; d_in];
                    let mut dz2 = vec![0.0_f64; hidden];
                    let mut dz1 = vec![0.0_f64; hidden];
                    let mut d_prev_delta_next = 0.0_f64;
                    for t in (0..n_steps).rev() {
                        let c = &caches[t];
                        let d_cur = c.out;
                        let d_prev = c.x[2];
                        let sgn = if d_cur >= d_prev { 1.0 } else { -1.0 };
                        let d_direct = dp * ((path[t + 1] - path[t]) - cost * sgn * path[t]);
                        pol.backward_step(
                            c,
                            d_direct + d_prev_delta_next,
                            local_g,
                            &mut d_in_buf,
                            &mut dz2,
                            &mut dz1,
                        );
                        d_prev_delta_next = d_in_buf[2] + dp * cost * sgn * path[t];
                    }
                    0.0
                });
            opt.step(&mut pol.params, &grad);
        }
        pol.params
    };
    let p1 = run_with(1);
    let p8 = run_with(8);
    let diffs = p1
        .iter()
        .zip(p8.iter())
        .filter(|(&a, &b)| a.to_bits() != b.to_bits())
        .count();
    (p1.len(), diffs)
}

/// Precomputed 3D table `(t, x = ln(S/K), v)` for fast Heston Fourier delta evaluation across 100k paths.
pub struct HestonDeltaTable {
    pub n_steps: usize,
    pub x_min: f64,
    pub x_max: f64,
    pub nx: usize,
    pub v_min: f64,
    pub v_max: f64,
    pub nv: usize,
    pub table: Vec<Vec<f64>>,
}

impl HestonDeltaTable {
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        k: f64,
        r: f64,
        kappa: f64,
        theta: f64,
        xi: f64,
        rho: f64,
        t_mat: f64,
        n_steps: usize,
        n_workers: usize,
    ) -> Self {
        let nx = 45_usize;
        let nv = 25_usize;
        let x_min = -0.55_f64;
        let x_max = 0.55_f64;
        let v_min = 0.002_f64;
        let v_max = 0.25_f64;

        let table = run_chunks(n_steps, n_workers, |t| {
            let tau = t_mat * (n_steps - t) as f64 / n_steps as f64;
            let mut slice = vec![0.0_f64; nx * nv];
            for ix in 0..nx {
                let x = x_min + (x_max - x_min) * ix as f64 / (nx - 1) as f64;
                let s = k * x.exp();
                for iv in 0..nv {
                    let v = v_min + (v_max - v_min) * iv as f64 / (nv - 1) as f64;
                    slice[ix * nv + iv] =
                        heston_delta_fourier(s, k, r, v, kappa, theta, xi, rho, tau);
                }
            }
            slice
        });

        Self {
            n_steps,
            x_min,
            x_max,
            nx,
            v_min,
            v_max,
            nv,
            table,
        }
    }

    pub fn lookup(&self, t: usize, s: f64, k: f64, v: f64) -> f64 {
        let x = (s / k).ln().clamp(self.x_min, self.x_max);
        let vc = v.clamp(self.v_min, self.v_max);
        let tx = ((x - self.x_min) / (self.x_max - self.x_min)) * (self.nx - 1) as f64;
        let tv = ((vc - self.v_min) / (self.v_max - self.v_min)) * (self.nv - 1) as f64;
        let ix0 = (tx.floor() as usize).min(self.nx - 2);
        let iv0 = (tv.floor() as usize).min(self.nv - 2);
        let fx = tx - ix0 as f64;
        let fv = tv - iv0 as f64;
        let sl = &self.table[t];
        let d00 = sl[ix0 * self.nv + iv0];
        let d01 = sl[ix0 * self.nv + iv0 + 1];
        let d10 = sl[(ix0 + 1) * self.nv + iv0];
        let d11 = sl[(ix0 + 1) * self.nv + iv0 + 1];
        (1.0 - fx) * ((1.0 - fv) * d00 + fv * d01) + fx * ((1.0 - fv) * d10 + fv * d11)
    }
}
