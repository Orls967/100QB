use crate::fourier::heston_delta_fourier;
use crate::hedge::FeatureSpec;
use crate::nn::{MlpPolicy, StepCache};
use crate::parallel::run_chunks;
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
/// in `[d_center - 0.35, d_center + 0.35]` and finding the half-width where `|delta_t - delta_{t-1}| <= tol_hold`.
pub fn measure_nn_notrade_half_width(
    policy: &MlpPolicy,
    spec: FeatureSpec,
    t: usize,
    n_steps: usize,
    s: f64,
    d_center: f64,
) -> (f64, Vec<f64>, Vec<f64>) {
    let m = 71_usize;
    let d_lo = (d_center - 0.25).max(-0.1);
    let d_hi = (d_center + 0.25).min(1.1);
    let tau = spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
    let mut d_prevs = Vec::with_capacity(m);
    let mut d_outs = Vec::with_capacity(m);
    let mut cache = StepCache::new(policy.d_in, policy.hidden);
    let mut feat = vec![0.0_f64; spec.d_in()];
    feat[0] = (s / spec.k).ln();
    feat[1] = tau;

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
        if (d_out - dp).abs() < 0.012 {
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
