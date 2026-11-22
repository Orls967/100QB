use crate::bs::{bs_delta, leland_sigma, whalley_wilmott_half_width};
use crate::dp::DpResult;
use crate::nn::{MlpPolicy, StepCache};
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[derive(Clone, Debug)]
pub struct PathData {
    pub s: Vec<f64>,
    pub v: Option<Vec<f64>>,
}

#[derive(Clone, Copy, Debug)]
pub struct GbmParams {
    pub s0: f64,
    pub k: f64,
    pub r: f64,
    pub sigma: f64,
    pub mu: f64,
    pub t_mat: f64,
    pub n_steps: usize,
}

pub fn generate_gbm_paths(
    p: &GbmParams,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    tag: u64,
) -> Vec<PathData> {
    let dt = p.t_mat / p.n_steps as f64;
    let drift = (p.mu - 0.5 * p.sigma * p.sigma) * dt;
    let vol_dt = p.sigma * dt.sqrt();
    run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, tag, r as u64));
        let mut s = Vec::with_capacity(p.n_steps + 1);
        let mut cur = p.s0;
        s.push(cur);
        for _ in 0..p.n_steps {
            cur *= (drift + vol_dt * rg.normal()).exp();
            s.push(cur);
        }
        PathData { s, v: None }
    })
}

#[derive(Clone, Copy, Debug)]
pub struct FeatureSpec {
    pub k: f64,
    pub t_mat: f64,
    pub use_v: bool,
    pub use_rv5: bool,
    pub use_prev_delta: bool,
    pub scale_inputs: bool,
    /// Optional `(r, sigma)` baseline for residual hedging `delta = bs_delta + NN(...)`.
    pub residual_bs: Option<(f64, f64)>,
}

impl FeatureSpec {
    pub fn basic(k: f64, t_mat: f64, use_v: bool, use_rv5: bool, use_prev_delta: bool) -> Self {
        Self {
            k,
            t_mat,
            use_v,
            use_rv5,
            use_prev_delta,
            scale_inputs: false,
            residual_bs: None,
        }
    }

    pub fn d_in(self) -> usize {
        2 + usize::from(self.use_v || self.use_rv5) + usize::from(self.use_prev_delta)
    }

    #[inline]
    pub fn prev_delta_idx(self) -> Option<usize> {
        if self.use_prev_delta {
            Some(2 + usize::from(self.use_v || self.use_rv5))
        } else {
            None
        }
    }

    pub fn fill_features(
        self,
        path: &PathData,
        t: usize,
        n_steps: usize,
        prev_delta: f64,
        out: &mut [f64],
    ) {
        let tau = self.t_mat * (n_steps - t) as f64 / n_steps as f64;
        let log_m = (path.s[t] / self.k).ln();
        if self.scale_inputs {
            out[0] = log_m / 0.1;
            out[1] = tau / self.t_mat;
        } else {
            out[0] = log_m;
            out[1] = tau;
        }
        let mut idx = 2;
        if self.use_v {
            let vt = path.v.as_ref().map_or(0.04, |vv| vv[t]);
            out[idx] = vt.max(0.0).sqrt();
            idx += 1;
        } else if self.use_rv5 {
            let dt = self.t_mat / n_steps as f64;
            let rv = if t == 0 {
                0.2
            } else {
                let start = t.saturating_sub(5);
                let cnt = (t - start) as f64;
                let mut ss = 0.0_f64;
                for j in start..t {
                    let lr = (path.s[j + 1] / path.s[j]).ln();
                    ss += lr * lr;
                }
                (ss / (cnt * dt)).sqrt()
            };
            out[idx] = rv;
            idx += 1;
        }
        if self.use_prev_delta {
            out[idx] = prev_delta;
        }
    }
}

pub enum Strategy<'a> {
    BsDelta {
        k: f64,
        r: f64,
        sigma: f64,
        t_mat: f64,
    },
    WhalleyWilmott {
        k: f64,
        r: f64,
        sigma: f64,
        t_mat: f64,
        c: f64,
        lambda: f64,
    },
    Leland {
        k: f64,
        r: f64,
        sigma: f64,
        t_mat: f64,
        c: f64,
    },
    HestonFourierDelta {
        deltas: &'a [Vec<f64>],
    },
    NeuralNet {
        policy: &'a MlpPolicy,
        spec: FeatureSpec,
    },
    Dp {
        dp: &'a DpResult,
    },
}

#[derive(Clone, Debug)]
pub struct EvalMetrics {
    pub mean_pnl: f64,
    pub sd_pnl: f64,
    pub se_sd_pnl: f64,
    pub mse: f64,
    pub entropic_risk: f64,
    pub indifference_price: f64,
    pub cvar_95: f64,
    pub mean_turnover: f64,
    pub mean_cost: f64,
    pub pnl: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
pub fn eval_strategy(
    paths: &[PathData],
    k: f64,
    t_mat: f64,
    cost: f64,
    p0: f64,
    lambda: f64,
    strat: Strategy<'_>,
    n_workers: usize,
) -> EvalMetrics {
    let n_paths = paths.len();
    let n_steps = paths[0].s.len() - 1;
    let dt = t_mat / n_steps as f64;

    let per_path: Vec<(f64, f64, f64)> = run_chunks(n_paths, n_workers, |r_idx| {
        let path = &paths[r_idx];
        let mut prev_d = 0.0_f64;
        let mut wealth = p0;
        let mut turnover = 0.0_f64;
        let mut total_cost = 0.0_f64;

        let mut cache = match &strat {
            Strategy::NeuralNet { policy, .. } => Some(StepCache::new(policy.d_in, policy.hidden)),
            _ => None,
        };
        let mut feat = match &strat {
            Strategy::NeuralNet { spec, .. } => vec![0.0_f64; spec.d_in()],
            _ => Vec::new(),
        };

        for t in 0..n_steps {
            let tau = t_mat * (n_steps - t) as f64 / n_steps as f64;
            let s_t = path.s[t];
            let d = match &strat {
                Strategy::BsDelta {
                    k: ks,
                    r,
                    sigma,
                    t_mat: _,
                } => bs_delta(s_t, *ks, *r, *sigma, tau),
                Strategy::WhalleyWilmott {
                    k: ks,
                    r,
                    sigma,
                    t_mat: _,
                    c,
                    lambda: lam,
                } => {
                    let db = bs_delta(s_t, *ks, *r, *sigma, tau);
                    let hw = whalley_wilmott_half_width(s_t, *ks, *r, *sigma, tau, *c, *lam);
                    if prev_d < db - hw {
                        db - hw
                    } else if prev_d > db + hw {
                        db + hw
                    } else {
                        prev_d
                    }
                }
                Strategy::Leland {
                    k: ks,
                    r,
                    sigma,
                    t_mat: _,
                    c,
                } => {
                    let sig_l = leland_sigma(*sigma, *c, dt);
                    bs_delta(s_t, *ks, *r, sig_l, tau)
                }
                Strategy::HestonFourierDelta { deltas } => deltas[r_idx][t],
                Strategy::NeuralNet { policy, spec } => {
                    spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                    let nn_d = policy.forward_step(&feat, cache.as_mut().unwrap());
                    if let Some((r_bs, sig_bs)) = spec.residual_bs {
                        bs_delta(s_t, spec.k, r_bs, sig_bs, tau) + nn_d
                    } else {
                        nn_d
                    }
                }
                Strategy::Dp { dp } => dp.act(t, s_t, prev_d),
            };
            let trade = (d - prev_d).abs();
            let step_cost = cost * trade * s_t;
            turnover += trade;
            total_cost += step_cost;
            wealth += d * (path.s[t + 1] - s_t) - step_cost;
            prev_d = d;
        }
        let payoff = (path.s[n_steps] - k).max(0.0);
        (wealth - payoff, turnover, total_cost)
    });

    let pnl: Vec<f64> = per_path.iter().map(|x| x.0).collect();
    let turnovers: Vec<f64> = per_path.iter().map(|x| x.1).collect();
    let costs: Vec<f64> = per_path.iter().map(|x| x.2).collect();

    let mean_pnl = stats::mean(&pnl);
    let sd_pnl = stats::sd(&pnl);
    let se_sd_pnl = sd_pnl / (2.0 * (n_paths.saturating_sub(1).max(1)) as f64).sqrt();
    let mse = pnl.iter().map(|&x| x * x).sum::<f64>() / n_paths as f64;
    let entropic_risk = stats::entropic_risk(&pnl, lambda);
    let indifference_price = p0 + entropic_risk;
    let cvar_95 = stats::cvar_loss(&pnl, 0.95);

    EvalMetrics {
        mean_pnl,
        sd_pnl,
        se_sd_pnl,
        mse,
        entropic_risk,
        indifference_price,
        cvar_95,
        mean_turnover: stats::mean(&turnovers),
        mean_cost: stats::mean(&costs),
        pnl,
    }
}

pub fn delta_l2_vs_bs(
    paths: &[PathData],
    policy: &MlpPolicy,
    spec: FeatureSpec,
    r: f64,
    sigma: f64,
    n_workers: usize,
) -> f64 {
    let n_paths = paths.len();
    let n_steps = paths[0].s.len() - 1;
    let sq_errs = run_chunks(n_paths, n_workers, |r_idx| {
        let path = &paths[r_idx];
        let mut cache = StepCache::new(policy.d_in, policy.hidden);
        let mut feat = vec![0.0_f64; spec.d_in()];
        let mut prev_d = 0.0_f64;
        let mut sum_sq = 0.0_f64;
        for t in 0..n_steps {
            let tau = spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
            spec.fill_features(path, t, n_steps, prev_d, &mut feat);
            let nn_d = policy.forward_step(&feat, &mut cache);
            let d_nn = if let Some((r_bs, sig_bs)) = spec.residual_bs {
                bs_delta(path.s[t], spec.k, r_bs, sig_bs, tau) + nn_d
            } else {
                nn_d
            };
            let d_bs = bs_delta(path.s[t], spec.k, r, sigma, tau);
            let diff = d_nn - d_bs;
            sum_sq += diff * diff;
            prev_d = d_nn;
        }
        sum_sq / n_steps as f64
    });
    (sq_errs.iter().sum::<f64>() / n_paths as f64).sqrt()
}
