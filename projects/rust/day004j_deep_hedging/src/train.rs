use crate::hedge::{FeatureSpec, PathData};
use crate::nn::{loss_and_d_pnl, Activation, LossType, MlpPolicy, StepCache};
use crate::optim::Adam;
use crate::parallel::{reduce_vec_fixed_blocks, run_chunks};

#[derive(Clone, Copy, Debug)]
pub struct TrainConfig {
    pub hidden: usize,
    pub act: Activation,
    pub n_iters: usize,
    pub batch_size: usize,
    pub lr: f64,
    pub lr_decay_step: usize,
    pub lr_decay_gamma: f64,
    pub cost: f64,
    pub loss_type: LossType,
    pub spec: FeatureSpec,
    pub init_p0: f64,
    pub train_p0: bool,
    pub seed: u64,
}

#[derive(Clone, Debug)]
pub struct TrainResult {
    pub policy: MlpPolicy,
    pub p0: f64,
    pub loss_history: Vec<(usize, f64)>,
}

pub fn train_policy(paths: &[PathData], cfg: TrainConfig, n_workers: usize) -> TrainResult {
    let n_total = paths.len();
    let n_steps = paths[0].s.len() - 1;
    let d_in = cfg.spec.d_in();
    let mut pol = MlpPolicy::new(d_in, cfg.hidden, cfg.act, cfg.seed);
    let mut p0 = cfg.init_p0;
    let mut opt = Adam::new(pol.n_params() + 1, cfg.lr);
    let mut loss_history = Vec::new();

    for it in 0..cfg.n_iters {
        if cfg.lr_decay_step > 0 && it > 0 && it % cfg.lr_decay_step == 0 {
            opt.lr *= cfg.lr_decay_gamma;
        }
        let start = (it * cfg.batch_size) % n_total;
        let bsize = cfg.batch_size.min(n_total);

        // Pass 1: forward PnL on mini-batch
        let pnl = run_chunks(bsize, n_workers, |b| {
            let path = &paths[(start + b) % n_total];
            let mut cache = StepCache::new(d_in, cfg.hidden);
            let mut feat = vec![0.0_f64; d_in];
            let mut prev_d = 0.0_f64;
            let mut wealth = p0;
            for t in 0..n_steps {
                cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                let d = pol.forward_step(&feat, &mut cache);
                wealth +=
                    d * (path.s[t + 1] - path.s[t]) - cfg.cost * (d - prev_d).abs() * path.s[t];
                prev_d = d;
            }
            wealth - (path.s[n_steps] - cfg.spec.k).max(0.0)
        });

        let (loss, d_pnl) = loss_and_d_pnl(&pnl, cfg.loss_type);
        if it == 0 || (it + 1) % 20 == 0 || it + 1 == cfg.n_iters {
            loss_history.push((it + 1, loss));
        }

        // Pass 2: deterministic fixed-block BPTT gradient accumulation
        let prev_idx = cfg.spec.prev_delta_idx();
        let (_s, mut grad) =
            reduce_vec_fixed_blocks(bsize, pol.n_params(), n_workers, |b, local_g| {
                let path = &paths[(start + b) % n_total];
                let mut caches = Vec::with_capacity(n_steps);
                let mut feat = vec![0.0_f64; d_in];
                let mut prev_d = 0.0_f64;
                for t in 0..n_steps {
                    cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                    let mut c = StepCache::new(d_in, cfg.hidden);
                    let d = pol.forward_step(&feat, &mut c);
                    prev_d = d;
                    caches.push(c);
                }
                let dp = d_pnl[b];
                let mut d_in_buf = vec![0.0_f64; d_in];
                let mut dz2 = vec![0.0_f64; cfg.hidden];
                let mut dz1 = vec![0.0_f64; cfg.hidden];
                let mut d_prev_delta_next = 0.0_f64;
                for t in (0..n_steps).rev() {
                    let c = &caches[t];
                    let d_cur = c.out;
                    let d_prev = prev_idx.map_or(0.0, |idx| c.x[idx]);
                    let sgn = if d_cur >= d_prev { 1.0 } else { -1.0 };
                    let d_direct = dp * ((path.s[t + 1] - path.s[t]) - cfg.cost * sgn * path.s[t]);
                    pol.backward_step(
                        c,
                        d_direct + d_prev_delta_next,
                        local_g,
                        &mut d_in_buf,
                        &mut dz2,
                        &mut dz1,
                    );
                    d_prev_delta_next = if let Some(idx) = prev_idx {
                        d_in_buf[idx] + dp * cfg.cost * sgn * path.s[t]
                    } else {
                        0.0
                    };
                }
                0.0
            });

        let dp0 = if cfg.train_p0 {
            d_pnl.iter().sum::<f64>()
        } else {
            0.0
        };
        grad.push(dp0);
        let mut params_ext = Vec::with_capacity(pol.n_params() + 1);
        params_ext.extend_from_slice(&pol.params);
        params_ext.push(p0);
        opt.step(&mut params_ext, &grad);
        let np = pol.n_params();
        pol.params.copy_from_slice(&params_ext[..np]);
        p0 = params_ext[np];
    }

    TrainResult {
        policy: pol,
        p0,
        loss_history,
    }
}
