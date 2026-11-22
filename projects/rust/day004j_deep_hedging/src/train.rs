use crate::bs::bs_delta;
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
    pub init_output_bias: f64,
    pub zero_output: bool,
    pub seed: u64,
}

#[derive(Clone, Debug)]
pub struct TrainResult {
    pub policy: MlpPolicy,
    pub p0: f64,
    pub loss_history: Vec<(usize, f64)>,
}

pub fn train_policy(paths: &[PathData], cfg: TrainConfig, n_workers: usize) -> TrainResult {
    train_policy_with_checkpoints(paths, cfg, &[], n_workers).0
}

pub fn train_policy_with_checkpoints(
    paths: &[PathData],
    cfg: TrainConfig,
    checkpoints: &[usize],
    n_workers: usize,
) -> (TrainResult, Vec<(usize, MlpPolicy, f64)>) {
    let n_total = paths.len();
    let n_steps = paths[0].s.len() - 1;
    let d_in = cfg.spec.d_in();
    let mut pol = MlpPolicy::new(d_in, cfg.hidden, cfg.act, cfg.seed);
    if cfg.zero_output {
        pol.zero_output_layer();
    } else if cfg.init_output_bias.abs() > 0.0 {
        pol.set_output_bias(cfg.init_output_bias);
    }
    let mut p0 = cfg.init_p0;
    let mut opt = Adam::new(pol.n_params() + 1, cfg.lr);
    let mut loss_history = Vec::new();
    let mut snaps = Vec::new();

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
                let tau = cfg.spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
                cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                let nn_d = pol.forward_step(&feat, &mut cache);
                let d = if let Some((r_bs, sig_bs)) = cfg.spec.residual_bs {
                    bs_delta(path.s[t], cfg.spec.k, r_bs, sig_bs, tau) + nn_d
                } else {
                    nn_d
                };
                wealth +=
                    d * (path.s[t + 1] - path.s[t]) - cfg.cost * (d - prev_d).abs() * path.s[t];
                prev_d = d;
            }
            wealth - (path.s[n_steps] - cfg.spec.k).max(0.0)
        });

        let (loss, d_pnl) = loss_and_d_pnl(&pnl, cfg.loss_type);
        let log_stride = (cfg.n_iters / 20).max(20);
        if it == 0 || (it + 1) % log_stride == 0 || it + 1 == cfg.n_iters {
            loss_history.push((it + 1, loss));
        }

        // Pass 2: deterministic fixed-block BPTT gradient accumulation
        let prev_idx = cfg.spec.prev_delta_idx();
        let (_s, mut grad) =
            reduce_vec_fixed_blocks(bsize, pol.n_params(), n_workers, |b, local_g| {
                let path = &paths[(start + b) % n_total];
                let mut caches = Vec::with_capacity(n_steps);
                let mut d_actuals = Vec::with_capacity(n_steps);
                let mut feat = vec![0.0_f64; d_in];
                let mut prev_d = 0.0_f64;
                for t in 0..n_steps {
                    let tau = cfg.spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
                    cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                    let mut c = StepCache::new(d_in, cfg.hidden);
                    let nn_d = pol.forward_step(&feat, &mut c);
                    let d = if let Some((r_bs, sig_bs)) = cfg.spec.residual_bs {
                        bs_delta(path.s[t], cfg.spec.k, r_bs, sig_bs, tau) + nn_d
                    } else {
                        nn_d
                    };
                    prev_d = d;
                    d_actuals.push(d);
                    caches.push(c);
                }
                let dp = d_pnl[b];
                let mut d_in_buf = vec![0.0_f64; d_in];
                let mut dz2 = vec![0.0_f64; cfg.hidden];
                let mut dz1 = vec![0.0_f64; cfg.hidden];
                let mut d_prev_delta_next = 0.0_f64;
                for t in (0..n_steps).rev() {
                    let c = &caches[t];
                    let d_cur = d_actuals[t];
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

        if checkpoints.contains(&(it + 1)) {
            snaps.push((it + 1, pol.clone(), p0));
        }
    }

    (
        TrainResult {
            policy: pol,
            p0,
            loss_history,
        },
        snaps,
    )
}
