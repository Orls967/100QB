use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Barrier, RwLock};

use crate::bs::bs_delta;
use crate::hedge::{FeatureSpec, PathData};
use crate::nn::{loss_and_d_pnl, Activation, LossType, MlpPolicy, StepCache};
use crate::optim::Adam;
#[cfg(test)]
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

    if cfg.n_iters == 0 {
        return (
            TrainResult {
                policy: pol,
                p0,
                loss_history,
            },
            snaps,
        );
    }

    let bsize = cfg.batch_size.min(n_total);
    let n_params = pol.n_params();
    const BLOCK: usize = 64;
    let n_blocks = bsize.div_ceil(BLOCK);
    let prev_idx = cfg.spec.prev_delta_idx();
    let diff_idx = cfg.spec.delta_diff_idx();

    let w = n_workers.max(1);
    if w == 1 {
        let mut pnl = vec![0.0_f64; bsize];
        let mut caches: Vec<StepCache> = (0..n_steps)
            .map(|_| StepCache::new(d_in, cfg.hidden))
            .collect();
        let mut d_actuals = vec![0.0_f64; n_steps];
        let mut feat = vec![0.0_f64; d_in];
        let mut d_in_buf = vec![0.0_f64; d_in];
        let mut dz2 = vec![0.0_f64; cfg.hidden];
        let mut dz1 = vec![0.0_f64; cfg.hidden];
        let mut local_g = vec![0.0_f64; n_params];
        let mut block_acc = vec![0.0_f64; n_params];
        let mut grad = Vec::with_capacity(n_params + 1);
        let mut params_ext = Vec::with_capacity(n_params + 1);

        for it in 0..cfg.n_iters {
            if cfg.lr_decay_step > 0 && it > 0 && it % cfg.lr_decay_step == 0 {
                opt.lr *= cfg.lr_decay_gamma;
            }
            let start = (it * cfg.batch_size) % n_total;

            for (b, pnl_slot) in pnl.iter_mut().enumerate() {
                let path = &paths[(start + b) % n_total];
                let mut prev_d = 0.0_f64;
                let mut wealth = p0;
                for t in 0..n_steps {
                    let tau = cfg.spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
                    cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                    let nn_d = pol.forward_step(&feat, &mut caches[0]);
                    let d = if let Some((r_bs, sig_bs)) = cfg.spec.residual_bs {
                        bs_delta(path.s[t], cfg.spec.k, r_bs, sig_bs, tau) + nn_d
                    } else {
                        nn_d
                    };
                    wealth +=
                        d * (path.s[t + 1] - path.s[t]) - cfg.cost * (d - prev_d).abs() * path.s[t];
                    prev_d = d;
                }
                *pnl_slot = wealth - (path.s[n_steps] - cfg.spec.k).max(0.0);
            }

            let (loss, d_pnl) = loss_and_d_pnl(&pnl, cfg.loss_type);
            let log_stride = (cfg.n_iters / 20).max(20);
            if it == 0 || (it + 1) % log_stride == 0 || it + 1 == cfg.n_iters {
                loss_history.push((it + 1, loss));
            }

            grad.clear();
            grad.resize(n_params, 0.0);
            for blk in 0..n_blocks {
                let b_start = blk * BLOCK;
                let b_end = ((blk + 1) * BLOCK).min(bsize);
                block_acc.fill(0.0);
                for (b, &dp) in d_pnl.iter().enumerate().take(b_end).skip(b_start) {
                    let path = &paths[(start + b) % n_total];
                    let mut prev_d = 0.0_f64;
                    for t in 0..n_steps {
                        let tau = cfg.spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
                        cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                        let nn_d = pol.forward_step(&feat, &mut caches[t]);
                        let d = if let Some((r_bs, sig_bs)) = cfg.spec.residual_bs {
                            bs_delta(path.s[t], cfg.spec.k, r_bs, sig_bs, tau) + nn_d
                        } else {
                            nn_d
                        };
                        prev_d = d;
                        d_actuals[t] = d;
                    }
                    local_g.fill(0.0);
                    let mut d_prev_delta_next = 0.0_f64;
                    for t in (0..n_steps).rev() {
                        let c = &caches[t];
                        let d_cur = d_actuals[t];
                        let d_prev = prev_idx.map_or(0.0, |idx| c.x[idx]);
                        let sgn = if d_cur >= d_prev { 1.0 } else { -1.0 };
                        let d_direct =
                            dp * ((path.s[t + 1] - path.s[t]) - cfg.cost * sgn * path.s[t]);
                        pol.backward_step(
                            c,
                            d_direct + d_prev_delta_next,
                            &mut local_g,
                            &mut d_in_buf,
                            &mut dz2,
                            &mut dz1,
                        );
                        let d_from_prev = prev_idx.map_or(0.0, |idx| d_in_buf[idx]);
                        let d_from_diff = diff_idx.map_or(0.0, |idx| d_in_buf[idx]);
                        d_prev_delta_next = if prev_idx.is_some() || diff_idx.is_some() {
                            d_from_prev + d_from_diff + dp * cfg.cost * sgn * path.s[t]
                        } else {
                            0.0
                        };
                    }
                    for j in 0..n_params {
                        block_acc[j] += local_g[j];
                    }
                }
                for j in 0..n_params {
                    grad[j] += block_acc[j];
                }
            }

            let dp0 = if cfg.train_p0 {
                d_pnl.iter().sum::<f64>()
            } else {
                0.0
            };
            grad.push(dp0);
            params_ext.clear();
            params_ext.extend_from_slice(&pol.params);
            params_ext.push(p0);
            opt.step(&mut params_ext, &grad);
            pol.params.copy_from_slice(&params_ext[..n_params]);
            p0 = params_ext[n_params];

            if checkpoints.contains(&(it + 1)) {
                snaps.push((it + 1, pol.clone(), p0));
            }
        }

        return (
            TrainResult {
                policy: pol,
                p0,
                loss_history,
            },
            snaps,
        );
    }

    // Persistent worker pool for n_workers > 1 (single thread::scope for all iterations)
    struct SharedBufs {
        params: Vec<f64>,
        d_pnl: Vec<f64>,
    }

    let w2 = w.min(n_blocks);
    let blocks_per_w = n_blocks.div_ceil(w2);

    let phase = AtomicUsize::new(0); // 1=pass1, 2=pass2, 3=exit
    let start_idx = AtomicUsize::new(0);
    let p0_bits = AtomicU64::new(p0.to_bits());
    let shared_bufs = RwLock::new(SharedBufs {
        params: pol.params.clone(),
        d_pnl: vec![0.0_f64; bsize],
    });
    let bar_start = Barrier::new(w + 1);
    let bar_end = Barrier::new(w + 1);

    let worker_pnl: Vec<RwLock<Vec<f64>>> = (0..w)
        .map(|tid| {
            let len = if tid < w2 {
                let bs = tid * blocks_per_w * BLOCK;
                let be = ((tid + 1) * blocks_per_w * BLOCK).min(bsize);
                be.saturating_sub(bs)
            } else {
                0
            };
            RwLock::new(vec![0.0_f64; len])
        })
        .collect();

    let worker_grads: Vec<RwLock<Vec<Vec<f64>>>> = (0..w)
        .map(|tid| {
            let cnt = if tid < w2 {
                let bs = tid * blocks_per_w;
                let be = ((tid + 1) * blocks_per_w).min(n_blocks);
                be.saturating_sub(bs)
            } else {
                0
            };
            RwLock::new((0..cnt).map(|_| vec![0.0_f64; n_params]).collect())
        })
        .collect();

    std::thread::scope(|scope| {
        let phase_ref = &phase;
        let start_ref = &start_idx;
        let p0_ref = &p0_bits;
        let bufs_ref = &shared_bufs;
        let bs_ref = &bar_start;
        let be_ref = &bar_end;
        let wpnl_ref = &worker_pnl;
        let wgrad_ref = &worker_grads;

        for tid in 0..w {
            scope.spawn(move || {
                let mut local_pol = MlpPolicy {
                    d_in,
                    hidden: cfg.hidden,
                    act: cfg.act,
                    params: vec![0.0_f64; n_params],
                };
                let max_paths_w = if tid < w2 {
                    let blk_start = tid * blocks_per_w;
                    let blk_end = ((tid + 1) * blocks_per_w).min(n_blocks);
                    (blk_end - blk_start) * BLOCK
                } else {
                    0
                };
                let mut path_caches: Vec<Vec<StepCache>> = (0..max_paths_w)
                    .map(|_| {
                        (0..n_steps)
                            .map(|_| StepCache::new(d_in, cfg.hidden))
                            .collect()
                    })
                    .collect();
                let mut path_d_actuals: Vec<Vec<f64>> =
                    (0..max_paths_w).map(|_| vec![0.0_f64; n_steps]).collect();
                let mut feat = vec![0.0_f64; d_in];
                let mut d_in_buf = vec![0.0_f64; d_in];
                let mut dz2 = vec![0.0_f64; cfg.hidden];
                let mut dz1 = vec![0.0_f64; cfg.hidden];
                let mut local_g = vec![0.0_f64; n_params];

                loop {
                    bs_ref.wait();
                    let ph = phase_ref.load(Ordering::Acquire);
                    if ph == 3 {
                        return;
                    }
                    let start = start_ref.load(Ordering::Acquire);
                    let cur_p0 = f64::from_bits(p0_ref.load(Ordering::Acquire));

                    if ph == 1 && tid < w2 {
                        {
                            let guard = bufs_ref.read().expect("bufs read");
                            local_pol.params.copy_from_slice(&guard.params);
                        }
                        let blk_start = tid * blocks_per_w;
                        let blk_end = ((tid + 1) * blocks_per_w).min(n_blocks);
                        let b_base = blk_start * BLOCK;
                        let mut out = wpnl_ref[tid].write().expect("pnl write");
                        for blk in blk_start..blk_end {
                            let idx_start = blk * BLOCK;
                            let idx_end = ((blk + 1) * BLOCK).min(bsize);
                            for b in idx_start..idx_end {
                                let loc = b - b_base;
                                let path = &paths[(start + b) % n_total];
                                let mut prev_d = 0.0_f64;
                                let mut wealth = cur_p0;
                                let caches = &mut path_caches[loc];
                                let d_actuals = &mut path_d_actuals[loc];
                                for t in 0..n_steps {
                                    let tau =
                                        cfg.spec.t_mat * (n_steps - t) as f64 / n_steps as f64;
                                    cfg.spec.fill_features(path, t, n_steps, prev_d, &mut feat);
                                    let nn_d = local_pol.forward_step(&feat, &mut caches[t]);
                                    let d = if let Some((r_bs, sig_bs)) = cfg.spec.residual_bs {
                                        bs_delta(path.s[t], cfg.spec.k, r_bs, sig_bs, tau) + nn_d
                                    } else {
                                        nn_d
                                    };
                                    wealth += d * (path.s[t + 1] - path.s[t])
                                        - cfg.cost * (d - prev_d).abs() * path.s[t];
                                    prev_d = d;
                                    d_actuals[t] = d;
                                }
                                out[loc] = wealth - (path.s[n_steps] - cfg.spec.k).max(0.0);
                            }
                        }
                    } else if ph == 2 && tid < w2 {
                        let guard = bufs_ref.read().expect("bufs read");
                        let d_pnl_slice = &guard.d_pnl;
                        let blk_start = tid * blocks_per_w;
                        let blk_end = ((tid + 1) * blocks_per_w).min(n_blocks);
                        let b_base = blk_start * BLOCK;
                        let mut bgrads = wgrad_ref[tid].write().expect("grad write");
                        for blk in blk_start..blk_end {
                            let acc_v = &mut bgrads[blk - blk_start];
                            acc_v.fill(0.0);
                            let idx_start = blk * BLOCK;
                            let idx_end = ((blk + 1) * BLOCK).min(bsize);
                            for b in idx_start..idx_end {
                                let loc = b - b_base;
                                let path = &paths[(start + b) % n_total];
                                let caches = &path_caches[loc];
                                let d_actuals = &path_d_actuals[loc];
                                let dp = d_pnl_slice[b];
                                local_g.fill(0.0);
                                let mut d_prev_delta_next = 0.0_f64;
                                for t in (0..n_steps).rev() {
                                    let c = &caches[t];
                                    let d_cur = d_actuals[t];
                                    let d_prev = prev_idx.map_or(0.0, |idx| c.x[idx]);
                                    let sgn = if d_cur >= d_prev { 1.0 } else { -1.0 };
                                    let d_direct = dp
                                        * ((path.s[t + 1] - path.s[t])
                                            - cfg.cost * sgn * path.s[t]);
                                    local_pol.backward_step(
                                        c,
                                        d_direct + d_prev_delta_next,
                                        &mut local_g,
                                        &mut d_in_buf,
                                        &mut dz2,
                                        &mut dz1,
                                    );
                                    let d_from_prev = prev_idx.map_or(0.0, |idx| d_in_buf[idx]);
                                    let d_from_diff = diff_idx.map_or(0.0, |idx| d_in_buf[idx]);
                                    d_prev_delta_next = if prev_idx.is_some() || diff_idx.is_some()
                                    {
                                        d_from_prev + d_from_diff + dp * cfg.cost * sgn * path.s[t]
                                    } else {
                                        0.0
                                    };
                                }
                                for j in 0..n_params {
                                    acc_v[j] += local_g[j];
                                }
                            }
                        }
                    }

                    be_ref.wait();
                }
            });
        }

        let mut pnl = vec![0.0_f64; bsize];
        let mut grad = Vec::with_capacity(n_params + 1);
        let mut params_ext = Vec::with_capacity(n_params + 1);

        for it in 0..cfg.n_iters {
            if cfg.lr_decay_step > 0 && it > 0 && it % cfg.lr_decay_step == 0 {
                opt.lr *= cfg.lr_decay_gamma;
            }
            let start = (it * cfg.batch_size) % n_total;

            // Pass 1
            start_ref.store(start, Ordering::Release);
            p0_ref.store(p0.to_bits(), Ordering::Release);
            {
                let mut guard = bufs_ref.write().expect("bufs write p1");
                guard.params.copy_from_slice(&pol.params);
            }
            phase_ref.store(1, Ordering::Release);
            bs_ref.wait();
            be_ref.wait();

            let mut off = 0;
            for chunk_rw in wpnl_ref.iter().take(w2) {
                let c = chunk_rw.read().expect("collect pnl");
                pnl[off..off + c.len()].copy_from_slice(&c);
                off += c.len();
            }

            let (loss, d_pnl) = loss_and_d_pnl(&pnl, cfg.loss_type);
            let log_stride = (cfg.n_iters / 20).max(20);
            if it == 0 || (it + 1) % log_stride == 0 || it + 1 == cfg.n_iters {
                loss_history.push((it + 1, loss));
            }

            // Pass 2
            {
                let mut guard = bufs_ref.write().expect("bufs write p2");
                guard.d_pnl.copy_from_slice(&d_pnl);
            }
            phase_ref.store(2, Ordering::Release);
            bs_ref.wait();
            be_ref.wait();

            grad.clear();
            grad.resize(n_params, 0.0);
            for bg_rw in wgrad_ref.iter().take(w2) {
                let bgs = bg_rw.read().expect("collect grad");
                for blk_v in bgs.iter() {
                    for j in 0..n_params {
                        grad[j] += blk_v[j];
                    }
                }
            }

            let dp0 = if cfg.train_p0 {
                d_pnl.iter().sum::<f64>()
            } else {
                0.0
            };
            grad.push(dp0);
            params_ext.clear();
            params_ext.extend_from_slice(&pol.params);
            params_ext.push(p0);
            opt.step(&mut params_ext, &grad);
            pol.params.copy_from_slice(&params_ext[..n_params]);
            p0 = params_ext[n_params];

            if checkpoints.contains(&(it + 1)) {
                snaps.push((it + 1, pol.clone(), p0));
            }
        }

        // Shutdown workers
        phase_ref.store(3, Ordering::Release);
        bs_ref.wait();
    });

    (
        TrainResult {
            policy: pol,
            p0,
            loss_history,
        },
        snaps,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hedge::{generate_gbm_paths, GbmParams};
    use std::time::Instant;

    fn train_policy_legacy_20(
        paths: &[PathData],
        cfg: TrainConfig,
        n_workers: usize,
    ) -> TrainResult {
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
        for it in 0..cfg.n_iters {
            if cfg.lr_decay_step > 0 && it > 0 && it % cfg.lr_decay_step == 0 {
                opt.lr *= cfg.lr_decay_gamma;
            }
            let start = (it * cfg.batch_size) % n_total;
            let bsize = cfg.batch_size.min(n_total);
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
                        let d_direct =
                            dp * ((path.s[t + 1] - path.s[t]) - cfg.cost * sgn * path.s[t]);
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

    #[test]
    fn test_persistent_pool_bit_identical_20_iters() {
        let gp = GbmParams {
            s0: 1.0,
            k: 1.0,
            r: 0.0,
            sigma: 0.2,
            mu: 0.0,
            t_mat: 0.25,
            n_steps: 30,
        };
        let paths = generate_gbm_paths(&gp, 1024, 4, 20261122, 0x9901);
        let cfg = TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 20,
            batch_size: 1024,
            lr: 5e-3,
            lr_decay_step: 10,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: FeatureSpec::basic(1.0, 0.25, false, false, true),
            init_p0: 0.04,
            train_p0: true,
            init_output_bias: 0.5,
            zero_output: false,
            seed: 20261122,
        };
        let t0 = Instant::now();
        let leg = train_policy_legacy_20(&paths, cfg, 4);
        let dt_leg = t0.elapsed().as_secs_f64();

        let t1 = Instant::now();
        let pool = train_policy(&paths, cfg, 4);
        let dt_pool = t1.elapsed().as_secs_f64();

        assert_eq!(leg.p0.to_bits(), pool.p0.to_bits());
        assert_eq!(leg.policy.params.len(), pool.policy.params.len());
        for (a, b) in leg.policy.params.iter().zip(pool.policy.params.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        eprintln!(
            "20-iter benchmark (w=4): legacy={:.4}s, pool={:.4}s, speedup={:.2}x",
            dt_leg,
            dt_pool,
            dt_leg / dt_pool.max(1e-9)
        );
    }
}
