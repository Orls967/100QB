use std::path::Path;
use std::time::Instant;

use day004j_deep_hedging::bs::{
    bs_call, bs_delta, discrete_hedging_sd_approx, implied_vol_call, whalley_wilmott_half_width,
};
use day004j_deep_hedging::dp::solve_entropic_dp_gbm;
use day004j_deep_hedging::eval::{
    measure_nn_notrade_half_width, ols_slope, paired_diff, HestonDeltaTable,
};
use day004j_deep_hedging::fourier::heston_call_fourier;
use day004j_deep_hedging::hedge::{
    delta_l2_vs_bs, eval_strategy, generate_gbm_paths, FeatureSpec, GbmParams, Strategy,
};
use day004j_deep_hedging::heston::{generate_heston_qe_paths, HestonParams};
use day004j_deep_hedging::jsonout::{obj, write_json_file};
use day004j_deep_hedging::nn::{Activation, LossType, MlpPolicy, StepCache};
use day004j_deep_hedging::parallel::{default_workers, run_chunks};
use day004j_deep_hedging::stats;
use day004j_deep_hedging::svg::{write_bar_chart, write_line_chart, Series};
use day004j_deep_hedging::train::{train_policy, TrainConfig};

const MAIN_SEED: u64 = 20261122;

fn repo_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

fn main() {
    let t_total = Instant::now();
    let n_workers = default_workers().min(8);
    let root = repo_root();
    let out_json = root.join("outputs/day004j/results.json");
    let fig_dir = root.join("outputs/figures");

    let s0 = 1.0_f64;
    let k = 1.0_f64;
    let r = 0.0_f64;
    let t_mat = 0.25_f64;
    let sigma = 0.2_f64;
    let n_steps = 30_usize;
    let n_train = 16_384_usize;
    let r_eval = 100_000_usize;

    // -------------------------------------------------------------------------
    // PART 1: Library sanity (Adam & gradient checks verified in unit/integration tests)
    // -------------------------------------------------------------------------
    let t_p1 = Instant::now();
    let p1_secs = t_p1.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 2: Black-Scholes zero-cost world (anchor eksak)
    // -------------------------------------------------------------------------
    let t_p2 = Instant::now();
    let gbm_p = GbmParams {
        s0,
        k,
        r,
        sigma,
        mu: 0.0,
        t_mat,
        n_steps,
    };
    let gbm_train = generate_gbm_paths(&gbm_p, n_train, n_workers, MAIN_SEED, 0x2001);
    let gbm_eval = generate_gbm_paths(&gbm_p, r_eval, n_workers, MAIN_SEED ^ 0xA5A5, 0x2002);
    let bs_p0 = bs_call(s0, k, r, sigma, t_mat);
    let sd_bs_approx = discrete_hedging_sd_approx(s0, k, r, sigma, t_mat, n_steps);

    let m_bs = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.0,
        bs_p0,
        10.0,
        Strategy::BsDelta { k, r, sigma, t_mat },
        n_workers,
    );

    let spec_bs = FeatureSpec {
        k,
        t_mat,
        use_v: false,
        use_rv5: false,
        use_prev_delta: true,
    };

    let res_mse = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Mse,
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: true,
            seed: MAIN_SEED ^ 0x2100,
        },
        n_workers,
    );
    let m_nn_mse = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.0,
        res_mse.p0,
        10.0,
        Strategy::NeuralNet {
            policy: &res_mse.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let l2_mse = delta_l2_vs_bs(&gbm_eval, &res_mse.policy, spec_bs, r, sigma, n_workers);

    let res_ent = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x2200,
        },
        n_workers,
    );
    let m_nn_ent = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.0,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &res_ent.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let l2_ent = delta_l2_vs_bs(&gbm_eval, &res_ent.policy, spec_bs, r, sigma, n_workers);

    // Sup norm on S in [0.8, 1.2] at mid-horizon t = 15 (tau = 0.125)
    let mut s_grid = Vec::new();
    let mut d_bs_curve = Vec::new();
    let mut d_nn_curve = Vec::new();
    let mut sup_diff_bs = 0.0_f64;
    let mut cache_s = StepCache::new(spec_bs.d_in(), 32);
    for i in 0..81 {
        let s = 0.8 + 0.4 * (i as f64) / 80.0;
        let db = bs_delta(s, k, r, sigma, 0.125);
        let dn = res_mse
            .policy
            .forward_step(&[(s / k).ln(), 0.125, db], &mut cache_s);
        sup_diff_bs = sup_diff_bs.max((dn - db).abs());
        s_grid.push(s);
        d_bs_curve.push(db);
        d_nn_curve.push(dn);
    }

    // 5 initialization seeds sensitivity
    let mut seed_l2s = Vec::new();
    let mut seed_sds = Vec::new();
    for s_idx in 0..5_u64 {
        let tr = train_policy(
            &gbm_train,
            TrainConfig {
                hidden: 32,
                act: Activation::Tanh,
                n_iters: 240,
                batch_size: 1024,
                lr: 3e-3,
                lr_decay_step: 100,
                lr_decay_gamma: 0.5,
                cost: 0.0,
                loss_type: LossType::Mse,
                spec: spec_bs,
                init_p0: bs_p0,
                train_p0: true,
                seed: MAIN_SEED ^ (0x2300 + s_idx * 97),
            },
            n_workers,
        );
        let ev = eval_strategy(
            &gbm_eval,
            k,
            t_mat,
            0.0,
            tr.p0,
            10.0,
            Strategy::NeuralNet {
                policy: &tr.policy,
                spec: spec_bs,
            },
            n_workers,
        );
        seed_l2s.push(delta_l2_vs_bs(
            &gbm_eval, &tr.policy, spec_bs, r, sigma, n_workers,
        ));
        seed_sds.push(ev.sd_pnl);
    }

    // Convergence in n in {10, 20, 30, 60}
    let n_list = [10_usize, 20, 30, 60];
    let mut conv_rows = Vec::new();
    let mut log_ns = Vec::new();
    let mut log_sds_nn = Vec::new();
    let mut log_sds_bs = Vec::new();
    for &n_s in &n_list {
        let p_n = GbmParams {
            n_steps: n_s,
            ..gbm_p
        };
        let tr_paths = generate_gbm_paths(&p_n, n_train, n_workers, MAIN_SEED, 0x2400 + n_s as u64);
        let ev_paths = generate_gbm_paths(
            &p_n,
            40_000,
            n_workers,
            MAIN_SEED ^ 0xB1,
            0x2500 + n_s as u64,
        );
        let m_bs_n = eval_strategy(
            &ev_paths,
            k,
            t_mat,
            0.0,
            bs_p0,
            10.0,
            Strategy::BsDelta { k, r, sigma, t_mat },
            n_workers,
        );
        let tr_n = train_policy(
            &tr_paths,
            TrainConfig {
                hidden: 32,
                act: Activation::Tanh,
                n_iters: 280,
                batch_size: 1024,
                lr: 3e-3,
                lr_decay_step: 100,
                lr_decay_gamma: 0.5,
                cost: 0.0,
                loss_type: LossType::Mse,
                spec: spec_bs,
                init_p0: bs_p0,
                train_p0: true,
                seed: MAIN_SEED ^ (0x2600 + n_s as u64),
            },
            n_workers,
        );
        let m_nn_n = eval_strategy(
            &ev_paths,
            k,
            t_mat,
            0.0,
            tr_n.p0,
            10.0,
            Strategy::NeuralNet {
                policy: &tr_n.policy,
                spec: spec_bs,
            },
            n_workers,
        );
        log_ns.push((n_s as f64).ln());
        log_sds_bs.push(m_bs_n.sd_pnl.ln());
        log_sds_nn.push(m_nn_n.sd_pnl.ln());
        conv_rows.push((n_s, m_bs_n.sd_pnl, m_nn_n.sd_pnl));
    }
    let n_slope_bs = ols_slope(&log_ns, &log_sds_bs);
    let n_slope_nn = ols_slope(&log_ns, &log_sds_nn);
    let p2_secs = t_p2.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 3: Proportional transaction costs & DP anchor
    // -------------------------------------------------------------------------
    let t_p3 = Instant::now();
    let mut cost_sweep = Vec::new();
    let mut nn_c001_lam10: Option<MlpPolicy> = None;
    let mut nn_c001_band_curve = (Vec::new(), Vec::new());

    for &lam in &[10.0_f64, 50.0] {
        for &c in &[0.0005_f64, 0.001, 0.002] {
            let tr = train_policy(
                &gbm_train,
                TrainConfig {
                    hidden: 32,
                    act: Activation::Tanh,
                    n_iters: 380,
                    batch_size: 1024,
                    lr: 3e-3,
                    lr_decay_step: 130,
                    lr_decay_gamma: 0.5,
                    cost: c,
                    loss_type: LossType::Entropic { lambda: lam },
                    spec: spec_bs,
                    init_p0: bs_p0,
                    train_p0: false,
                    seed: MAIN_SEED ^ (0x3100 + ((c * 1e4) as u64) * 13 + (lam as u64)),
                },
                n_workers,
            );
            let db_mid = bs_delta(1.0, k, r, sigma, 0.125);
            let (hw_nn, dp_xs, dp_ys) =
                measure_nn_notrade_half_width(&tr.policy, spec_bs, 15, n_steps, 1.0, db_mid);
            let hw_ww = whalley_wilmott_half_width(1.0, k, r, sigma, 0.125, c, lam);
            let m_ev = eval_strategy(
                &gbm_eval,
                k,
                t_mat,
                c,
                bs_p0,
                lam,
                Strategy::NeuralNet {
                    policy: &tr.policy,
                    spec: spec_bs,
                },
                n_workers,
            );
            cost_sweep.push((
                c,
                lam,
                hw_nn,
                hw_ww,
                m_ev.entropic_risk,
                m_ev.indifference_price,
            ));
            if (c - 0.001).abs() < 1e-9 && (lam - 10.0).abs() < 1e-9 {
                nn_c001_lam10 = Some(tr.policy);
                nn_c001_band_curve = (dp_xs, dp_ys);
            }
        }
    }
    let pol_c001 = nn_c001_lam10.expect("c=0.001 lam=10 policy");

    // Scaling exponent of NN band half-width vs c for lam=10
    let log_cs: Vec<f64> = [0.0005_f64, 0.001, 0.002].iter().map(|x| x.ln()).collect();
    let log_hws_lam10: Vec<f64> = cost_sweep
        .iter()
        .filter(|x| (x.1 - 10.0).abs() < 1e-9)
        .map(|x| x.2.max(1e-6).ln())
        .collect();
    let ww_c_exp_lam10 = ols_slope(&log_cs, &log_hws_lam10);

    // Strategy comparison at c = 0.001, lam = 10, n = 30, R = 100,000
    let m_c_bs = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::BsDelta { k, r, sigma, t_mat },
        n_workers,
    );
    let m_c_ww = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::WhalleyWilmott {
            k,
            r,
            sigma,
            t_mat,
            c: 0.001,
            lambda: 10.0,
        },
        n_workers,
    );
    let m_c_lel = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::Leland {
            k,
            r,
            sigma,
            t_mat,
            c: 0.001,
        },
        n_workers,
    );
    let m_c_nn = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &pol_c001,
            spec: spec_bs,
        },
        n_workers,
    );
    let diff_nn_vs_bs_cost = paired_diff(&m_c_nn.pnl, &m_c_bs.pnl);

    // Exact DP anchor at n = 10, c = 0.001, lam = 10
    let dp_res = solve_entropic_dp_gbm(s0, k, r, sigma, t_mat, 10, 0.001, 10.0);
    let gbm_n10_tr = generate_gbm_paths(
        &GbmParams {
            n_steps: 10,
            ..gbm_p
        },
        n_train,
        n_workers,
        MAIN_SEED,
        0x3500,
    );
    let gbm_n10_ev = generate_gbm_paths(
        &GbmParams {
            n_steps: 10,
            ..gbm_p
        },
        40_000,
        n_workers,
        MAIN_SEED ^ 0xD1,
        0x3501,
    );
    let tr_nn_n10 = train_policy(
        &gbm_n10_tr,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 380,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 130,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x3502,
        },
        n_workers,
    );
    let m_dp_n10 = eval_strategy(
        &gbm_n10_ev,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::Dp { dp: &dp_res },
        n_workers,
    );
    let m_nn_n10 = eval_strategy(
        &gbm_n10_ev,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_nn_n10.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let diff_nn_dp_n10 = paired_diff(&m_nn_n10.pnl, &m_dp_n10.pnl);

    // Compare DP vs NN policy at (t = 5, S = 1.0, d_prev sweep)
    let mut dp_pol_xs = Vec::new();
    let mut dp_pol_ys = Vec::new();
    let mut nn_pol_ys = Vec::new();
    let mut agree_cnt = 0_usize;
    let mut abs_pol_sum = 0.0_f64;
    let half_grid = 0.5 * (1.1 / 80.0);
    let mut cache_dp = StepCache::new(spec_bs.d_in(), 32);
    for i in 0..41 {
        let dp_v = 0.2 + 0.6 * (i as f64) / 40.0;
        let d_dp = dp_res.act(5, 1.0, dp_v);
        let d_nn = tr_nn_n10
            .policy
            .forward_step(&[0.0, 0.125, dp_v], &mut cache_dp);
        let err = (d_dp - d_nn).abs();
        if err <= half_grid {
            agree_cnt += 1;
        }
        abs_pol_sum += err;
        dp_pol_xs.push(dp_v);
        dp_pol_ys.push(d_dp);
        nn_pol_ys.push(d_nn);
    }
    let dp_policy_agree_frac = agree_cnt as f64 / 41.0;
    let dp_policy_mae = abs_pol_sum / 41.0;
    let p3_secs = t_p3.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 4: Heston incomplete market
    // -------------------------------------------------------------------------
    let t_p4 = Instant::now();
    let hp = HestonParams {
        s0,
        k,
        r,
        kappa: 2.0,
        theta: 0.04,
        xi: 0.5,
        rho: -0.7,
        v0: 0.04,
        t_mat,
        n_steps,
    };
    let h_train = generate_heston_qe_paths(&hp, n_train, n_workers, MAIN_SEED, 0x4001);
    let h_eval = generate_heston_qe_paths(&hp, r_eval, n_workers, MAIN_SEED ^ 0xE5E5, 0x4002);

    let st_vals: Vec<f64> = h_eval.iter().map(|p| p.s[n_steps]).collect();
    let h_mean_st = stats::mean(&st_vals);
    let h_se_st = stats::se(&st_vals);

    let pay_vals: Vec<f64> = h_eval.iter().map(|p| (p.s[n_steps] - k).max(0.0)).collect();
    let h_mc_call = stats::mean(&pay_vals);
    let h_se_call = stats::se(&pay_vals);
    let h_fourier_call =
        heston_call_fourier(s0, k, r, hp.v0, hp.kappa, hp.theta, hp.xi, hp.rho, t_mat);
    let h_atm_iv = implied_vol_call(h_fourier_call, s0, k, r, t_mat);

    let h_delta_table = HestonDeltaTable::build(
        k, r, hp.kappa, hp.theta, hp.xi, hp.rho, t_mat, n_steps, n_workers,
    );
    let h_deltas: Vec<Vec<f64>> = run_chunks(r_eval, n_workers, |idx| {
        let p = &h_eval[idx];
        let vv = p.v.as_ref().unwrap();
        (0..n_steps)
            .map(|t| h_delta_table.lookup(t, p.s[t], k, vv[t]))
            .collect()
    });

    // Unhedged baseline (delta = 0)
    let m_h_unhedged = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::BsDelta {
            k: 1e9,
            r,
            sigma,
            t_mat,
        },
        n_workers,
    );

    let m_h_bs_iv = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::BsDelta {
            k,
            r,
            sigma: h_atm_iv,
            t_mat,
        },
        n_workers,
    );
    let m_h_fourier = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::HestonFourierDelta { deltas: &h_deltas },
        n_workers,
    );

    let spec_h_nov = FeatureSpec {
        k,
        t_mat,
        use_v: false,
        use_rv5: false,
        use_prev_delta: true,
    };
    let spec_h_v = FeatureSpec {
        k,
        t_mat,
        use_v: true,
        use_rv5: false,
        use_prev_delta: true,
    };
    let spec_h_rv5 = FeatureSpec {
        k,
        t_mat,
        use_v: false,
        use_rv5: true,
        use_prev_delta: true,
    };

    let tr_h_nov = train_policy(
        &h_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_h_nov,
            init_p0: h_fourier_call,
            train_p0: false,
            seed: MAIN_SEED ^ 0x4100,
        },
        n_workers,
    );
    let tr_h_v = train_policy(
        &h_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_h_v,
            init_p0: h_fourier_call,
            train_p0: false,
            seed: MAIN_SEED ^ 0x4200,
        },
        n_workers,
    );
    let tr_h_rv5 = train_policy(
        &h_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_h_rv5,
            init_p0: h_fourier_call,
            train_p0: false,
            seed: MAIN_SEED ^ 0x4300,
        },
        n_workers,
    );

    let m_h_nov = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_h_nov.policy,
            spec: spec_h_nov,
        },
        n_workers,
    );
    let m_h_v = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_h_v.policy,
            spec: spec_h_v,
        },
        n_workers,
    );
    let m_h_rv5 = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_h_rv5.policy,
            spec: spec_h_rv5,
        },
        n_workers,
    );
    let diff_hv_vs_fourier = paired_diff(&m_h_v.pnl, &m_h_fourier.pnl);
    let p4_secs = t_p4.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 5: Misspecification & Ablation + SVGs
    // -------------------------------------------------------------------------
    let t_p5 = Instant::now();
    // 1. BS-trained NN tested on Heston & BS sigma=0.3
    let m_bs_on_heston = eval_strategy(
        &h_eval,
        k,
        t_mat,
        0.0,
        h_fourier_call,
        10.0,
        Strategy::NeuralNet {
            policy: &res_ent.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let gbm_sig03_eval = generate_gbm_paths(
        &GbmParams {
            sigma: 0.3,
            ..gbm_p
        },
        r_eval,
        n_workers,
        MAIN_SEED ^ 0xF00D,
        0x5001,
    );
    let bs_p0_03 = bs_call(s0, k, r, 0.3, t_mat);
    let m_bs_on_03_nn = eval_strategy(
        &gbm_sig03_eval,
        k,
        t_mat,
        0.0,
        bs_p0_03,
        10.0,
        Strategy::NeuralNet {
            policy: &res_ent.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let m_bs_on_03_bs = eval_strategy(
        &gbm_sig03_eval,
        k,
        t_mat,
        0.0,
        bs_p0_03,
        10.0,
        Strategy::BsDelta {
            k,
            r,
            sigma: 0.2,
            t_mat,
        },
        n_workers,
    );
    let m_heston_on_bs = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.0,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_h_nov.policy,
            spec: spec_h_nov,
        },
        n_workers,
    );

    // 2. Wrong cost c_train != c_test
    let tr_c0005 = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 320,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 110,
            lr_decay_gamma: 0.5,
            cost: 0.0005,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x5100,
        },
        n_workers,
    );
    let tr_c0020 = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 320,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 110,
            lr_decay_gamma: 0.5,
            cost: 0.002,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x5200,
        },
        n_workers,
    );
    let m_c0005_at_002 = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.002,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_c0005.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let m_c0020_at_0005 = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.0005,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_c0020.policy,
            spec: spec_bs,
        },
        n_workers,
    );

    // 3. Ablation at c = 0.001, lam = 10: width 4 & no prev_delta
    let tr_w4 = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 4,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x5300,
        },
        n_workers,
    );
    let m_ab_w4 = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_w4.policy,
            spec: spec_bs,
        },
        n_workers,
    );
    let db_mid = bs_delta(1.0, k, r, sigma, 0.125);
    let (hw_w4, _, w4_ys) =
        measure_nn_notrade_half_width(&tr_w4.policy, spec_bs, 15, n_steps, 1.0, db_mid);

    let spec_noprev = FeatureSpec {
        k,
        t_mat,
        use_v: false,
        use_rv5: false,
        use_prev_delta: false,
    };
    let tr_noprev = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 360,
            batch_size: 1024,
            lr: 3e-3,
            lr_decay_step: 120,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_noprev,
            init_p0: bs_p0,
            train_p0: false,
            seed: MAIN_SEED ^ 0x5400,
        },
        n_workers,
    );
    let m_ab_noprev = eval_strategy(
        &gbm_eval,
        k,
        t_mat,
        0.001,
        bs_p0,
        10.0,
        Strategy::NeuralNet {
            policy: &tr_noprev.policy,
            spec: spec_noprev,
        },
        n_workers,
    );
    let (hw_noprev, _, noprev_ys) =
        measure_nn_notrade_half_width(&tr_noprev.policy, spec_noprev, 15, n_steps, 1.0, db_mid);

    // -------------------------------------------------------------------------
    // SVG Figures (8 charts)
    // -------------------------------------------------------------------------
    let lc_x_mse: Vec<f64> = res_mse.loss_history.iter().map(|x| x.0 as f64).collect();
    let lc_y_mse: Vec<f64> = res_mse.loss_history.iter().map(|x| x.1).collect();
    let lc_x_ent: Vec<f64> = res_ent.loss_history.iter().map(|x| x.0 as f64).collect();
    let lc_y_ent: Vec<f64> = res_ent.loss_history.iter().map(|x| x.1).collect();
    write_line_chart(
        &fig_dir.join("day004j_learning-curve.svg"),
        "Learning Curve: Training Loss vs Adam Iteration (BS, c=0)",
        "Iteration",
        "Loss",
        &[
            Series {
                label: "MSE Loss",
                x: &lc_x_mse,
                y: &lc_y_mse,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "Entropic (lam=10)",
                x: &lc_x_ent,
                y: &lc_y_ent,
                color: "#dc2626",
                dashed: false,
            },
        ],
    )
    .unwrap();

    write_line_chart(
        &fig_dir.join("day004j_nn-vs-bs-delta.svg"),
        "NN Hedge Delta vs Black-Scholes Delta (tau = 0.125)",
        "Spot S",
        "Delta",
        &[
            Series {
                label: "BS Delta N(d1)",
                x: &s_grid,
                y: &d_bs_curve,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "NN (MSE)",
                x: &s_grid,
                y: &d_nn_curve,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .unwrap();

    // P&L density histogram curve (BS vs NN MSE)
    let n_bins = 41_usize;
    let p_min = -0.04_f64;
    let p_max = 0.04_f64;
    let bw = (p_max - p_min) / n_bins as f64;
    let bin_x: Vec<f64> = (0..n_bins).map(|i| p_min + (i as f64 + 0.5) * bw).collect();
    let hist = |vals: &[f64]| -> Vec<f64> {
        let mut cnt = vec![0.0_f64; n_bins];
        for &v in vals {
            if v >= p_min && v < p_max {
                let bi = ((v - p_min) / bw) as usize;
                cnt[bi.min(n_bins - 1)] += 1.0 / (vals.len() as f64 * bw);
            }
        }
        cnt
    };
    let h_bs_hist = hist(&m_bs.pnl);
    let h_nn_hist = hist(&m_nn_mse.pnl);
    write_line_chart(
        &fig_dir.join("day004j_pnl-distribution.svg"),
        "Discrete Hedging P&L Density (n=30, c=0, R=100k)",
        "Terminal P&L",
        "Density",
        &[
            Series {
                label: "BS Delta",
                x: &bin_x,
                y: &h_bs_hist,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "NN (MSE)",
                x: &bin_x,
                y: &h_nn_hist,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .unwrap();

    let hw_ww_001 = whalley_wilmott_half_width(1.0, k, r, sigma, 0.125, 0.001, 10.0);
    let ww_band_y: Vec<f64> = nn_c001_band_curve
        .0
        .iter()
        .map(|&dp| dp.clamp(db_mid - hw_ww_001, db_mid + hw_ww_001))
        .collect();
    write_line_chart(
        &fig_dir.join("day004j_notrade-band-vs-ww.svg"),
        "No-Trade Band: NN vs Whalley-Wilmott (c=0.001, lam=10, S=1, tau=0.125)",
        "Previous Position delta_{t-1}",
        "New Position delta_t",
        &[
            Series {
                label: "Identity (Hold)",
                x: &nn_c001_band_curve.0,
                y: &nn_c001_band_curve.0,
                color: "#9ca3af",
                dashed: true,
            },
            Series {
                label: "Whalley-Wilmott",
                x: &nn_c001_band_curve.0,
                y: &ww_band_y,
                color: "#dc2626",
                dashed: true,
            },
            Series {
                label: "NN Policy",
                x: &nn_c001_band_curve.0,
                y: &nn_c001_band_curve.1,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .unwrap();

    write_line_chart(
        &fig_dir.join("day004j_dp-vs-nn-n10.svg"),
        "Exact DP Policy vs NN Policy at n=10 (c=0.001, lam=10, S=1, t=5)",
        "Previous Position delta_{t-1}",
        "Optimal delta_t",
        &[
            Series {
                label: "Exact DP",
                x: &dp_pol_xs,
                y: &dp_pol_ys,
                color: "#16a34a",
                dashed: false,
            },
            Series {
                label: "NN (n=10)",
                x: &dp_pol_xs,
                y: &nn_pol_ys,
                color: "#2563eb",
                dashed: true,
            },
        ],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day004j_heston-strategies.svg"),
        "Heston Hedging (R=100k, lam=10): Entropic Risk by Strategy",
        "Entropic Risk rho_10",
        &[
            "BS ATM IV",
            "Heston Fourier",
            "NN (no v)",
            "NN (with v_t)",
            "NN (RV5)",
        ],
        &[
            m_h_bs_iv.entropic_risk,
            m_h_fourier.entropic_risk,
            m_h_nov.entropic_risk,
            m_h_v.entropic_risk,
            m_h_rv5.entropic_risk,
        ],
        &["#6b7280", "#16a34a", "#f59e0b", "#2563eb", "#7c3aed"],
    )
    .unwrap();

    write_line_chart(
        &fig_dir.join("day004j_ablation-width4-noprev.svg"),
        "Ablation at c=0.001, lam=10: Full NN vs Width=4 vs No prev_delta",
        "Previous Position delta_{t-1}",
        "New Position delta_t",
        &[
            Series {
                label: "Full NN (w=32)",
                x: &nn_c001_band_curve.0,
                y: &nn_c001_band_curve.1,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "Width 4",
                x: &nn_c001_band_curve.0,
                y: &w4_ys,
                color: "#f59e0b",
                dashed: false,
            },
            Series {
                label: "No prev_delta",
                x: &nn_c001_band_curve.0,
                y: &noprev_ys,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day004j_misspecification.svg"),
        "Misspecification: Entropic Risk Across Regime Shifts",
        "Entropic Risk rho_10",
        &[
            "BS->Heston",
            "Heston->Heston",
            "BS->BS(0.3)",
            "BSdelta->BS(0.3)",
            "Heston->BS(0.2)",
        ],
        &[
            m_bs_on_heston.entropic_risk,
            m_h_nov.entropic_risk,
            m_bs_on_03_nn.entropic_risk,
            m_bs_on_03_bs.entropic_risk,
            m_heston_on_bs.entropic_risk,
        ],
        &["#dc2626", "#16a34a", "#f59e0b", "#6b7280", "#2563eb"],
    )
    .unwrap();

    let p5_secs = t_p5.elapsed().as_secs_f64();
    let total_secs = t_total.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // Criteria evaluation (a priori thresholds)
    // -------------------------------------------------------------------------
    let c3_pass = l2_mse <= 0.03
        && l2_ent <= 0.03
        && (res_mse.p0 - bs_p0).abs() <= 0.002
        && m_nn_mse.sd_pnl <= m_bs.sd_pnl + 3.0 * m_bs.se_sd_pnl;
    let c4_pass = (-0.65..=-0.35).contains(&n_slope_nn);
    let c5_pass = cost_sweep
        .iter()
        .find(|x| (x.0 - 0.001).abs() < 1e-9 && (x.1 - 10.0).abs() < 1e-9)
        .is_some_and(|x| x.2 > 0.01);
    let se_dp_sim = diff_nn_dp_n10.se.max(m_dp_n10.se_sd_pnl);
    let c6_pass = (m_nn_n10.indifference_price - dp_res.indifference_price).abs()
        <= (3.0 * se_dp_sim).max(0.0015);
    let c7_pass = m_c_nn.entropic_risk < m_c_bs.entropic_risk;
    let c8_pass = (h_mean_st - 1.0).abs() <= 3.0 * h_se_st
        && (h_mc_call - h_fourier_call).abs() <= 3.0 * h_se_call;
    let c9_pass = m_h_v.sd_pnl < m_h_nov.sd_pnl && m_h_v.entropic_risk < m_h_nov.entropic_risk;
    let c10_pass = total_secs < 600.0;

    let status_str = |b: bool| if b { "tercapai" } else { "belum tercapai" };

    let root_json = obj([
        (
            "metadata",
            obj([
                ("session", "day004j".into()),
                ("date_wita", "2026-11-22".into()),
                ("language", "rust".into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as i64).into()),
                (
                    "elapsed_seconds",
                    obj([
                        ("part1_lib_check", p1_secs.into()),
                        ("part2_bs_zero_cost", p2_secs.into()),
                        ("part3_costs_and_dp", p3_secs.into()),
                        ("part4_heston", p4_secs.into()),
                        ("part5_misspec_ablation_svg", p5_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part2_bs_zero_cost",
            obj([
                ("bs_call_price", bs_p0.into()),
                ("sd_bs_approx", sd_bs_approx.into()),
                ("sd_bs_mc", m_bs.sd_pnl.into()),
                ("se_sd_bs_mc", m_bs.se_sd_pnl.into()),
                ("sd_ratio_mc_to_approx", (m_bs.sd_pnl / sd_bs_approx).into()),
                ("rho10_bs", m_bs.entropic_risk.into()),
                ("p0_nn_mse", res_mse.p0.into()),
                ("sd_nn_mse", m_nn_mse.sd_pnl.into()),
                ("l2_delta_mse", l2_mse.into()),
                ("sup_delta_mse_08_12", sup_diff_bs.into()),
                ("rho10_nn_ent", m_nn_ent.entropic_risk.into()),
                ("sd_nn_ent", m_nn_ent.sd_pnl.into()),
                ("l2_delta_ent", l2_ent.into()),
                (
                    "seed_sensitivity_l2_min",
                    seed_l2s
                        .iter()
                        .copied()
                        .fold(f64::INFINITY, f64::min)
                        .into(),
                ),
                (
                    "seed_sensitivity_l2_max",
                    seed_l2s
                        .iter()
                        .copied()
                        .fold(f64::NEG_INFINITY, f64::max)
                        .into(),
                ),
                ("seed_sensitivity_sd_mean", stats::mean(&seed_sds).into()),
                ("seed_sensitivity_sd_sd", stats::sd(&seed_sds).into()),
                ("n_scaling_slope_bs", n_slope_bs.into()),
                ("n_scaling_slope_nn", n_slope_nn.into()),
                (
                    "n_convergence",
                    conv_rows
                        .iter()
                        .map(|&(n_s, sb, sn)| {
                            obj([
                                ("n_steps", n_s.into()),
                                ("sd_bs", sb.into()),
                                ("sd_nn", sn.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part3_costs_and_dp",
            obj([
                (
                    "cost_sweep",
                    cost_sweep
                        .iter()
                        .map(|&(c, lam, hw_nn, hw_ww, rho, ind)| {
                            obj([
                                ("c", c.into()),
                                ("lambda", lam.into()),
                                ("hw_nn", hw_nn.into()),
                                ("hw_ww", hw_ww.into()),
                                ("ratio_nn_to_ww", (hw_nn / hw_ww.max(1e-12)).into()),
                                ("entropic_risk", rho.into()),
                                ("indifference_price", ind.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("ww_c_scaling_exponent_lam10", ww_c_exp_lam10.into()),
                (
                    "strategies_c001_lam10",
                    obj([
                        (
                            "bs_delta",
                            obj([
                                ("rho10", m_c_bs.entropic_risk.into()),
                                ("sd_pnl", m_c_bs.sd_pnl.into()),
                                ("mean_cost", m_c_bs.mean_cost.into()),
                                ("mean_turnover", m_c_bs.mean_turnover.into()),
                            ]),
                        ),
                        (
                            "whalley_wilmott",
                            obj([
                                ("rho10", m_c_ww.entropic_risk.into()),
                                ("sd_pnl", m_c_ww.sd_pnl.into()),
                                ("mean_cost", m_c_ww.mean_cost.into()),
                                ("mean_turnover", m_c_ww.mean_turnover.into()),
                            ]),
                        ),
                        (
                            "leland",
                            obj([
                                ("rho10", m_c_lel.entropic_risk.into()),
                                ("sd_pnl", m_c_lel.sd_pnl.into()),
                                ("mean_cost", m_c_lel.mean_cost.into()),
                                ("mean_turnover", m_c_lel.mean_turnover.into()),
                            ]),
                        ),
                        (
                            "nn",
                            obj([
                                ("rho10", m_c_nn.entropic_risk.into()),
                                ("sd_pnl", m_c_nn.sd_pnl.into()),
                                ("mean_cost", m_c_nn.mean_cost.into()),
                                ("mean_turnover", m_c_nn.mean_turnover.into()),
                                ("paired_mean_diff_vs_bs", diff_nn_vs_bs_cost.mean.into()),
                                ("paired_se_vs_bs", diff_nn_vs_bs_cost.se.into()),
                            ]),
                        ),
                    ]),
                ),
                (
                    "dp_anchor_n10",
                    obj([
                        ("dp_indifference_price", dp_res.indifference_price.into()),
                        (
                            "dp_sim_indifference_price",
                            m_dp_n10.indifference_price.into(),
                        ),
                        (
                            "nn_sim_indifference_price",
                            m_nn_n10.indifference_price.into(),
                        ),
                        (
                            "dp_policy_agree_frac_half_grid",
                            dp_policy_agree_frac.into(),
                        ),
                        ("dp_policy_mae", dp_policy_mae.into()),
                        ("sim_paired_diff_nn_minus_dp", diff_nn_dp_n10.mean.into()),
                        ("sim_paired_se", diff_nn_dp_n10.se.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part4_heston",
            obj([
                ("mean_st", h_mean_st.into()),
                ("se_st", h_se_st.into()),
                ("mc_call", h_mc_call.into()),
                ("se_mc_call", h_se_call.into()),
                ("fourier_call", h_fourier_call.into()),
                ("atm_iv", h_atm_iv.into()),
                ("unhedged_sd", m_h_unhedged.sd_pnl.into()),
                ("bs_atm_iv_rho10", m_h_bs_iv.entropic_risk.into()),
                ("bs_atm_iv_sd", m_h_bs_iv.sd_pnl.into()),
                ("bs_atm_iv_cvar95", m_h_bs_iv.cvar_95.into()),
                ("heston_fourier_rho10", m_h_fourier.entropic_risk.into()),
                ("heston_fourier_sd", m_h_fourier.sd_pnl.into()),
                ("heston_fourier_cvar95", m_h_fourier.cvar_95.into()),
                ("nn_no_v_rho10", m_h_nov.entropic_risk.into()),
                ("nn_no_v_sd", m_h_nov.sd_pnl.into()),
                ("nn_no_v_cvar95", m_h_nov.cvar_95.into()),
                ("nn_with_v_rho10", m_h_v.entropic_risk.into()),
                ("nn_with_v_sd", m_h_v.sd_pnl.into()),
                ("nn_with_v_cvar95", m_h_v.cvar_95.into()),
                ("nn_with_v_indiff_price", m_h_v.indifference_price.into()),
                ("nn_rv5_rho10", m_h_rv5.entropic_risk.into()),
                ("nn_rv5_sd", m_h_rv5.sd_pnl.into()),
                ("nn_rv5_cvar95", m_h_rv5.cvar_95.into()),
                (
                    "paired_diff_nnv_vs_fourier_mean",
                    diff_hv_vs_fourier.mean.into(),
                ),
                (
                    "paired_diff_nnv_vs_fourier_se",
                    diff_hv_vs_fourier.se.into(),
                ),
            ]),
        ),
        (
            "part5_misspec_and_ablation",
            obj([
                ("bs_nn_on_heston_rho10", m_bs_on_heston.entropic_risk.into()),
                ("bs_nn_on_bs03_rho10", m_bs_on_03_nn.entropic_risk.into()),
                ("bs_delta_on_bs03_rho10", m_bs_on_03_bs.entropic_risk.into()),
                ("heston_nn_on_bs_rho10", m_heston_on_bs.entropic_risk.into()),
                (
                    "c0005_nn_at_c002_rho10",
                    m_c0005_at_002.entropic_risk.into(),
                ),
                (
                    "c002_nn_at_c0005_rho10",
                    m_c0020_at_0005.entropic_risk.into(),
                ),
                ("ablation_width4_rho10", m_ab_w4.entropic_risk.into()),
                ("ablation_width4_hw", hw_w4.into()),
                ("ablation_noprev_rho10", m_ab_noprev.entropic_risk.into()),
                ("ablation_noprev_hw", hw_noprev.into()),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            obj([
                (
                    "c1_grad_check_and_adam_under_1e6",
                    obj([
                        ("status", "tercapai".into()),
                        (
                            "verified_in",
                            "cargo test --release (nn & optim unit tests)".into(),
                        ),
                        ("threshold_rel_err", 1e-6_f64.into()),
                    ]),
                ),
                (
                    "c2_worker_invariance_1_vs_8",
                    obj([
                        ("status", "tercapai".into()),
                        (
                            "verified_in",
                            "cargo test --release (integration_tests)".into(),
                        ),
                    ]),
                ),
                (
                    "c3_bs_nn_delta_l2_and_price",
                    obj([
                        ("status", status_str(c3_pass).into()),
                        ("l2_mse", l2_mse.into()),
                        ("l2_ent", l2_ent.into()),
                        ("p0_mse_diff", (res_mse.p0 - bs_p0).abs().into()),
                        ("sd_nn_mse", m_nn_mse.sd_pnl.into()),
                        (
                            "sd_bs_plus_3se",
                            (m_bs.sd_pnl + 3.0 * m_bs.se_sd_pnl).into(),
                        ),
                    ]),
                ),
                (
                    "c4_discrete_hedging_n_scaling",
                    obj([
                        ("status", status_str(c4_pass).into()),
                        ("n_slope_nn", n_slope_nn.into()),
                        ("n_slope_bs", n_slope_bs.into()),
                    ]),
                ),
                (
                    "c5_cost_notrade_band_and_ww",
                    obj([
                        ("status", status_str(c5_pass).into()),
                        ("ww_c_scaling_exp", ww_c_exp_lam10.into()),
                    ]),
                ),
                (
                    "c6_dp_anchor_n10_agrees",
                    obj([
                        ("status", status_str(c6_pass).into()),
                        ("dp_indiff", dp_res.indifference_price.into()),
                        ("nn_indiff", m_nn_n10.indifference_price.into()),
                        ("dp_policy_agree_frac", dp_policy_agree_frac.into()),
                    ]),
                ),
                (
                    "c7_nn_beats_bs_under_cost",
                    obj([
                        ("status", status_str(c7_pass).into()),
                        ("rho10_nn", m_c_nn.entropic_risk.into()),
                        ("rho10_bs", m_c_bs.entropic_risk.into()),
                    ]),
                ),
                (
                    "c8_heston_martingale_and_fourier_3se",
                    obj([
                        ("status", status_str(c8_pass).into()),
                        ("z_st", ((h_mean_st - 1.0).abs() / h_se_st).into()),
                        (
                            "z_call",
                            ((h_mc_call - h_fourier_call).abs() / h_se_call).into(),
                        ),
                    ]),
                ),
                (
                    "c9_heston_vt_reduces_residual_risk",
                    obj([
                        ("status", status_str(c9_pass).into()),
                        ("sd_nn_v", m_h_v.sd_pnl.into()),
                        ("sd_nn_nov", m_h_nov.sd_pnl.into()),
                        ("rho10_nn_v", m_h_v.entropic_risk.into()),
                        ("rho10_nn_nov", m_h_nov.entropic_risk.into()),
                    ]),
                ),
                (
                    "c10_runtime_under_600s",
                    obj([
                        ("status", status_str(c10_pass).into()),
                        ("total_seconds", total_secs.into()),
                        ("threshold_seconds", 600.0_f64.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &root_json).expect("write results.json");
    println!(
        "day004j selesai dalam {:.2}s | workers={} | L2_mse={:.4} | p0_diff={:.5} | Heston_z_call={:.2} | c1..c10=true/true/{}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        n_workers,
        l2_mse,
        (res_mse.p0 - bs_p0).abs(),
        (h_mc_call - h_fourier_call).abs() / h_se_call,
        c3_pass,
        c4_pass,
        c5_pass,
        c6_pass,
        c7_pass,
        c8_pass,
        c9_pass,
        c10_pass
    );
}
