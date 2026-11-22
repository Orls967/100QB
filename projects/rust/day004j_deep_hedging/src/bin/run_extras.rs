use std::path::Path;
use std::time::Instant;

use day004j_deep_hedging::bs::{bs_call, bs_delta, whalley_wilmott_half_width};
use day004j_deep_hedging::eval::measure_nn_notrade_half_width;
use day004j_deep_hedging::hedge::{
    delta_l2_vs_bs, eval_strategy, generate_gbm_paths, FeatureSpec, GbmParams, Strategy,
};
use day004j_deep_hedging::jsonout::{obj, write_json_file};
use day004j_deep_hedging::nn::{Activation, LossType};
use day004j_deep_hedging::parallel::default_workers;
use day004j_deep_hedging::svg::{write_line_chart, Series};
use day004j_deep_hedging::train::{train_policy, train_policy_with_checkpoints, TrainConfig};

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
    let out_json = root.join("outputs/day004j/extras.json");
    let fig_dir = root.join("outputs/figures");

    let s0 = 1.0_f64;
    let k = 1.0_f64;
    let r = 0.0_f64;
    let t_mat = 0.25_f64;
    let sigma = 0.2_f64;
    let n_steps = 30_usize;
    let n_train = 16_384_usize;
    let r_eval = 100_000_usize;

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
    let db_mid = bs_delta(1.0, k, r, sigma, 0.125);

    let spec_bs_raw = FeatureSpec::basic(k, t_mat, false, false, true);
    let spec_bs_scaled = FeatureSpec {
        scale_inputs: true,
        ..spec_bs_raw
    };

    // 1. Batch=4096 convergence curve
    let (_res_b4096, snaps_4096) = train_policy_with_checkpoints(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 3000,
            batch_size: 4096,
            lr: 5e-3,
            lr_decay_step: 1000,
            lr_decay_gamma: 0.5,
            cost: 0.0,
            loss_type: LossType::Mse,
            spec: spec_bs_scaled,
            init_p0: bs_p0,
            train_p0: true,
            init_output_bias: 0.5,
            zero_output: false,
            seed: MAIN_SEED ^ 0x2150,
        },
        &[360, 1000, 3000],
        n_workers,
    );
    let mut curve_b4096 = Vec::new();
    for (it_c, pol_c, p0_c) in &snaps_4096 {
        let ev = eval_strategy(
            &gbm_eval,
            k,
            t_mat,
            0.0,
            *p0_c,
            10.0,
            Strategy::NeuralNet {
                policy: pol_c,
                spec: spec_bs_scaled,
            },
            n_workers,
        );
        let l2 = delta_l2_vs_bs(&gbm_eval, pol_c, spec_bs_scaled, r, sigma, n_workers);
        curve_b4096.push((*it_c, ev.sd_pnl, l2, *p0_c));
    }

    // 2. Zero-cost diagnostics & residual BS+NN
    let run_diag = |name: &str,
                    spec_d: FeatureSpec,
                    hidden: usize,
                    bias: f64,
                    lr: f64,
                    decay: usize,
                    zero_out: bool,
                    iters: usize,
                    tag: u64| {
        let tr = train_policy(
            &gbm_train,
            TrainConfig {
                hidden,
                act: Activation::Tanh,
                n_iters: iters,
                batch_size: 1024,
                lr,
                lr_decay_step: decay,
                lr_decay_gamma: 0.5,
                cost: 0.0,
                loss_type: LossType::Mse,
                spec: spec_d,
                init_p0: bs_p0,
                train_p0: true,
                init_output_bias: bias,
                zero_output: zero_out,
                seed: MAIN_SEED ^ tag,
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
                spec: spec_d,
            },
            n_workers,
        );
        let l2 = delta_l2_vs_bs(&gbm_eval, &tr.policy, spec_d, r, sigma, n_workers);
        (name.to_string(), ev.sd_pnl, l2, tr.p0)
    };
    let diag_experiments = [
        run_diag(
            "raw_inputs_1000",
            spec_bs_raw,
            32,
            0.0,
            3e-3,
            400,
            false,
            1000,
            0x2181,
        ),
        run_diag(
            "scaled_inputs_bias0_1000",
            spec_bs_scaled,
            32,
            0.0,
            3e-3,
            400,
            false,
            1000,
            0x2182,
        ),
        run_diag(
            "scaled_inputs_bias05_1000",
            spec_bs_scaled,
            32,
            0.5,
            3e-3,
            400,
            false,
            1000,
            0x2183,
        ),
        run_diag(
            "scaled_lr1e2_1000",
            spec_bs_scaled,
            32,
            0.5,
            1e-2,
            350,
            false,
            1000,
            0x2184,
        ),
        run_diag(
            "scaled_width64_1000",
            spec_bs_scaled,
            64,
            0.5,
            3e-3,
            400,
            false,
            1000,
            0x2185,
        ),
    ];
    let spec_residual = FeatureSpec {
        residual_bs: Some((r, sigma)),
        ..spec_bs_scaled
    };
    let diag_residual = run_diag(
        "residual_bs_plus_nn_1000",
        spec_residual,
        32,
        0.0,
        1e-3,
        400,
        true,
        1000,
        0x2186,
    );

    // 3. Activation & derivative feature test (A4 hypothesis: tanh rounds the band?)
    let hw_ww_001 = whalley_wilmott_half_width(1.0, k, r, sigma, 0.125, 0.001, 10.0);
    let spec_diff = FeatureSpec {
        delta_diff_bs: Some((r, sigma)),
        ..spec_bs_scaled
    };
    let mut act_tests = Vec::new();
    for (nm, act, sp, tag) in [
        ("tanh_10k", Activation::Tanh, spec_bs_scaled, 0x6101_u64),
        ("relu_10k", Activation::Relu, spec_bs_scaled, 0x6102),
        ("softplus_10k", Activation::Softplus, spec_bs_scaled, 0x6103),
        ("tanh_diff_10k", Activation::Tanh, spec_diff, 0x6104),
        ("relu_diff_10k", Activation::Relu, spec_diff, 0x6105),
    ] {
        let tr = train_policy(
            &gbm_train,
            TrainConfig {
                hidden: 32,
                act,
                n_iters: 10_000,
                batch_size: 1024,
                lr: 5e-3,
                lr_decay_step: 2500,
                lr_decay_gamma: 0.5,
                cost: 0.001,
                loss_type: LossType::Entropic { lambda: 10.0 },
                spec: sp,
                init_p0: bs_p0,
                train_p0: false,
                init_output_bias: 0.5,
                zero_output: false,
                seed: MAIN_SEED ^ tag,
            },
            n_workers,
        );
        let (hw, _, _) = measure_nn_notrade_half_width(&tr.policy, sp, 15, n_steps, 1.0, db_mid);
        let ev = eval_strategy(
            &gbm_eval,
            k,
            t_mat,
            0.001,
            bs_p0,
            10.0,
            Strategy::NeuralNet {
                policy: &tr.policy,
                spec: sp,
            },
            n_workers,
        );
        act_tests.push((
            nm,
            hw,
            hw / hw_ww_001,
            ev.entropic_risk,
            ev.indifference_price,
        ));
    }

    // 4. Cost misspecification & width4/noprev ablations
    let tr_full = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 2000,
            batch_size: 1024,
            lr: 5e-3,
            lr_decay_step: 700,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs_scaled,
            init_p0: bs_p0,
            train_p0: false,
            init_output_bias: 0.5,
            zero_output: false,
            seed: MAIN_SEED ^ 0x5001,
        },
        n_workers,
    );
    let (_hw_full, xs_band, full_ys) =
        measure_nn_notrade_half_width(&tr_full.policy, spec_bs_scaled, 15, n_steps, 1.0, db_mid);

    let tr_w4 = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 4,
            act: Activation::Tanh,
            n_iters: 2000,
            batch_size: 1024,
            lr: 5e-3,
            lr_decay_step: 700,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_bs_scaled,
            init_p0: bs_p0,
            train_p0: false,
            init_output_bias: 0.5,
            zero_output: false,
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
            spec: spec_bs_scaled,
        },
        n_workers,
    );
    let (hw_w4, _, w4_ys) =
        measure_nn_notrade_half_width(&tr_w4.policy, spec_bs_scaled, 15, n_steps, 1.0, db_mid);

    let spec_noprev = FeatureSpec {
        scale_inputs: true,
        ..FeatureSpec::basic(k, t_mat, false, false, false)
    };
    let tr_noprev = train_policy(
        &gbm_train,
        TrainConfig {
            hidden: 32,
            act: Activation::Tanh,
            n_iters: 2000,
            batch_size: 1024,
            lr: 5e-3,
            lr_decay_step: 700,
            lr_decay_gamma: 0.5,
            cost: 0.001,
            loss_type: LossType::Entropic { lambda: 10.0 },
            spec: spec_noprev,
            init_p0: bs_p0,
            train_p0: false,
            init_output_bias: 0.5,
            zero_output: false,
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

    write_line_chart(
        &fig_dir.join("day004j_ablation-width4-noprev.svg"),
        "Ablation at c=0.001, lam=10: Full NN vs Width=4 vs No prev_delta",
        "Previous Position delta_{t-1}",
        "New Position delta_t",
        &[
            Series {
                label: "Full NN (w=32)",
                x: &xs_band,
                y: &full_ys,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "Width 4",
                x: &xs_band,
                y: &w4_ys,
                color: "#f59e0b",
                dashed: false,
            },
            Series {
                label: "No prev_delta",
                x: &xs_band,
                y: &noprev_ys,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let total_secs = t_total.elapsed().as_secs_f64();
    let doc = obj([
        ("total_seconds", total_secs.into()),
        (
            "iter_curve_b4096",
            curve_b4096
                .iter()
                .map(|&(it, sd_v, l2_v, p0_v)| {
                    obj([
                        ("iters", it.into()),
                        ("sd_pnl", sd_v.into()),
                        ("l2_delta", l2_v.into()),
                        ("p0", p0_v.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
        (
            "diagnostics_ablation",
            diag_experiments
                .iter()
                .map(|(nm, sd_v, l2_v, p0_v)| {
                    obj([
                        ("name", nm.clone().into()),
                        ("sd_pnl", (*sd_v).into()),
                        ("l2_delta", (*l2_v).into()),
                        ("p0", (*p0_v).into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
        (
            "residual_diagnostic",
            obj([
                ("name", diag_residual.0.into()),
                ("sd_pnl", diag_residual.1.into()),
                ("l2_delta", diag_residual.2.into()),
                ("p0", diag_residual.3.into()),
            ]),
        ),
        (
            "activation_and_diff_feature",
            act_tests
                .iter()
                .map(|&(nm, hw, ratio, rho, ind)| {
                    obj([
                        ("name", nm.into()),
                        ("hw_nn", hw.into()),
                        ("hw_ww", hw_ww_001.into()),
                        ("ratio_vs_ww", ratio.into()),
                        ("entropic_risk", rho.into()),
                        ("indifference_price", ind.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
        ("ablation_width4_rho10", m_ab_w4.entropic_risk.into()),
        ("ablation_width4_hw", hw_w4.into()),
        ("ablation_noprev_rho10", m_ab_noprev.entropic_risk.into()),
        ("ablation_noprev_hw", hw_noprev.into()),
    ]);
    write_json_file(&out_json, &doc).expect("write extras.json");
    println!("day004j run_extras selesai dalam {:.2}s", total_secs);
}
