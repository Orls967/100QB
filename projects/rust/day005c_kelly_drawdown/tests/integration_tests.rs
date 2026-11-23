use day005c_kelly_drawdown::estimation::run_estimation_loss_study;
use day005c_kelly_drawdown::gbm::GbmParams;
use day005c_kelly_drawdown::heavy_tail::run_heavy_tail_study;
use day005c_kelly_drawdown::hitting::{prob_hit_finite, prob_hit_infinite, run_hitting_grid};
use day005c_kelly_drawdown::kelly_binary::run_binary_study;

#[test]
fn integration_worker_invariance_and_determinism() {
    let gbm = GbmParams::default();
    let alt = GbmParams {
        mu: 0.11,
        sigma: 0.30,
        r: 0.02,
    };

    // Binary Kelly: 1 vs 4 workers must be bit-identical
    let b1 = run_binary_study(0.55, 1.0, 200, 2048, &[0.05, 0.10], 20261123, 1);
    let b4 = run_binary_study(0.55, 1.0, 200, 2048, &[0.05, 0.10], 20261123, 4);
    for (p1, p4) in b1.points.iter().zip(b4.points.iter()) {
        assert_eq!(p1.e_ln_w_mc.to_bits(), p4.e_ln_w_mc.to_bits());
        assert_eq!(p1.arith_mean_mc.to_bits(), p4.arith_mean_mc.to_bits());
        assert_eq!(
            p1.prob_touch_half_mc.to_bits(),
            p4.prob_touch_half_mc.to_bits()
        );
    }

    // Hitting grid: 1 vs 4 workers bit-identical
    let h1 = run_hitting_grid(
        &gbm,
        &[1.0 / 12.0],
        &[1.0, 5.0],
        &[0.5],
        &[1.25],
        4096,
        20261123,
        1,
    );
    let h4 = run_hitting_grid(
        &gbm,
        &[1.0 / 12.0],
        &[1.0, 5.0],
        &[0.5],
        &[1.25],
        4096,
        20261123,
        4,
    );
    for (r1, r4) in h1.iter().zip(h4.iter()) {
        assert_eq!(r1.naive_mc.to_bits(), r4.naive_mc.to_bits());
        assert_eq!(r1.bridge_mc.to_bits(), r4.bridge_mc.to_bits());
    }

    // Estimation loss & heavy tail: 1 vs 4 workers bit-identical
    let e1 = run_estimation_loss_study(&gbm, &alt, &[5.0], 4096, 20261123, 1);
    let e4 = run_estimation_loss_study(&gbm, &alt, &[5.0], 4096, 20261123, 4);
    assert_eq!(e1[0].loss_mc.to_bits(), e4[0].loss_mc.to_bits());

    let t1 = run_heavy_tail_study(&gbm, 4.0, &[0.5, 1.0], 2.0, 1024, 20261123, 1);
    let t4 = run_heavy_tail_study(&gbm, 4.0, &[0.5, 1.0], 2.0, 1024, 20261123, 4);
    for (x1, x4) in t1.iter().zip(t4.iter()) {
        assert_eq!(x1.g_t4_mc.to_bits(), x4.g_t4_mc.to_bits());
        assert_eq!(
            x1.prob_hit_half_t4_mc.to_bits(),
            x4.prob_hit_half_t4_mc.to_bits()
        );
    }
}

#[test]
fn integration_exact_anchors_and_3se() {
    let gbm = GbmParams::default();
    let f_star = gbm.kelly_fraction();
    let g_star = gbm.kelly_growth();
    let s_star = gbm.vol_log(f_star);

    // Infinite vs finite horizon at T = 5000
    let p_inf = prob_hit_infinite(g_star, s_star, 0.5);
    let p_fin = prob_hit_finite(g_star, s_star, 5000.0, 0.5);
    assert!((p_inf - p_fin).abs() < 1e-8);

    // Binary MC anchor within 3 SE
    let bin = run_binary_study(0.55, 1.0, 500, 20_000, &[0.10], 20261123, 2);
    let pt = &bin.points[0];
    let z = (pt.e_ln_w_mc - pt.e_ln_w_exact).abs() / pt.e_ln_w_se;
    assert!(z < 3.0, "binary E[ln W] z = {z}");
}

#[test]
fn integration_seed_determinism_and_sensitivity() {
    let gbm = GbmParams::default();
    let a1 = run_hitting_grid(
        &gbm,
        &[1.0 / 12.0],
        &[5.0],
        &[0.5],
        &[1.25],
        2048,
        20261123,
        2,
    );
    let a2 = run_hitting_grid(
        &gbm,
        &[1.0 / 12.0],
        &[5.0],
        &[0.5],
        &[1.25],
        2048,
        20261123,
        2,
    );
    let b1 = run_hitting_grid(
        &gbm,
        &[1.0 / 12.0],
        &[5.0],
        &[0.5],
        &[1.25],
        2048,
        99999999,
        2,
    );
    assert_eq!(a1[0].bridge_mc.to_bits(), a2[0].bridge_mc.to_bits());
    assert_ne!(a1[0].bridge_mc.to_bits(), b1[0].bridge_mc.to_bits());
}
