use day004h_rough_bergomi::calib::{eval_smile_sse, nelder_mead, CalibCrnBank};
use day004h_rough_bergomi::rbergomi::{
    run_call_estimators_comparison, run_martingale_anchors, RBergomiParams,
};
use day004h_rough_bergomi::rng::MAIN_SEED;
use day004h_rough_bergomi::smile::compute_smile;
use day004h_rough_bergomi::volterra::{hybrid_bias_study, run_cholesky_reference_check};

#[test]
fn test_volterra_cholesky_and_worker_invariance() {
    let c1 = run_cholesky_reference_check(10, 1.0, 0.1, 4000, 1, MAIN_SEED, 0x101);
    let c4 = run_cholesky_reference_check(10, 1.0, 0.1, 4000, 4, MAIN_SEED, 0x101);
    assert_eq!(
        c1.var_w_terminal_emp.to_bits(),
        c4.var_w_terminal_emp.to_bits()
    );
    assert!(c1.max_chol_recon_err < 1e-13);
    assert!(c1.max_cov_z_score < 4.0);

    let h1 = hybrid_bias_study(0.1, 1.0, &[50], 4000, 1, MAIN_SEED);
    let h4 = hybrid_bias_study(0.1, 1.0, &[50], 4000, 4, MAIN_SEED);
    assert_eq!(
        h1[0].var_w_hybrid_mc.to_bits(),
        h4[0].var_w_hybrid_mc.to_bits()
    );
}

#[test]
fn test_rbergomi_pricing_smile_and_calib_smoke() {
    let p = RBergomiParams::default_params();
    let m = run_martingale_anchors(p, &[0.5], 50, 5000, 4, MAIN_SEED);
    assert!(m[0].z_exp_w < 3.5);
    assert!(m[0].z_s_t < 3.5);

    let cmp1 = run_call_estimators_comparison(p, &[0.5], &[0.0], 50, 4000, 1, MAIN_SEED);
    let cmp4 = run_call_estimators_comparison(p, &[0.5], &[0.0], 50, 4000, 4, MAIN_SEED);
    assert_eq!(
        cmp1.rows[0].cond_call.to_bits(),
        cmp4.rows[0].cond_call.to_bits()
    );
    assert!(cmp1.rows[0].vrf > 1.5);

    let sm = compute_smile(p, &[0.25], &[-0.1, 0.0, 0.1], 50, 4000, 4, MAIN_SEED, 0x99);
    assert!(
        sm[0].iv > sm[2].iv,
        "negative rho should produce downward skew"
    );

    // Small calibration smoke test on 2000 CRN paths
    let bank = CalibCrnBank::build(2000, 50, 4, MAIN_SEED, 0x88);
    let mkt = compute_smile(
        p,
        &[0.5, 1.0],
        &[-0.1, 0.0, 0.1],
        50,
        4000,
        4,
        MAIN_SEED ^ 1,
        0x77,
    );
    let (x, _, _) = nelder_mead(&[0.2, 1.5, -0.5], &[-0.05, 0.2, -0.1], 25, 1e-6, |v| {
        eval_smile_sse(
            &bank,
            p.xi0,
            p.s0,
            v[0],
            v[1],
            v[2],
            &mkt,
            &[0.5, 1.0],
            &[-0.1, 0.0, 0.1],
            4,
        )
    });
    assert!(x[0] > 0.02 && x[0] < 0.45);
}
