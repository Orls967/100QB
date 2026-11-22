use day004h_rough_bergomi::rng::MAIN_SEED;
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
