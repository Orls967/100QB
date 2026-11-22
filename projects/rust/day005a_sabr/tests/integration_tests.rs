use day005a_sabr::cev::{cev_absorb_prob, cev_call_exact, cev_call_pde_cn, cev_mc_exact};
use day005a_sabr::sabr_mc::{sabr_cond_rho0_mc, sabr_euler_mc, EulerScheme};

#[test]
fn test_worker_invariance_1_vs_8_bit_identical() {
    let ks = [0.02, 0.04, 0.06];
    let c1 = cev_mc_exact(0.04, 0.04, 0.0325, 0.5, 5.0, 2048, 1, 20261123);
    let c8 = cev_mc_exact(0.04, 0.04, 0.0325, 0.5, 5.0, 2048, 8, 20261123);
    assert_eq!(c1.0.to_bits(), c8.0.to_bits());
    assert_eq!(c1.1.to_bits(), c8.1.to_bits());
    assert_eq!(c1.3.to_bits(), c8.3.to_bits());

    let r1 = sabr_cond_rho0_mc(0.04, &ks, 1.0, 0.0325, 0.5, 0.4, 50, 1024, 1, 20261123);
    let r8 = sabr_cond_rho0_mc(0.04, &ks, 1.0, 0.0325, 0.5, 0.4, 50, 1024, 8, 20261123);
    for i in 0..ks.len() {
        assert_eq!(r1.prices[i].to_bits(), r8.prices[i].to_bits());
        assert_eq!(r1.ses[i].to_bits(), r8.ses[i].to_bits());
    }

    let e1 = sabr_euler_mc(
        0.04,
        &ks,
        1.0,
        0.0325,
        0.5,
        -0.3,
        0.4,
        50,
        1024,
        EulerScheme::Absorbing,
        1,
        20261123,
    );
    let e8 = sabr_euler_mc(
        0.04,
        &ks,
        1.0,
        0.0325,
        0.5,
        -0.3,
        0.4,
        50,
        1024,
        EulerScheme::Absorbing,
        8,
        20261123,
    );
    for i in 0..ks.len() {
        assert_eq!(e1.prices[i].to_bits(), e8.prices[i].to_bits());
        assert_eq!(e1.ses[i].to_bits(), e8.ses[i].to_bits());
    }
}

#[test]
fn test_cev_three_ways_and_absorption() {
    let f0 = 0.04;
    let alpha = 0.0325;
    let t = 10.0;
    let p_abs_ex = cev_absorb_prob(f0, alpha, 0.5, t);
    let p_abs_analytic = (-2.0 * f0 / (alpha * alpha * t)).exp();
    assert!((p_abs_ex - p_abs_analytic).abs() < 1e-14);

    let c_ex = cev_call_exact(f0, 0.04, alpha, 0.5, t);
    let c_pde = cev_call_pde_cn(f0, 0.04, alpha, 0.5, t, 1200, 800);
    assert!((c_ex - c_pde).abs() / c_ex < 5e-6);

    let (c_mc, se_mc, p_mc, _, _) = cev_mc_exact(f0, 0.04, alpha, 0.5, t, 100_000, 4, 20261123);
    assert!((c_mc - c_ex).abs() <= 3.0 * se_mc);
    assert!((p_mc - p_abs_ex).abs() < 0.005);
}
