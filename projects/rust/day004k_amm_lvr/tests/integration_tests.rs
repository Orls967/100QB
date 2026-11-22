use day004k_amm_lvr::blocks::{exact_block_expectations_gl, mc_block_annualized};
use day004k_amm_lvr::cpmm::{expected_pool_value, mc_pool_moments, variance_pool_value};
use day004k_amm_lvr::parallel::reduce_tuple4_fixed_blocks;
use day004k_amm_lvr::rng::{derive_seed, Rng};

#[test]
fn test_worker_invariance_1_vs_8_bit_identical() {
    let run_with = |w: usize| {
        reduce_tuple4_fixed_blocks(2048, w, |r| {
            let mut rg = Rng::new(derive_seed(20261122, 0x9001, r as u64));
            let z = rg.normal();
            let v = (-0.5 * 0.64 / 365.0 + 0.8 * (1.0 / 365.0_f64).sqrt() * z).exp();
            [v, v.sqrt(), v * v, 1.0]
        })
    };
    let a = run_with(1);
    let b = run_with(8);
    for i in 0..4 {
        assert_eq!(a[i].to_bits(), b[i].to_bits());
    }
}

#[test]
fn test_gbm_pool_moments_3se_and_block_gl() {
    let (m, se, var) = mc_pool_moments(1.0, 0.8, 30.0 / 365.0, 200_000, 4, 20261122);
    let exp_m = expected_pool_value(1.0, 0.8, 30.0 / 365.0);
    let exp_v = variance_pool_value(1.0, 0.8, 30.0 / 365.0);
    assert!((m - exp_m).abs() <= 3.0 * se, "m={m} exp={exp_m} se={se}");
    assert!((var - exp_v).abs() < 0.002, "var={var} exp_v={exp_v}");

    let dt = 60.0 / (365.0 * 86_400.0);
    let (gp, _, _) = exact_block_expectations_gl(0.8, dt, 0.0005, 0.5);
    let (mp, se_p, _, _, _, _) = mc_block_annualized(0.8, dt, 0.0005, 300_000, 4, 20261122);
    assert!(
        (gp / dt - mp).abs() <= 3.0 * se_p,
        "gl={} mc={} se={se_p}",
        gp / dt,
        mp
    );
}
