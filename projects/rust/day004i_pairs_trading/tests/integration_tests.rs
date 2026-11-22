use day004i_pairs_trading::fpt::{expected_fpt_exact, simulate_fpt_mc, solve_fpt_bvp};
use day004i_pairs_trading::ou::{run_ou_replication_study, OuParams};
use day004i_pairs_trading::rng::MAIN_SEED;
use day004i_pairs_trading::strategy::{
    optimal_threshold_direct, optimal_threshold_foc, verify_strategy_mc_200yr,
};

#[test]
fn test_ou_replication_worker_invariance_and_fpt_anchor() {
    let p = OuParams::default_params();
    let r1 = run_ou_replication_study(p, &[2.0], 400, 1, MAIN_SEED);
    let r4 = run_ou_replication_study(p, &[2.0], 400, 4, MAIN_SEED);
    assert_eq!(r1[0].bias_theta.to_bits(), r4[0].bias_theta.to_bits());

    let sd_stat = p.stationary_sd();
    let ex = expected_fpt_exact(-sd_stat, sd_stat, p.theta, p.sigma);
    let (bvp, _) = solve_fpt_bvp(-sd_stat, sd_stat, p.theta, p.sigma, 1500);
    assert!((ex - bvp).abs() / ex < 1e-3);

    let mc1 = simulate_fpt_mc(
        -sd_stat,
        sd_stat,
        p.theta,
        p.sigma,
        1.0 / 2520.0,
        8000,
        1,
        MAIN_SEED,
        0x201,
    );
    let mc4 = simulate_fpt_mc(
        -sd_stat,
        sd_stat,
        p.theta,
        p.sigma,
        1.0 / 2520.0,
        8000,
        4,
        MAIN_SEED,
        0x201,
    );
    assert_eq!(mc1.mean_tau.to_bits(), mc4.mean_tau.to_bits());
    assert!((mc4.mean_tau - ex).abs() <= 3.5 * mc4.se_tau);

    let (a_dir, _) = optimal_threshold_direct(p.theta, p.sigma, 0.002);
    let a_foc = optimal_threshold_foc(p.theta, p.sigma, 0.002);
    assert!((a_dir - a_foc).abs() < 1e-6);

    let chk = verify_strategy_mc_200yr(p.theta, p.mu, p.sigma, a_foc, 0.002, MAIN_SEED);
    assert!(chk.z_trades < 3.5, "z_trades={}", chk.z_trades);
    assert!(chk.z_profit < 3.5, "z_profit={}", chk.z_profit);
}
