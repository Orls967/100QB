use day004g_market_making::glt::{
    eigendecompose_m, solve_v_jacobi, theta_from_v, GltParams, QuoteTable,
};
use day004g_market_making::hjb_fd::solve_nonlinear_hjb_theta;
use day004g_market_making::kolmogorov::solve_forward_kolmogorov;
use day004g_market_making::rng::MAIN_SEED;
use day004g_market_making::sim::{run_mc_sim, MarketSpec};
use day004g_market_making::strategies::Strategy;

#[test]
fn test_worker_invariance_and_determinism() {
    let p = GltParams::default_params();
    let table = QuoteTable::build(p, 500);
    let mkt = MarketSpec {
        t_horizon: p.t_horizon,
        sigma: p.sigma,
        mu: 0.0,
        a_intens: p.a_intens,
        k: p.k,
        gamma: p.gamma,
        q_max: p.q_max,
        eps_adverse: 0.0,
    };
    let s1 = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table },
        2000,
        500,
        1,
        MAIN_SEED,
        0x1111,
    );
    let s4 = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table },
        2000,
        500,
        4,
        MAIN_SEED,
        0x1111,
    );
    let s4b = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table },
        2000,
        500,
        4,
        MAIN_SEED,
        0x1111,
    );
    assert_eq!(s1.mean_wealth.to_bits(), s4.mean_wealth.to_bits());
    assert_eq!(s1.certainty_equiv.to_bits(), s4.certainty_equiv.to_bits());
    assert_eq!(s4.mean_wealth.to_bits(), s4b.mean_wealth.to_bits());
}

#[test]
fn test_small_anchors_part1_and_part2() {
    let p = GltParams::default_params();
    let eig = eigendecompose_m(&p);
    let v0 = solve_v_jacobi(&p, &eig, p.t_horizon);
    let th_exact = theta_from_v(&p, &v0);
    let exact_ce = th_exact[p.idx(0)];

    // Nonlinear HJB check
    let (th_fd, _, _) = solve_nonlinear_hjb_theta(&p, 5e-4, 0.0);
    assert!((th_fd[p.idx(0)] - exact_ce).abs() < 1e-2);

    // Forward Kolmogorov vs MC (small 10k paths within 4 SE)
    let table = QuoteTable::build(p, 1000);
    let kol = solve_forward_kolmogorov(&table, 1000);
    let mkt = MarketSpec {
        t_horizon: p.t_horizon,
        sigma: p.sigma,
        mu: 0.0,
        a_intens: p.a_intens,
        k: p.k,
        gamma: p.gamma,
        q_max: p.q_max,
        eps_adverse: 0.0,
    };
    let sim_res = run_mc_sim(
        mkt,
        &Strategy::GltOptimal { table: &table },
        10_000,
        1000,
        4,
        MAIN_SEED,
        0x2222,
    );
    assert!((sim_res.certainty_equiv - exact_ce).abs() < 4.0 * sim_res.se_ce.max(0.05));
    assert!((sim_res.mean_executions - kol.expected_executions).abs() < 1.5);
}
