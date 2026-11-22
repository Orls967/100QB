use day004g_market_making::glt::{
    eigendecompose_m, solve_v_jacobi, theta_from_v, GltParams, QuoteTable,
};
use day004g_market_making::hjb_fd::{solve_nonlinear_hjb_theta, solve_nonlinear_hjb_theta_table};
use day004g_market_making::kolmogorov::solve_forward_kolmogorov;
use day004g_market_making::rng::MAIN_SEED;
use day004g_market_making::sim::{run_continuous_mc_sim, run_mc_sim, MarketSpec};
use day004g_market_making::stats;
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

#[test]
fn test_richardson_order_and_continuous_sim_and_eps0_reduction() {
    // 1. Richardson order helper check on synthetic f(dt) = 10.0 + 2.0 * dt^1
    let w2 = vec![10.1; 100];
    let w3 = vec![10.05; 100];
    let (ce_r, se_crn, se_indep) =
        stats::richardson_ce_and_se(&w2, 10.1, 0.01, &w3, 10.05, 0.005, 0.1, 1.0);
    assert!((ce_r - 10.0).abs() < 1e-12);
    assert!(se_crn <= se_indep + 1e-12);

    // 2. eps=0 FD table reduction matches QuoteTable
    let p = GltParams::default_params();
    let (_, fd_tab) = solve_nonlinear_hjb_theta_table(&p, 1e-3, 0.0);
    let glt_tab = QuoteTable::build(p, 1000);
    let (b1, a1) = fd_tab.interp_quotes(0.5, 5);
    let (b2, a2) = glt_tab.interp_quotes(0.5, 5);
    assert!((b1 - b2).abs() < 1e-4);
    assert!((a1 - a2).abs() < 1e-4);

    // 3. Continuous-time simulator determinism & agreement with Kolmogorov (small 5k paths)
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
    let c1 = run_continuous_mc_sim(mkt, &glt_tab, 5000, 1, MAIN_SEED, 0x3333);
    let c4 = run_continuous_mc_sim(mkt, &glt_tab, 5000, 4, MAIN_SEED, 0x3333);
    assert_eq!(c1.certainty_equiv.to_bits(), c4.certainty_equiv.to_bits());
    let kol = solve_forward_kolmogorov(&glt_tab, 1000);
    assert!((c4.var_inventory - kol.var_inventory).abs() <= 4.0 * c4.se_var_inventory.max(0.5));
}
