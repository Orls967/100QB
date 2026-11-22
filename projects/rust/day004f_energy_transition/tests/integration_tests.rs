use day004f_energy_transition::afns::{verify_transition_mc, AfnsParams};
use day004f_energy_transition::battery::BatterySpec;
use day004f_energy_transition::dp::{build_kernel, solve_stochastic_dp, PriceGrid};
use day004f_energy_transition::integration;
use day004f_energy_transition::panel_a;
use day004f_energy_transition::panel_b;
use day004f_energy_transition::panel_c;
use day004f_energy_transition::panel_d;
use day004f_energy_transition::panel_e;
use day004f_energy_transition::parallel;
use day004f_energy_transition::power::{stationary_distribution, PowerModelSpec};
use day004f_energy_transition::rng::{Rng, SEED_DAY004A, SEED_DAY004E_KERNEL};

#[test]
fn test_worker_invariance_1_vs_4_and_8() {
    // 1. parallel::run_chunks deterministic across 1 vs 4 vs 8 workers
    let r1 = parallel::run_chunks(200, 1, |i| (i * 17 + 3) % 97);
    let r4 = parallel::run_chunks(200, 4, |i| (i * 17 + 3) % 97);
    let r8 = parallel::run_chunks(200, 8, |i| (i * 17 + 3) % 97);
    assert_eq!(r1, r4);
    assert_eq!(r1, r8);

    // 2. AFNS MC transition invariance across 1 vs 4 workers
    let p = AfnsParams::default_params();
    let m1 = verify_transition_mc(&p, [0.05, -0.02, 0.01], 1.0, 20_000, 8, 1, SEED_DAY004A);
    let m4 = verify_transition_mc(&p, [0.05, -0.02, 0.01], 1.0, 20_000, 8, 4, SEED_DAY004A);
    for i in 0..3 {
        assert_eq!(m1.mc_mean[i].to_bits(), m4.mc_mean[i].to_bits());
    }

    // 3. DP kernel invariance across 1 vs 4 workers
    let spec = PowerModelSpec::default_spec();
    let g17 = PriceGrid::coarse_17();
    let k1 = build_kernel(
        &spec,
        g17.clone(),
        100,
        1,
        SEED_DAY004E_KERNEL,
        0x4E17,
        false,
    );
    let k4 = build_kernel(&spec, g17, 100, 4, SEED_DAY004E_KERNEL, 0x4E17, false);
    for (a, b) in k1.probs.iter().zip(k4.probs.iter()) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn test_determinism_same_seed_two_runs() {
    let mut r1 = Rng::new(SEED_DAY004A);
    let mut r2 = Rng::new(SEED_DAY004A);
    for _ in 0..500 {
        assert_eq!(r1.next_u64(), r2.next_u64());
        assert_eq!(r1.normal().to_bits(), r2.normal().to_bits());
    }
}

#[test]
fn test_all_panels_small_smoke() {
    let pa = panel_a::run_panel_a(4, true);
    assert!(pa.max_integral_err_ns < 1e-12);
    assert!(pa.max_bootstrap_df_err < 1e-12);
    assert!(pa.anchor_t5.log_lik_abs_diff < 1e-9);

    let pb = panel_b::run_panel_b(4, true);
    assert_eq!(pb.baseline.spec.num_reps, 4);
    assert_eq!(pb.longspan.len(), 4);
    assert_eq!(pb.infill.len(), 3);

    let pc = panel_c::run_panel_c(4, true);
    assert!(pc.algebraic_c_over_h_max_err < 1e-12);
    assert_eq!(pc.bootstrap_20y.total_refits, 18);

    let pd = panel_d::run_panel_d(4, true);
    assert_eq!(pd.total_inner_refits, 24);

    let pe = panel_e::run_panel_e(4, true);
    assert!(pe.stationary_residual < 1e-13);
    assert!(pe.res_dp33.soc_violation_count == 0);
    assert!(pe.res_dp33.pf_dom_violations == 0);

    let pf = integration::run_integration(
        &pe.res_dp33.daily_profits,
        &pb.baseline.replications,
        &pc.jackknife_20y.full_fits,
        &pd.outer_fits,
        true,
    );
    assert!(pf.duration_rel_err < 1e-6);
    assert!(pf.krd.sum_vs_par_closed_err < 1e-10);
    assert!(pf.hedge.residual_node1.abs() < 1e-9);
    assert!(pf.hedge.residual_node2.abs() < 1e-9);

    let spec = PowerModelSpec::default_spec();
    let pi = stationary_distribution(&spec.transition_matrix).unwrap();
    assert!((pi.iter().sum::<f64>() - 1.0).abs() < 1e-13);
    let bat = BatterySpec::default_spec();
    let pol17 = solve_stochastic_dp(&pe.kernel_17, bat);
    assert!(pol17.max_bellman_residual < 1e-11);
}
