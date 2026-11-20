use day003b_optimal_execution::ac::{
    ac_continuous_sup_error, ac_efficient_frontier, ac_error_decomposition, evaluate_ac_moments,
    simulate_ac_mc, simulate_ac_mc_plain, solve_ac_qp, solve_ac_sinh_discrete, AcParams,
};
use day003b_optimal_execution::dp_liquidity::{
    compute_adaptivity_value, simulate_dp_paths, solve_liquidity_dp, verify_dp_bruteforce,
    LiquidityDpParams,
};
use day003b_optimal_execution::linalg::{cholesky, jacobi_min_eigenvalue, solve_tridiagonal};
use day003b_optimal_execution::transient::{
    analyze_non_pd_kernel, analyze_oscillating_kernel, build_kernel_matrix,
    scan_oscillating_kernels, solve_ow_exponential_analytic, solve_transient_general,
    solve_transient_with_risk, OwParams,
};

#[test]
fn test_tridiagonal_solver() {
    let lower = vec![0.0, -1.0, -1.0, -1.0];
    let diag = vec![4.0, 4.0, 4.0, 4.0];
    let upper = vec![-1.0, -1.0, -1.0, 0.0];
    let rhs = vec![3.0, 2.0, 2.0, 3.0];
    let x = solve_tridiagonal(&lower, &diag, &upper, &rhs);
    for xi in x {
        assert!((xi - 1.0).abs() < 1e-12);
    }
}

#[test]
fn test_ac_qp_vs_sinh_and_limits() {
    let p = AcParams::default();
    let qp = solve_ac_qp(&p);
    let sinh = solve_ac_sinh_discrete(&p);
    let max_rel =
        qp.x.iter()
            .zip(sinh.x.iter())
            .map(|(a, b)| (a - b).abs() / p.x_total)
            .fold(0.0f64, f64::max);
    assert!(max_rel < 1e-9, "max_rel={max_rel}");

    // Limit lambda -> 0 (TWAP)
    let mut p_twap = p;
    p_twap.lambda = 1e-16;
    let s_twap = solve_ac_qp(&p_twap);
    let twap_x = p.x_total / (p.n_steps as f64);
    for xk in &s_twap.x {
        assert!((xk - twap_x).abs() / p.x_total < 1e-7);
    }

    // Limit lambda -> inf (hampir seketika di trade 1)
    let mut p_fast = p;
    p_fast.lambda = 1.0;
    let s_fast = solve_ac_qp(&p_fast);
    assert!(s_fast.x[0] / p.x_total > 0.99);

    // Limit kontinu N -> inf konvergen (orde 1 karena eta_tilde = eta - 0.5*gamma*dt)
    let err_50 = ac_continuous_sup_error(&p, 50);
    let err_200 = ac_continuous_sup_error(&p, 200);
    assert!(
        err_200 < err_50 * 0.35,
        "err_50={err_50}, err_200={err_200}"
    );
}

#[test]
fn test_ac_frontier_monotone_convex_and_mc() {
    let p = AcParams::default();
    let frontier = ac_efficient_frontier(&p, 20);
    for w in frontier.windows(2) {
        assert!(w[1].1.var_c <= w[0].1.var_c + 1e-6);
        assert!(w[1].1.expected_c >= w[0].1.expected_c - 1e-6);
    }

    let qp = solve_ac_qp(&p);
    let mc = simulate_ac_mc(&p, &qp.x, 40_000, 4, 20261121);
    assert!(mc.z_mean < 3.0, "z_mean={}", mc.z_mean);
    assert!(mc.z_var < 3.0, "z_var={}", mc.z_var);

    // Paritas evaluasi langsung
    let (ec, vc) = evaluate_ac_moments(&p, &qp.x);
    assert!((ec - qp.expected_c).abs() < 1e-9);
    assert!((vc - qp.var_c).abs() < 1e-9);
}

#[test]
fn test_transient_ow_and_kernels() {
    let ow = OwParams::default();
    let g = build_kernel_matrix(ow.n_steps, ow.t_horizon, |d| {
        ow.lambda_ow * (-ow.rho * d).exp()
    });
    let (x_gen, c_gen) = solve_transient_general(&g, ow.x_total).unwrap();
    let (x_ana, c_ana) = solve_ow_exponential_analytic(&ow);

    let max_err = x_gen
        .iter()
        .zip(x_ana.iter())
        .map(|(a, b)| (a - b).abs() / ow.x_total)
        .fold(0.0f64, f64::max);
    assert!(max_err < 1e-10, "max_err={max_err}");
    assert!((c_gen - c_ana).abs() / c_ana < 1e-10);

    // Power-law kernel: semua x* >= 0, cost > 0
    let g_pl = build_kernel_matrix(30, 1.0, |d| ow.lambda_ow / (1.0 + d / 0.1).sqrt());
    let (x_pl, c_pl) = solve_transient_general(&g_pl, ow.x_total).unwrap();
    assert!(x_pl.iter().all(|&v| v >= -1e-12));
    assert!(c_pl > 0.0);

    // Kernel osilasi: PD namun punya komponen negatif (manipulasi harga)
    let (is_pd, min_eig, _x_osc, min_x, _c_osc) =
        analyze_oscillating_kernel(20, 1.0, ow.lambda_ow, 0.2, 6.0, 0.001, ow.x_total);
    assert!(is_pd && min_eig > 0.0);
    assert!(min_x < 0.0, "min_x={min_x}");

    // Kernel non-PD: Cholesky gagal & ada arah penurunan negatif
    let (chol_failed, min_eig_bad, d_quad, cost_pert) =
        analyze_non_pd_kernel(20, 1.0, ow.lambda_ow, ow.x_total);
    assert!(chol_failed);
    assert!(min_eig_bad < 0.0);
    assert!(d_quad < 0.0);
    assert!(cost_pert < 0.0);

    // Gabungan: G diagonal (impact sementara murni) cocok dengan Almgren-Chriss
    let p_ac = AcParams::default();
    let dt = p_ac.t_horizon / (p_ac.n_steps as f64);
    let eta_tilde = p_ac.eta - 0.5 * p_ac.gamma * dt;
    let mut g_diag = vec![vec![0.0; p_ac.n_steps]; p_ac.n_steps];
    for (i, row) in g_diag.iter_mut().enumerate() {
        row[i] = 2.0 * eta_tilde / dt;
    }
    let (x_comb, _) = solve_transient_with_risk(&g_diag, &p_ac).unwrap();
    let qp = solve_ac_qp(&p_ac);
    let diff_comb = x_comb
        .iter()
        .zip(qp.x.iter())
        .map(|(a, b)| (a - b).abs() / p_ac.x_total)
        .fold(0.0f64, f64::max);
    assert!(diff_comb < 1e-12);
    let _ = (cholesky(&g_diag), jacobi_min_eigenvalue(&g_pl));
}

#[test]
fn test_dp_liquidity_bruteforce_sim_adaptivity_and_determinism() {
    // Verifikasi brute-force N=5
    let p5 = LiquidityDpParams {
        n_steps: 5,
        ..Default::default()
    };
    let rel_err = verify_dp_bruteforce(&p5, 401);
    assert!(rel_err < 1e-8, "rel_err={rel_err}");

    // Simulasi vs V0 & invariansi jumlah thread (1 vs 4)
    let p = LiquidityDpParams::default();
    let sol = solve_liquidity_dp(&p);
    let sim1 = simulate_dp_paths(&p, &sol, 20_000, 1, 20261121);
    let sim4 = simulate_dp_paths(&p, &sol, 20_000, 4, 20261121);
    assert!((sim1.mean_j - sim4.mean_j).abs() == 0.0);
    assert!(sim1.z_score < 3.0, "z_score={}", sim1.z_score);

    // Nilai adaptivitas: DP selalu <= strategi statis
    let rep = compute_adaptivity_value(&p);
    assert!(rep.j_dp <= rep.j_static_stationary + 1e-6);
    assert!(rep.j_dp <= rep.j_static_liquid + 1e-6);
    assert!(rep.j_dp <= rep.j_static_illiquid + 1e-6);
    assert!(rep.adaptivity_pct_vs_stationary > 0.0);
}

// Kriteria dinyatakan sebelum dijalankan: MC TANPA antitetik (jalur independen) -> z mean harus berukuran
// O(1): 1e-4 < z < 3 untuk jadwal optimal dan TWAP (z ~ 1e-11 berarti uji vakum).
#[test]
fn test_ac_plain_mc_mean_is_not_vacuous() {
    let p = AcParams::default();
    let qp = solve_ac_qp(&p);
    let twap = vec![p.x_total / (p.n_steps as f64); p.n_steps];
    for (sched, seed) in [(&qp.x, 20261131u64), (&twap, 20261132u64)] {
        let r = simulate_ac_mc_plain(&p, sched, 100_000, 4, seed);
        assert!(r.z_mean > 1e-4 && r.z_mean < 3.0, "z_mean={}", r.z_mean);
        assert!(r.z_var < 3.0, "z_var={}", r.z_var);
    }
}

// Dekomposisi galat limit kontinu: total = d1 + d2 secara identitas; sup total <= sup d1 + sup d2.
#[test]
fn test_ac_error_decomposition_identity() {
    let p = AcParams::default();
    for n in [10usize, 50, 200] {
        let (e1, e2, et, _opp) = ac_error_decomposition(&p, n);
        let tot = ac_continuous_sup_error(&p, n);
        assert!((et - tot).abs() < 1e-14, "n={n}: et={et} tot={tot}");
        assert!(et <= e1 + e2 + 1e-15);
    }
}

// Kelas "PD tetapi x* negatif" harus ada pada scan yang dinyatakan (grid kecil).
#[test]
fn test_oscillating_scan_finds_pd_manipulation_class() {
    let rows = scan_oscillating_kernels(
        20,
        1.0,
        5e-7,
        1_000_000.0,
        &[0.2, 0.5],
        &[6.0, 10.0],
        &[0.0, 0.05],
    );
    let found = rows.iter().filter(|r| r.is_pd && r.min_x < 0.0).count();
    assert!(found > 0);
}
