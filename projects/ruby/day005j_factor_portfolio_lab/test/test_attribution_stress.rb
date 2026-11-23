# frozen_string_literal: true

require_relative 'test_helper'

class TestAttributionStress < Minitest::Test
  include FactorPortfolioLab

  def test_euler_decomposition_exact_sum
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 2101)
    w = RiskParity.solve_newton(world[:sigma])[:weights]
    dec = Attribution.euler_decomposition(w, world[:sigma])
    assert_operator dec[:euler_error], :<, 1e-12
    assert_in_delta 1.0, dec[:rel_rc].sum, 1e-12
  end

  def test_factor_risk_decomposition_matches_total_var
    world = FactorWorld.build_world(n_assets: 10, n_factors: 3, seed: 2102)
    w = Array.new(10, 0.1)
    dec = Attribution.factor_risk_decomposition(w, world[:betas], world[:sigma_f], world[:idio_var])
    exact_var = LinAlg.quad_form(w, world[:sigma])
    assert_in_delta exact_var, dec[:total_var], 1e-12
    assert_in_delta 1.0, dec[:systematic_share] + dec[:idiosyncratic_share], 1e-12
  end

  def test_misspecified_k_understates_systematic_share
    study = Attribution.misspecified_k_study(master_seed: MAIN_SEED, n_assets: 15, t_len: 300)
    k1 = study[:k_models].find { |m| m[:k_fit] == 1 }[:systematic_share]
    k3 = study[:k_models].find { |m| m[:k_fit] == 3 }[:systematic_share]
    assert_operator k1, :<, k3
  end

  def test_stress_scenario_linear_and_quadratic
    world = FactorWorld.build_world(n_assets: 10, n_factors: 3, seed: 2103)
    w = Array.new(10, 0.1)
    df = [-0.08, 0.02, -0.03]
    gamma = [[-2.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]
    res = Stress.scenario_pnl(weights: w, betas: world[:betas], delta_f: df, gamma_f: gamma)
    assert_operator res[:quadratic_pnl], :<, 0.0
    assert_in_delta res[:linear_pnl] + res[:quadratic_pnl], res[:total_pnl], 1e-12
  end

  def test_reverse_stress_test_closed_form_optimality
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 2104)
    w = RiskParity.solve_newton(world[:sigma])[:weights]
    rev = Stress.reverse_stress_test(weights: w, betas: world[:betas], sigma_f: world[:sigma_f], target_loss: 0.12)
    assert_in_delta(-0.12, rev[:achieved_pnl], 1e-12)
    assert_in_delta rev[:theory_mahalanobis_sq], rev[:mahalanobis_sq], 1e-12

    # Perturb along nullspace of g' delta_f = 0 and verify Mahalanobis norm strictly increases
    g = LinAlg.mat_vec(LinAlg.transpose(world[:betas]), w)
    orth = [-g[1], g[0], 0.0]
    pert = rev[:shock_df].each_with_index.map { |v, i| v + (0.05 * orth[i]) }
    assert_in_delta(-0.12, LinAlg.dot(g, pert), 1e-12)
    pert_m = LinAlg.quad_form(pert, LinAlg.invert(world[:sigma_f]))
    assert_operator pert_m, :>, rev[:mahalanobis_sq]
  end

  def test_stress_study_runs_cleanly
    st = Stress.run_study(master_seed: MAIN_SEED, n_assets: 12)
    assert_equal 3, st[:scenarios].keys.length
    assert_in_delta(-0.15, st[:reverse_rp][:achieved_pnl], 1e-12)
  end

  def test_robustness_t4_increases_cvar
    rob = Robustness.run(master_seed: MAIN_SEED, n_grs_reps: 40, n_port_worlds: 12)
    assert_operator rob[:portfolio_t4]['EqualWeight'][:cvar_95], :>, rob[:portfolio_gaussian]['EqualWeight'][:cvar_95]
  end

  def test_complexity_scaling_positive_slopes
    cmp = Complexity.run(master_seed: MAIN_SEED, sizes: [10, 20, 40])
    assert_operator cmp[:slope_invert], :>, 1.0
    assert_operator cmp[:rows].length, :==, 3
  end

  def test_json_out_handles_nan_and_inf
    sanitized = JsonOut.sanitize({ a: 1.5, b: Float::NAN, c: [Float::INFINITY, 2.0] })
    assert_equal({ 'a' => 1.5, 'b' => nil, 'c' => [nil, 2.0] }, sanitized)
  end

  def test_svg_line_and_bar_write_valid_xml
    dir = File.join(Dir.tmpdir, "day005j_svg_test_#{Process.pid}")
    p1 = File.join(dir, 'line.svg')
    p2 = File.join(dir, 'bar.svg')
    Svg.line_chart(path: p1, title: 'T', x_label: 'X', y_label: 'Y', series: [{ name: 'A', x: [1, 2], y: [3, 4] }])
    Svg.bar_chart(path: p2, title: 'B', y_label: 'Y', labels: %w[A B], values: [1.0, 2.0])
    assert File.exist?(p1) && File.read(p1).include?('<svg')
    assert File.exist?(p2) && File.read(p2).include?('</svg>')
  ensure
    FileUtils.rm_rf(dir)
  end
end
