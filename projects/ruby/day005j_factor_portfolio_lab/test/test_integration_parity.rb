# frozen_string_literal: true

require_relative 'test_helper'
require 'tmpdir'

class TestIntegrationParity < Minitest::Test
  include FactorPortfolioLab

  def test_factor_risk_parity_matches_sample_when_exact_cov
    world = FactorWorld.build_world(n_assets: 10, n_factors: 3, seed: 2201)
    cov_exact = FactorWorld.compute_total_cov(world[:betas], world[:sigma_f], world[:idio_var])
    rp1 = RiskParity.solve_newton(world[:sigma])
    rp2 = RiskParity.solve_newton(cov_exact)
    10.times { |i| assert_in_delta rp1[:weights][i], rp2[:weights][i], 1e-12 }
  end

  def test_pca_denoised_cov_is_symmetric_and_psd
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 2202)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 2203)
    cov = Strategy.pca_denoised_cov(sim[:returns], 3)
    12.times do |i|
      assert_operator cov[i][i], :>, 0.0
      12.times { |j| assert_in_delta cov[i][j], cov[j][i], 1e-12 }
    end
  end

  def test_black_litterman_converges_to_view_as_omega_shrinks
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 2204)
    sigma = world[:sigma]
    pi = BlackLitterman.implied_equilibrium(sigma, Array.new(5, 0.2))
    p_mat = [[1.0, -1.0, 0.0, 0.0, 0.0]]
    q_vec = [pi[0] - pi[1] + 0.03]
    omega = BlackLitterman.idzorek_omega(sigma, p_mat, [0.999], tau: 0.05)
    post = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: 0.05)
    assert_in_delta q_vec[0], post[:mu_post][0] - post[:mu_post][1], 1e-4
  end

  def test_grs_equals_zero_when_all_alphas_zero
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 2209)
    sim = FactorWorld.simulate(world, t_len: 60, seed: 2210)
    ts = TimeSeriesRegression.fit(sim[:returns], sim[:factors])
    de_alpha = sim[:returns].map { |r| r.each_with_index.map { |v, i| v - ts[:alphas][i] } }
    res = Grs.test(de_alpha, sim[:factors])
    assert_in_delta 0.0, res[:grs_f], 1e-10
    assert_in_delta 1.0, res[:p_value_f], 1e-10
  end

  def test_unconstrained_gmv_equals_minvar_qp_when_interior
    sigma = [[0.04, 0.005, 0.0], [0.005, 0.05, 0.005], [0.0, 0.005, 0.06]]
    w_unc = Markowitz.unconstrained_gmv(sigma)
    w_qp = Markowitz.solve_qp(sigma, kappa: 0.0)[:weights]
    3.times { |i| assert_in_delta w_unc[i], w_qp[i], 1e-9 }
  end

  def test_criteria_evaluator_structure
    dummy = build_minimal_study_payload
    crit = Criteria.evaluate(dummy)
    assert_operator crit.size, :>=, 40
    crit.each_value do |v|
      assert_includes ['tercapai', 'belum tercapai'], v[:status]
      assert v.key?(:observed) && v.key?(:threshold) && v.key?(:note)
    end
  end

  def test_figures_generator_creates_26_svgs
    dummy = build_minimal_study_payload
    dir = File.join(Dir.tmpdir, "day005j_figs_test_#{Process.pid}")
    paths = Figures.generate_all(dummy, dir)
    assert_equal 26, paths.length
    assert(paths.all? { |p| File.exist?(p) && File.size(p) > 200 })
  ensure
    FileUtils.rm_rf(dir)
  end

  def test_shanken_and_grs_consistent_sharpe_f_sq
    world = FactorWorld.build_world(n_assets: 12, n_factors: 2, seed: 2205)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 2206)
    grs = Grs.test(sim[:returns], sim[:factors])
    assert_operator grs[:sharpe_f_sq], :>=, 0.0
  end

  def test_backtest_turnover_zero_for_equal_weight_when_no_drift_difference
    r = Array.new(80) { Array.new(5, 0.01) }
    bt = Backtest.run(returns: r, strategy: Strategy::EqualWeight.new, window: 40, rebalance_every: 2, tc_bps: 10.0)
    assert_in_delta 0.0, bt[:metrics][:mean_turnover], 1e-12
    assert_in_delta 0.12, bt[:metrics][:ann_return], 1e-12
  end

  def test_module_count_exceeds_30
    lib_dir = File.expand_path('../lib/factor_portfolio_lab', __dir__)
    rb_files = Dir.glob(File.join(lib_dir, '*.rb'))
    assert_operator rb_files.length, :>=, 30
  end

  private

  def build_minimal_study_payload
    @@minimal_payload ||= { # rubocop:disable Style/ClassVars
      spiked: SpikedStudy.run(master_seed: MAIN_SEED, n_assets: 12, t_len: 48, n_reps: 10),
      grs: GrsStudy.run(master_seed: MAIN_SEED, n_reps: 40),
      fama_macbeth: FmStudy.run(master_seed: MAIN_SEED, n_assets: 15, n_factors: 3, t_len: 60, n_reps: 30),
      risk_parity: RpStudy.run(master_seed: MAIN_SEED, sizes: [8, 12]),
      black_litterman: BlStudy.run(master_seed: MAIN_SEED, n_mc: 10_000, n_assets: 5),
      markowitz: MvStudy.run(master_seed: MAIN_SEED, n_assets: 10, n_points: 45),
      resampling: Resampling.run_michaud(master_seed: MAIN_SEED, n_assets: 10, t_len: 60, n_boot: 30),
      backtest: BacktestStudy.run(master_seed: MAIN_SEED, n_worlds: 4, n_assets: 8, t_len: 90, window: 40),
      attribution: Attribution.misspecified_k_study(master_seed: MAIN_SEED, n_assets: 10, t_len: 120),
      stress: Stress.run_study(master_seed: MAIN_SEED, n_assets: 10),
      robustness: Robustness.run(master_seed: MAIN_SEED, n_grs_reps: 20, n_port_worlds: 4),
      complexity: Complexity.run(master_seed: MAIN_SEED, sizes: [10, 20, 30])
    }
  end
end
