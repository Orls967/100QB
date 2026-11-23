# frozen_string_literal: true

require_relative 'test_helper'

class TestFmRp < Minitest::Test
  include FactorPortfolioLab

  def test_fama_macbeth_recovers_factor_premia
    world = FactorWorld.build_world(n_assets: 30, n_factors: 3, seed: 1001)
    sim = FactorWorld.simulate(world, t_len: 1200, seed: 1002)
    fm = FamaMacbeth.fit(sim[:returns], sim[:factors])
    assert_in_delta 0.0, fm[:lambda_hat][0], 0.005
    3.times do |k|
      assert_in_delta world[:mu_f][k], fm[:lambda_hat][k + 1], 0.005
    end
    assert_operator fm[:cs_r2], :>, 0.3
  end

  def test_fama_macbeth_no_intercept_mode
    world = FactorWorld.build_world(n_assets: 15, n_factors: 2, seed: 1003)
    sim = FactorWorld.simulate(world, t_len: 300, seed: 1004)
    fm = FamaMacbeth.fit(sim[:returns], sim[:factors], add_intercept: false, hac_lags: 3)
    assert_equal 2, fm[:lambda_hat].length
    assert_operator fm[:se_hac][0], :>, 0.0
  end

  def test_shanken_multiplier_exceeds_one_and_inflates_se
    world = FactorWorld.build_world(n_assets: 20, n_factors: 2, seed: 1005)
    scaled_mu = world[:mu_f].map { |m| m * 4.0 }
    w2 = world.merge(mu_f: scaled_mu, mu: FactorWorld.compute_expected_returns(world[:alpha], world[:betas], scaled_mu))
    sim = FactorWorld.simulate(w2, t_len: 120, seed: 1006)
    fm = FamaMacbeth.fit(sim[:returns], sim[:factors])
    sh = Shanken.correct(fm, sim[:factors])

    assert_operator sh[:c_shanken], :>, 0.05
    assert_operator sh[:se_shanken][0], :>, fm[:se_fm][0]
    assert_operator sh[:se_shanken][1], :>=, fm[:se_fm][1]
  end

  def test_fm_study_shanken_improves_coverage
    study = FmStudy.run(master_seed: MAIN_SEED, n_assets: 35, n_factors: 3, t_len: 80, n_reps: 200, mu_f_scale: 5.0)
    assert_operator study[:mean_c_shanken], :>, 0.10
    assert_operator study[:coverage_shanken_95][0], :>=, study[:coverage_fm_95][0]
    assert_operator study[:mean_se_shanken][0], :>, study[:mean_se_fm][0]
  end

  def test_risk_parity_newton_equal_risk_contributions
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 1101)
    sol = RiskParity.solve_newton(world[:sigma])
    assert sol[:converged]
    assert_in_delta 1.0, sol[:weights].sum, 1e-12
    assert sol[:weights].all?(&:positive?)
    assert_operator sol[:rc_max_dev], :<, 1e-10
  end

  def test_risk_parity_ccd_matches_newton
    world = FactorWorld.build_world(n_assets: 15, n_factors: 3, seed: 1102)
    sol_n = RiskParity.solve_newton(world[:sigma])
    sol_c = RiskParity.solve_ccd(world[:sigma])
    assert sol_c[:converged]
    assert_operator sol_c[:rc_max_dev], :<, 1e-9
    15.times do |i|
      assert_in_delta sol_n[:weights][i], sol_c[:weights][i], 1e-8
    end
  end

  def test_risk_parity_custom_budgets
    sigma = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
    budgets = [0.5, 0.3, 0.2]
    sol = RiskParity.solve_newton(sigma, budgets: budgets)
    assert sol[:converged]
    3.times do |i|
      assert_in_delta budgets[i], sol[:rel_risk_contributions][i], 1e-10
    end
  end

  def test_risk_parity_diagonal_closed_form
    sigma = [[0.01, 0.0, 0.0], [0.0, 0.04, 0.0], [0.0, 0.0, 0.16]]
    sol = RiskParity.solve_newton(sigma)
    inv_vols = [10.0, 5.0, 2.5]
    expected = inv_vols.map { |v| v / inv_vols.sum }
    3.times { |i| assert_in_delta expected[i], sol[:weights][i], 1e-10 }
  end

  def test_rp_study_runs_cleanly
    study = RpStudy.run(master_seed: MAIN_SEED, sizes: [8, 12])
    assert_equal 2, study[:comparisons].length
    study[:comparisons].each do |row|
      assert_operator row[:max_weight_diff], :<, 1e-8
      assert_operator row[:newton_rc_max_dev], :<, 1e-10
    end
  end

  def test_risk_parity_vol_beats_equal_weight_when_vols_differ
    sigma = [[0.0025, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.09]]
    rp = RiskParity.solve_newton(sigma)
    eq = RiskParity.summarize_solution(sigma, [1.0 / 3, 1.0 / 3, 1.0 / 3], 0, true)
    assert_operator rp[:portfolio_vol], :<, eq[:portfolio_vol]
  end

  def test_fama_macbeth_t_stat_consistency
    world = FactorWorld.build_world(n_assets: 10, n_factors: 2, seed: 1105)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 1106)
    fm = FamaMacbeth.fit(sim[:returns], sim[:factors])
    3.times do |j|
      assert_in_delta fm[:lambda_hat][j] / fm[:se_fm][j], fm[:t_fm][j], 1e-10
    end
  end

  def test_shanken_zero_lambda_gives_c_zero
    world = FactorWorld.build_world(n_assets: 10, n_factors: 2, seed: 1107)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 1108)
    fm = FamaMacbeth.fit(sim[:returns], sim[:factors])
    fm_zero = fm.merge(lambda_hat: [0.0, 0.0, 0.0])
    sh = Shanken.correct(fm_zero, sim[:factors])
    assert_in_delta 0.0, sh[:c_shanken], 1e-12
    3.times { |j| assert_in_delta fm[:se_fm][j], sh[:se_shanken][j], 1e-12 }
  end

  def test_risk_parity_hhi_bounded
    world = FactorWorld.build_world(n_assets: 10, n_factors: 2, seed: 1109)
    sol = RiskParity.solve_newton(world[:sigma])
    hhi = sol[:weights].sum { |w| w * w }
    assert_operator hhi, :>=, 0.1 - 1e-12
    assert_operator hhi, :<=, 1.0
  end

  def test_risk_parity_ccd_convergence_within_budget
    world = FactorWorld.build_world(n_assets: 20, n_factors: 3, seed: 1110)
    sol = RiskParity.solve_ccd(world[:sigma], tol: 1e-11, max_iter: 300)
    assert sol[:converged]
    assert_operator sol[:iterations], :<=, 300
  end
end
