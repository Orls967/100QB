# frozen_string_literal: true

require_relative 'test_helper'

class TestTsGrs < Minitest::Test
  include FactorPortfolioLab

  def test_ols_exact_recovery_and_se_types
    rng = Rng::Sampler.new(404)
    x = Array.new(200) { [1.0, rng.next_normal, rng.next_normal] }
    y_clean = x.map { |r| 0.5 + (1.2 * r[1]) - (0.7 * r[2]) }
    fit_clean = Ols.fit(x, y_clean)
    assert_in_delta 0.5, fit_clean[:beta][0], 1e-10
    assert_in_delta 1.2, fit_clean[:beta][1], 1e-10
    assert_in_delta(-0.7, fit_clean[:beta][2], 1e-10)
    assert_in_delta 1.0, fit_clean[:r2], 1e-10

    y_noisy = y_clean.each_with_index.map { |v, t| v + ((1.0 + x[t][1].abs) * 0.2 * rng.next_normal) }
    fit_noisy = Ols.fit(x, y_noisy, hac_lags: 4)
    3.times do |j|
      assert_operator fit_noisy[:se_hc1][j], :>, 0.0
      assert_operator fit_noisy[:se_hac][j], :>, 0.0
    end
  end

  def test_ols_newey_west_inflates_for_ar1_errors
    rng = Rng::Sampler.new(505)
    x = Array.new(400) { [1.0, rng.next_normal] }
    e = 0.0
    y = x.map do |r|
      e = (0.75 * e) + (0.2 * rng.next_normal)
      1.0 + (0.5 * r[1]) + e
    end
    fit = Ols.fit(x, y, hac_lags: 6)
    assert_operator fit[:se_hac][0], :>, fit[:se_ols][0] * 1.4
  end

  def test_ts_regression_matches_single_asset_ols
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 606)
    sim = FactorWorld.simulate(world, t_len: 120, seed: 607)
    ts_fit = TimeSeriesRegression.fit(sim[:returns], sim[:factors])
    x_mat = sim[:factors].map { |f| [1.0] + f }

    5.times do |i|
      y_i = sim[:returns].map { |r| r[i] }
      ols_i = Ols.fit(x_mat, y_i)
      assert_in_delta ols_i[:beta][0], ts_fit[:alphas][i], 1e-10
      assert_in_delta ols_i[:beta][1], ts_fit[:betas][i][0], 1e-10
      assert_in_delta ols_i[:se_ols][0], ts_fit[:se_alpha][i], 1e-10
      assert_in_delta ols_i[:se_hc1][0], ts_fit[:se_alpha_hc1][i], 1e-10
      assert_in_delta ols_i[:r2], ts_fit[:r2][i], 1e-10
    end
  end

  def test_grs_dimensions_and_wald_relation
    world = FactorWorld.build_world(n_assets: 10, n_factors: 3, seed: 707)
    sim = FactorWorld.simulate(world, t_len: 120, seed: 708)
    res = Grs.test(sim[:returns], sim[:factors])
    assert_equal 10, res[:df1]
    assert_equal 107, res[:df2]
    assert_operator res[:wald_chi2], :>, res[:grs_f] * 10.0
  end

  def test_grs_raises_when_insufficient_dof
    world = FactorWorld.build_world(n_assets: 15, n_factors: 3, seed: 808)
    sim = FactorWorld.simulate(world, t_len: 18, seed: 809)
    assert_raises(ArgumentError) { Grs.test(sim[:returns], sim[:factors]) }
  end

  def test_grs_h0_calibration_vs_wald_overrejection
    regimes = [{ label: 'N25_T60_H0', n_assets: 25, t_len: 60, n_factors: 3, alpha_scale: 0.0 }]
    study = GrsStudy.run(master_seed: MAIN_SEED, n_reps: 250, regimes: regimes)
    r = study[:regimes][0]
    assert_in_delta 0.05, r[:grs_rejection_05], 0.04
    assert_operator r[:wald_rejection_05], :>, r[:grs_rejection_05] + 0.10
  end

  def test_grs_h1_power
    regimes = [{ label: 'N10_T120_H1', n_assets: 10, t_len: 120, n_factors: 3, alpha_scale: 0.006 }]
    study = GrsStudy.run(master_seed: MAIN_SEED, n_reps: 150, regimes: regimes)
    assert_operator study[:regimes][0][:grs_rejection_05], :>, 0.70
  end

  def test_two_regime_simulation_splits_volatility
    world = FactorWorld.build_world(n_assets: 6, n_factors: 2, seed: 909)
    sim = FactorWorld.simulate_two_regime(world, t_len: 400, seed: 910)
    first = sim[:returns][0...200].map { |r| r[0] }
    second = sim[:returns][200...400].map { |r| r[0] }
    var1 = first.sum { |x| x * x } / 200.0
    var2 = second.sum { |x| x * x } / 200.0
    assert_operator var2, :>, var1 * 2.5
  end
end
