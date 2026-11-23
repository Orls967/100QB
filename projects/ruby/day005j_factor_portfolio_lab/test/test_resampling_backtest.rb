# frozen_string_literal: true

require_relative 'test_helper'

class TestResamplingBacktest < Minitest::Test
  include FactorPortfolioLab

  def test_ledoit_wolf_shrinkage_delta_in_unit_interval
    world = FactorWorld.build_world(n_assets: 15, n_factors: 3, seed: 2001)
    sim = FactorWorld.simulate(world, t_len: 80, seed: 2002)
    lw = Shrinkage.shrink_cov_lw(sim[:returns])
    assert lw[:delta].between?(0.0, 1.0)
    assert_in_delta LinAlg.trace(lw[:sample_cov]), LinAlg.trace(lw[:sigma]), 1e-10
  end

  def test_james_stein_shrinks_toward_gmv_mean
    world = FactorWorld.build_world(n_assets: 10, n_factors: 2, seed: 2003)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 2004)
    mu_hat = Pca.column_means(sim[:returns], 10)
    sigma_hat = LinAlg.cov_matrix(sim[:returns])
    mu_js = Shrinkage.shrink_mu_js(mu_hat, sigma_hat, 100)
    spread_raw = mu_hat.max - mu_hat.min
    spread_js = mu_js.max - mu_js.min
    assert_operator spread_js, :<=, spread_raw + 1e-12
  end

  def test_michaud_reduces_concentration_vs_plugin
    res = Resampling.run_michaud(master_seed: MAIN_SEED, n_assets: 12, t_len: 90, n_boot: 60, kappa: 0.5)
    assert_operator res[:resampled][:hhi], :<, res[:plugin][:hhi]
    assert_in_delta 1.0, res[:resampled][:weights].sum, 1e-10
  end

  def test_all_8_strategies_produce_valid_weights
    world = FactorWorld.build_world(n_assets: 10, n_factors: 3, seed: 2005)
    sim = FactorWorld.simulate(world, t_len: 90, seed: 2006)
    Strategy.all_strategies.each do |strat|
      w = strat.weights(sim[:returns], factors_window: sim[:factors])
      assert_equal 10, w.length, "#{strat.name} length"
      assert_in_delta 1.0, w.sum, 1e-9, "#{strat.name} sum"
      assert w.all? { |wi| wi >= -1e-9 }, "#{strat.name} non-negative"
    end
  end

  def test_backtest_strictly_no_lookahead
    world = FactorWorld.build_world(n_assets: 8, n_factors: 2, seed: 2007)
    sim = FactorWorld.simulate(world, t_len: 100, seed: 2008)
    r_orig = sim[:returns].map(&:dup)
    r_mut = sim[:returns].map(&:dup)
    (70...100).each { |t| r_mut[t].map! { |v| v + 0.5 } }

    strat = Strategy::MinVarShrink.new
    bt1 = Backtest.run(returns: r_orig, strategy: strat, window: 60, rebalance_every: 5, tc_bps: 5.0)
    bt2 = Backtest.run(returns: r_mut, strategy: strat, window: 60, rebalance_every: 5, tc_bps: 5.0)

    # First two rebalance weight vectors (t=60, t=65) depend only on t < 70!
    assert_equal bt1[:weight_history][0], bt2[:weight_history][0]
    assert_equal bt1[:weight_history][1], bt2[:weight_history][1]
    assert_in_delta bt1[:port_returns][0], bt2[:port_returns][0], 1e-12
  end

  def test_backtest_transaction_costs_reduce_return
    world = FactorWorld.build_world(n_assets: 8, n_factors: 2, seed: 2009)
    sim = FactorWorld.simulate(world, t_len: 120, seed: 2010)
    strat = Strategy::TangencyShrink.new
    bt0 = Backtest.run(returns: sim[:returns], strategy: strat, window: 60, rebalance_every: 3, tc_bps: 0.0)
    bt20 = Backtest.run(returns: sim[:returns], strategy: strat, window: 60, rebalance_every: 3, tc_bps: 20.0)
    assert_operator bt20[:metrics][:ann_return], :<, bt0[:metrics][:ann_return]
  end

  def test_metrics_max_drawdown_and_cvar
    rets = [0.05, -0.10, -0.10, 0.04, 0.02, -0.05]
    mdd = Metrics.max_drawdown(rets)
    assert_in_delta 0.19, mdd, 1e-10
    var, cvar = Metrics.var_cvar(rets, alpha: 0.2)
    assert_operator cvar, :>=, var - 1e-12
  end

  def test_diversification_ratio_exceeds_one
    world = FactorWorld.build_world(n_assets: 10, n_factors: 2, seed: 2011)
    w = Array.new(10, 0.1)
    dr = Metrics.diversification_ratio(w, world[:sigma])
    assert_operator dr, :>, 1.0
  end

  def test_backtest_study_smoke
    study = BacktestStudy.run(master_seed: MAIN_SEED, n_worlds: 6, n_assets: 8, t_len: 100, window: 40)
    assert_equal 8, study[:summary].keys.length
    assert_operator study[:summary]['RiskParitySample'][5.0][:ann_vol], :>, 0.0
  end

  def test_shrink_cov_diag_preserves_diagonal
    cov = [[0.04, 0.02], [0.02, 0.09]]
    shrunk = Shrinkage.shrink_cov_diag(cov, delta: 0.5)
    assert_in_delta 0.04, shrunk[0][0], 1e-12
    assert_in_delta 0.01, shrunk[0][1], 1e-12
  end
end
