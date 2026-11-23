# frozen_string_literal: true

require 'minitest/autorun'
require 'backtest_overfitting'

class TestEngine < Minitest::Test
  def test_strategy_library_has_300_unique_strategies
    lib = BacktestOverfitting::Strategy.build_library
    assert_equal 300, lib.length
    assert_equal 300, lib.map(&:name).uniq.length
  end

  def test_fast_positions_match_causal_signal_at_method
    rng = BacktestOverfitting::Rng::Sampler.new(101)
    rets = BacktestOverfitting::DataGenerator.iid_gaussian(80, sigma: 0.015, rng: rng)
    strats = [
      BacktestOverfitting::Strategy::Momentum.new(lookback: 5),
      BacktestOverfitting::Strategy::SmaCrossover.new(fast: 4, slow: 15),
      BacktestOverfitting::Strategy::MeanReversionZScore.new(window: 10, threshold: 1.0)
    ]
    strats.each do |s|
      fast_pos = s.generate_positions(rets)
      slow_pos = Array.new(rets.length) { |t| s.signal_at(rets[0...t]) }
      assert_equal slow_pos, fast_pos, "Mismatch pada #{s.name}"
    end
  end

  def test_no_lookahead_across_all_strategy_families
    rng = BacktestOverfitting::Rng::Sampler.new(202)
    rets = BacktestOverfitting::DataGenerator.iid_gaussian(200, sigma: 0.01, rng: rng)
    lib = BacktestOverfitting::Strategy.build_library
    [lib[0], lib[99], lib[100], lib[150], lib[200], lib[299]].each do |strat|
      assert BacktestOverfitting::BacktestEngine.verify_no_lookahead(strat, rets, 120, rng),
             "Look-ahead terdeteksi pada #{strat.name}"
    end
  end

  def test_backtest_engine_equity_zero_cost_and_turnover_identities
    rng = BacktestOverfitting::Rng::Sampler.new(303)
    rets = BacktestOverfitting::DataGenerator.iid_gaussian(120, sigma: 0.012, rng: rng)
    pos = BacktestOverfitting::Strategy::Momentum.new(lookback: 3).generate_positions(rets)
    r0 = BacktestOverfitting::BacktestEngine.run(rets, pos, cost_bps: 0.0)
    r5 = BacktestOverfitting::BacktestEngine.run(rets, pos, cost_bps: 5.0)
    assert_operator BacktestOverfitting::BacktestEngine.verify_equity_identity(r0), :<, 1.0e-14
    assert_operator BacktestOverfitting::BacktestEngine.verify_equity_identity(r5), :<, 1.0e-14
    rets.each_index { |t| assert_in_delta pos[t] * rets[t], r0[:net_returns][t], 1.0e-15 }
    exp_to = pos.each_with_index.sum(0.0) { |p, t| (p - (t.zero? ? 0.0 : pos[t - 1])).abs }
    assert_in_delta exp_to, r5[:total_turnover], 1.0e-12
    assert_operator r5[:final_equity], :<, r0[:final_equity] if exp_to.positive?
  end

  def test_ar1_sign_strategy_exact_expectation_anchor
    rng = BacktestOverfitting::Rng::Sampler.new(404)
    rets = BacktestOverfitting::DataGenerator.ar1_gaussian(200_000, phi: 0.1, sigma: 0.01, rng: rng)
    pos = BacktestOverfitting::Strategy::Momentum.new(lookback: 1).generate_positions(rets)
    net = BacktestOverfitting::BacktestEngine.run(rets, pos, cost_bps: 0.0)[:net_returns]
    mean_ret = net[1..].sum(0.0) / (rets.length - 1)
    exact = 0.1 * 0.01 * Math.sqrt(2.0 / Math::PI)
    assert_in_delta exact, mean_ret, 3.5 * 0.01 / Math.sqrt(rets.length - 1)
  end

  def test_gbm_buy_and_hold_log_equity_anchor
    rng = BacktestOverfitting::Rng::Sampler.new(505)
    dt = 1.0 / 252.0
    pos = Array.new(252, 1.0)
    logs = Array.new(4_000) do
      rets = BacktestOverfitting::DataGenerator.gbm_simple_returns(252, mu: 0.08, sigma: 0.20, dt: dt, rng: rng)
      Math.log(BacktestOverfitting::BacktestEngine.run(rets, pos, cost_bps: 0.0)[:final_equity])
    end
    mean_log = logs.sum(0.0) / logs.length
    assert_in_delta 0.08 - (0.5 * 0.04), mean_log, 3.0 * 0.20 / Math.sqrt(logs.length)
  end

  def test_metrics_moments_and_max_drawdown
    sum = BacktestOverfitting::Metrics.summary([0.10, -0.20, 0.05])
    assert_in_delta 0.20, sum[:max_drawdown], 1.0e-12
    assert_equal 3, sum[:n_obs]
  end

  def test_two_state_regime_switches_volatility
    rng = BacktestOverfitting::Rng::Sampler.new(606)
    rets = BacktestOverfitting::DataGenerator.two_state_regime(
      2_000, p00: 0.98, p11: 0.95, sigma0: 0.005, sigma1: 0.03, rng: rng
    )
    sd = BacktestOverfitting::Metrics.summary(rets)[:std]
    assert_operator sd, :>, 0.006
    assert_operator sd, :<, 0.03
  end
end
