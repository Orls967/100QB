# frozen_string_literal: true

require 'minitest/autorun'
require 'garch_family'

class TestRiskBacktest < Minitest::Test
  def test_kupiec_exact_calibration_and_christoffersen_clustering
    hits_ok = Array.new(1000) { |i| (i % 20).zero? }
    k_ok = GarchFamily::RiskBacktest.kupiec_test(hits_ok, 0.05)
    assert_in_delta 0.0, k_ok[:lr_uc], 1.0e-12
    assert_in_delta 1.0, k_ok[:p_uc], 1.0e-12

    hits_cluster = Array.new(1000) { |i| i < 50 }
    c_clust = GarchFamily::RiskBacktest.christoffersen_test(hits_cluster, 0.05)
    assert_operator c_clust[:lr_ind], :>, 100.0
    assert_operator c_clust[:p_ind], :<, 1.0e-12
  end

  def test_acerbi_szekely_z2_zero_expectation_under_calibrated_normal
    rng = GarchFamily::Rng::Sampler.new(1006)
    n = 10_000
    losses = Array.new(n) { -rng.next_normal }
    alpha = 0.05
    nd = GarchFamily::Innovations::Normal.new
    vars = Array.new(n, nd.var_quantile(alpha))
    ess = Array.new(n, nd.expected_shortfall(alpha))
    res = GarchFamily::RiskBacktest.acerbi_szekely_z2(losses, vars, ess, alpha)
    assert_in_delta 0.0, res[:z2], 0.08
    assert_in_delta 0.05, res[:hits].count(true).to_f / n, 0.01
  end

  def test_acerbi_szekely_z2_negative_when_es_underestimated
    rng = GarchFamily::Rng::Sampler.new(1007)
    n = 4000
    losses = Array.new(n) { -rng.next_std_t(5.0) }
    nd = GarchFamily::Innovations::Normal.new
    vars = Array.new(n, nd.var_quantile(0.01))
    ess = Array.new(n, nd.expected_shortfall(0.01))
    res = GarchFamily::RiskBacktest.acerbi_szekely_z2(losses, vars, ess, 0.01)
    assert_operator res[:z2], :<, -0.25
  end

  def test_ewma_and_garch_variance_forecasts_positive_and_stationary
    ew = GarchFamily::Ewma.new(lambda_decay: 0.94)
    g = GarchFamily::Garch11.new
    p = GarchFamily::Mle::TRUE_GARCH
    [1, 5, 20, 100].each do |h|
      assert_in_delta 2.0e-4, ew.forecast_variance({}, 2.0e-4, h), 1.0e-15
      fh = g.forecast_variance(p, 2.0e-4, h)
      assert_operator fh, :>, 0.0
      assert_operator fh, :<=, 2.0e-4
    end
  end

  def test_risk_backtest_rolling_oos_smoke_no_errors
    out = GarchFamily::RiskBacktest.rolling_oos_study(seed: 99, n_total: 420, win: 300, refit_every: 60)
    assert_equal 10, out.length
    out.each_value do |row|
      assert_equal 120, row[:n_obs]
      assert row[:acerbi_z2].finite?
    end
  end
end
