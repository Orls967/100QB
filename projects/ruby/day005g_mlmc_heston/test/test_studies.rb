# frozen_string_literal: true

require 'minitest/autorun'
require 'mlmc_heston'

class TestStudies < Minitest::Test
  GBM_BASE = { spot: 100.0, rate: 0.05, sigma: 0.20, strike: 100.0, maturity: 1.0 }.freeze
  HESTON_BASE = {
    spot: 100.0, rate: 0.02, v0: 0.04, kappa: 1.5, theta: 0.04, xi: 0.5, rho: -0.7,
    strike: 100.0, maturity: 1.0
  }.freeze

  def test_gbm_milstein_has_higher_beta_than_euler
    r_eu = MlmcHeston::Mlmc.estimate_rates(
      model: :gbm, params: GBM_BASE.merge(scheme: :euler),
      levels: [1, 2, 3, 4], n_probe: 4_000, seed: 11, workers: 2
    )
    r_mi = MlmcHeston::Mlmc.estimate_rates(
      model: :gbm, params: GBM_BASE.merge(scheme: :milstein),
      levels: [1, 2, 3, 4], n_probe: 4_000, seed: 11, workers: 2
    )
    assert_operator r_mi[:beta_strong], :>, r_eu[:beta_strong]
    assert_operator r_mi[:beta_strong], :>, 1.5
  end

  def test_mlmc_matches_bs_and_heston_cf_anchors
    bs = MlmcHeston::Special.bs_call(100.0, 100.0, 0.05, 0.20, 1.0)
    res_g = MlmcHeston::Mlmc.run_adaptive(
      model: :gbm, params: GBM_BASE.merge(scheme: :milstein), eps: 0.04, seed: 22, workers: 2
    )
    assert_in_delta bs, res_g[:price], 0.15

    hc = MlmcHeston::HestonCf.call_price(HESTON_BASE, 100.0, 1.0)
    res_h = MlmcHeston::Mlmc.run_adaptive(
      model: :heston, params: HESTON_BASE.merge(scheme: :euler), eps: 0.04, seed: 33, workers: 2
    )
    assert_in_delta hc, res_h[:price], 0.15
  end

  def test_mlmc_seed_determinism_1_vs_2_workers
    a = MlmcHeston::Mlmc.run_adaptive(model: :heston, params: HESTON_BASE, eps: 0.05, seed: 77, workers: 1)
    b = MlmcHeston::Mlmc.run_adaptive(model: :heston, params: HESTON_BASE, eps: 0.05, seed: 77, workers: 2)
    assert_equal [a[:price]].pack('D'), [b[:price]].pack('D')
  end

  def test_complexity_slope_mlmc_lower_than_standard_mc
    comp = MlmcHeston::Complexity.compare_mlmc_vs_standard(
      model: :gbm, params: GBM_BASE.merge(scheme: :milstein),
      eps_grid: [0.08, 0.04, 0.02], seed: 55, workers: 2
    )
    assert_operator comp[:mlmc_cost_slope], :<, comp[:std_mc_cost_slope]
  end

  def test_antithetic_coupling_reduces_level_variance
    res = MlmcHeston::Antithetic.heston_antithetic_var(
      HESTON_BASE, 100.0, 1.0, level: 3, n_paths: 4_000, seed: 88, workers: 2
    )
    assert_operator res[:reduction_ratio], :>, 1.0
  end
end
