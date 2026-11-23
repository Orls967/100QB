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

  def test_telescoping_expectation_identity_gbm
    params = GBM_BASE.merge(scheme: :milstein)
    sum_levels = (0..3).sum(0.0) do |l|
      MlmcHeston::Mlmc.level_stats(model: :gbm, params: params, level: l, n_paths: 4_000,
                                   seed: 101, workers: 2)[:mean_y]
    end
    fine = MlmcHeston::Mlmc.level_stats(model: :gbm, params: params, level: 3, n_paths: 4_000, seed: 101, workers: 2)
    assert_in_delta fine[:sum_pf] / 4_000.0, sum_levels, 0.25
  end

  def test_telescoping_expectation_identity_heston
    params = HESTON_BASE.merge(scheme: :euler)
    sum_levels = (0..2).sum(0.0) do |l|
      MlmcHeston::Mlmc.level_stats(model: :heston, params: params, level: l, n_paths: 4_000,
                                   seed: 102, workers: 2)[:mean_y]
    end
    hc = MlmcHeston::HestonCf.call_price(HESTON_BASE, 100.0, 1.0)
    assert_in_delta hc, sum_levels, 0.35
  end

  def test_level_kurtosis_finite_and_bounded
    rng = MlmcHeston::Rng::Sampler.new(303)
    ys = Array.new(2_000) { MlmcHeston::Simulators.heston_coupled_pair(HESTON_BASE, 100.0, 1.0, 2, scheme: :euler, payoff_type: :european, rng: rng)[1] }
    m = ys.sum(0.0) / ys.length
    m2 = ys.sum(0.0) { |y| (y - m)**2 } / ys.length
    m4 = ys.sum(0.0) { |y| (y - m)**4 } / ys.length
    kurt = m4 / (m2 * m2)
    assert_operator kurt, :>, 1.0
    assert_operator kurt, :<, 50.0
  end

  def test_giles_theoretical_slope_prediction_formula
    assert_in_delta 2.0, MlmcHeston::Complexity.giles_asymptotic_slope(1.9, 1.0), 1.0e-12
    assert_in_delta 2.2, MlmcHeston::Complexity.giles_asymptotic_slope(0.8, 1.0), 1.0e-12
  end

  def test_antithetic_pair_symmetry_at_zero_vol_of_vol
    p_det = HESTON_BASE.merge(xi: 0.0)
    res = MlmcHeston::Antithetic.heston_antithetic_var(
      p_det, 100.0, 1.0, level: 2, n_paths: 1_000, seed: 44, workers: 1
    )
    assert_operator res[:var_antithetic], :<, res[:var_standard]
  end

  def test_antithetic_gbm_euler_ratio_is_unity_and_cir_decays_across_levels
    tbl = MlmcHeston::Antithetic.multi_model_level_table(
      GBM_BASE, HESTON_BASE, levels: [1, 6], n_paths: 6_000, seed: 52, workers: 2
    )
    assert_in_delta 1.0, tbl.first[:gbm_euler_ratio], 1.0e-9
    assert_in_delta 1.0, tbl.last[:gbm_euler_ratio], 1.0e-9
    assert_operator tbl.last[:cir_euler_ratio], :<=, 1.15
    assert_operator tbl.last[:heston_euler_ratio], :<=, tbl.last[:heston_theoretical_ceiling]
  end
end
