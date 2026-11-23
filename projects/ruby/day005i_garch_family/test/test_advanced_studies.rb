# frozen_string_literal: true

require 'minitest/autorun'
require 'garch_family'

class TestAdvancedStudies < Minitest::Test
  def test_patton_qlike_minimized_at_zero_when_forecast_equals_proxy
    assert_in_delta 0.0, GarchFamily::ModelComparison.qlike_loss(0.01, 0.01), 1.0e-15
    assert_operator GarchFamily::ModelComparison.qlike_loss(0.015, 0.01), :>, 0.0
    assert_operator GarchFamily::ModelComparison.qlike_loss(0.005, 0.01), :>, 0.0
  end

  def test_diebold_mariano_detects_superior_forecast_with_nw_se
    rng = GarchFamily::Rng::Sampler.new(2001)
    loss_bad = Array.new(500) { 1.5 + (0.2 * rng.next_normal) }
    loss_good = Array.new(500) { 1.0 + (0.2 * rng.next_normal) }
    dm = GarchFamily::ModelComparison.diebold_mariano(loss_bad, loss_good)
    assert_operator dm[:dm_stat], :>, 10.0
    assert_operator dm[:p_value], :<, 1.0e-10
  end

  def test_dcc_zero_params_reduces_identically_to_ccc
    rng = GarchFamily::Rng::Sampler.new(2002)
    z_rows = Array.new(100) { [rng.next_normal, rng.next_normal] }
    ccc = GarchFamily::Dcc.filter_corr(0.0, 0.0, z_rows)
    ccc[:corr_seq].each do |r|
      assert_in_delta ccc[:q_bar][0][1], r[0][1], 1.0e-15
    end
  end

  def test_dcc_cholesky_and_filtered_corr_remain_positive_definite
    anch = GarchFamily::Dcc.ccc_reduction_and_pd_anchor(GarchFamily::MAIN_SEED)
    assert_operator anch[:max_ccc_err], :<, 1.0e-13
    assert_operator anch[:min_det], :>, 0.0
    assert_operator anch[:max_abs_rho], :<, 1.0
    assert_operator anch[:sample_corr_true_vs_filt], :>, 0.5
  end

  def test_dcc_fit_stage2_recovers_positive_persistence
    sim = GarchFamily::Dcc.simulate_bivariate_dcc(1500, seed: 2003)
    g = GarchFamily::Garch11.new
    v1 = g.filter(GarchFamily::Dcc::M1_PARAMS, sim[:r1])
    v2 = g.filter(GarchFamily::Dcc::M2_PARAMS, sim[:r2])
    zs = Array.new(1500) { |t| [sim[:r1][t] / Math.sqrt(v1[t]), sim[:r2][t] / Math.sqrt(v2[t])] }
    fit = GarchFamily::Dcc.fit_stage2(zs)
    assert fit[:converged]
    assert_operator fit[:persistence], :<, 1.0
    assert_operator fit[:a], :>, 0.01
  end

  def test_har_rv_infill_anchor_unbiased_and_ols_design
    sim = GarchFamily::HarRv.simulate_intraday_sv(n_days: 1000, m_bars: 78, seed: 2004)
    diffs = sim[:rv].each_with_index.map { |rv, i| rv - sim[:iv][i] }
    mean_d = diffs.sum(0.0) / diffs.length
    se_d = Math.sqrt(diffs.sum(0.0) { |x| (x - mean_d)**2 } / (diffs.length - 1)) / Math.sqrt(diffs.length)
    assert_operator mean_d.abs / se_d, :<=, 3.0
    x, y = GarchFamily::HarRv.har_design(sim[:rv], 22, 500)
    assert_equal 478, x.length
    assert_equal 4, x.first.length
    assert_equal 478, y.length
  end

  def test_har_vs_garch_outperforms_daily_garch_on_qlike
    res = GarchFamily::HarRv.har_vs_garch_study(GarchFamily::MAIN_SEED)
    assert_operator res[:qlike_har_vs_iv], :<, res[:qlike_garch_vs_iv]
    assert_operator res[:dm_garch_minus_har][:p_value], :<, 0.05
  end

  def test_structural_break_cusum_detects_omega_jump_near_midpoint
    sim = GarchFamily::MisspecPortfolio.simulate_break_series(n_obs: 3000, seed: 2005)
    cusum = GarchFamily::MisspecPortfolio.inclan_tiao_cusum(sim[:returns])
    assert cusum[:break_detected]
    assert_in_delta 1500.0, cusum[:detected_idx].to_f, 300.0
  end

  def test_five_asset_euler_risk_decomposition_sums_exactly_to_port_sd
    port = GarchFamily::MisspecPortfolio.run_five_asset_portfolio(seed: GarchFamily::MAIN_SEED)
    calm = port[:stress_test][:calm]
    stress = port[:stress_test][:stress]
    assert_operator calm[:euler_sum_err], :<, 1.0e-13
    assert_operator stress[:euler_sum_err], :<, 1.0e-13
    assert_operator calm[:div_ratio], :>, stress[:div_ratio]
  end
end
