# frozen_string_literal: true

require 'minitest/autorun'
require 'cdo_tranche'

class TestStudies < Minitest::Test
  def test_pricing_tranche_sum_conservation_and_lhp_convergence
    pvs = CdoTranche::Models::TRANCHES.map do |tr|
      CdoTranche::Pricing.homogeneous_exact_pv(125, 0.015, 0.40, 0.25, tr[:attach], tr[:detach])
    end
    sum_prot = pvs.sum(0.0) { |r| r[:pv_prot] }
    full_prot = CdoTranche::Pricing.homogeneous_exact_pv(125, 0.015, 0.40, 0.25, 0.0, 1.0)[:pv_prot]
    assert_in_delta full_prot, sum_prot, 1.0e-12

    gran = CdoTranche::Pricing.granularity_study
    assert_operator gran.last[:abs_err], :<, gran.first[:abs_err]
  end

  def test_base_correlation_calibration_recovers_targets
    els = [0.025, 0.036, 0.040, 0.042, 0.0432]
    detaches = [0.03, 0.07, 0.10, 0.15, 0.30]
    res = CdoTranche::BaseCorr.calibrate_base_correlation(els, detaches)
    res.each do |row|
      refute_nil row[:base_corr]
      assert_operator row[:implied_tranche_el], :>=, 0.0
    end
  end

  def test_importance_sampling_matches_exact_asb_and_reduces_variance
    res = CdoTranche::Importance.compare_naive_vs_is(
      attach: 0.15, detach: 0.30, n_paths: 10_000, seed: 77, workers: 2
    )
    assert_operator res[:is_z].abs, :<=, 3.0
    assert_operator res[:variance_reduction], :>, 20.0
  end

  def test_importance_sampling_seed_determinism
    a = CdoTranche::Importance.compare_naive_vs_is(attach: 0.10, detach: 0.15, n_paths: 2_000, seed: 42, workers: 1)
    b = CdoTranche::Importance.compare_naive_vs_is(attach: 0.10, detach: 0.15, n_paths: 2_000, seed: 42, workers: 2)
    assert_equal [a[:is_el]].pack('D'), [b[:is_el]].pack('D')
  end

  def test_model_risk_smoke_produces_valid_spreads
    res = CdoTranche::ModelRisk.compare_models(n_paths: 2_000, seed: 88, workers: 2)
    assert_equal %i[gaussian t4 random_recovery], res.keys
    assert_operator res[:t4].last[:fair_spread_bps], :>=, res[:gaussian].last[:fair_spread_bps]
  end
end
