# frozen_string_literal: true

require 'minitest/autorun'
require 'backtest_overfitting'

class TestSharpeAndSelection < Minitest::Test
  def test_delta_method_se_matches_gaussian_mc
    res = BacktestOverfitting::SharpeStats.gaussian_coverage_mc(
      n_obs: 500, sr_true: 0.10, n_reps: 4_000, seed: 11, workers: 2
    )
    assert_in_delta res[:theo_se], res[:emp_sd], 0.005
    assert_in_delta 0.95, res[:coverage_95], 0.025
  end

  def test_nongaussian_student_t_kurtosis_widens_se_above_gaussian
    se_g = BacktestOverfitting::SharpeStats.se_sharpe(0.30, 1_000, skew: 0.0, kurt: 3.0)
    se_t5 = BacktestOverfitting::SharpeStats.se_sharpe(0.30, 1_000, skew: 0.0, kurt: 9.0)
    assert_operator se_t5, :>, se_g * 1.05
  end

  def test_lo_autocorr_factor_identity_and_monotone
    assert_in_delta 1.0, BacktestOverfitting::SharpeStats.lo_autocorr_factor(21, 0.0), 1.0e-12
    eta_pos = BacktestOverfitting::SharpeStats.lo_autocorr_factor(21, 0.3)
    eta_neg = BacktestOverfitting::SharpeStats.lo_autocorr_factor(21, -0.3)
    assert_operator eta_pos, :<, 0.7
    assert_operator eta_neg, :>, 1.3
  end

  def test_stationary_bootstrap_ci_covers_ar1_sharpe
    rng = BacktestOverfitting::Rng::Sampler.new(22)
    rets = BacktestOverfitting::DataGenerator.ar1_gaussian(400, phi: 0.3, sigma: 0.01, mu: 0.0005, rng: rng)
    lo, hi = BacktestOverfitting::SharpeStats.stationary_bootstrap_sharpe_ci(
      rets, mean_block_len: 10.0, n_boot: 300, rng: rng
    )
    sr = BacktestOverfitting::Metrics.sharpe_per_period(rets)
    assert_operator lo, :<, sr
    assert_operator sr, :<, hi
  end

  def test_expected_max_exact_matches_euler_mascheroni_asymptotically
    e300 = BacktestOverfitting::SelectionBias.expected_max_exact(300)
    a300 = BacktestOverfitting::SelectionBias.expected_max_approx(300)
    assert_in_delta 3.0, e300, 0.25
    assert_in_delta e300, a300, 0.08
  end

  def test_deflated_sharpe_penalizes_threshold_and_negative_skew
    dsr_base = BacktestOverfitting::SelectionBias.deflated_sharpe_ratio(0.20, 0.15, 252, skew: 0.0, kurt: 3.0)
    dsr_skew = BacktestOverfitting::SelectionBias.deflated_sharpe_ratio(0.20, 0.15, 252, skew: -1.0, kurt: 8.0)
    assert_operator dsr_base, :>, 0.70
    assert_operator dsr_skew, :<, dsr_base
  end

  def test_bonferroni_holm_and_bh_ordering_on_pvalues
    pvals = [0.0001, 0.004, 0.012, 0.04, 0.20, 0.50]
    bonf = BacktestOverfitting::MultipleTesting.bonferroni_rejections(pvals, alpha: 0.05)
    holm = BacktestOverfitting::MultipleTesting.holm_rejections(pvals, alpha: 0.05)
    bh = BacktestOverfitting::MultipleTesting.benjamini_hochberg_rejections(pvals, alpha: 0.05)
    assert_operator bonf.length, :<=, holm.length
    assert_operator holm.length, :<=, bh.length
    assert_equal [0, 1], bonf
    assert_equal [0, 1, 2], bh
  end

  def test_romano_wolf_rejects_strong_signal_and_controls_null
    rng = BacktestOverfitting::Rng::Sampler.new(33)
    centered = Array.new(120) { Array.new(10) { rng.next_normal.abs } }
    obs = [5.5, 0.2, 0.3, 0.1, 0.4, 0.5, 0.2, 0.1, 0.3, 0.2]
    rej = BacktestOverfitting::MultipleTesting.romano_wolf_rejections(centered, obs, alpha: 0.05)
    assert_equal [0], rej
  end
end
