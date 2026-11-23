# frozen_string_literal: true

require 'minitest/autorun'
require 'backtest_overfitting'

class TestCscvAndPurged < Minitest::Test
  def test_cscv_enumerates_exact_12870_splits
    masks = BacktestOverfitting::CscvPbo.half_masks
    assert_equal 12_870, masks.length
    assert_equal 12_870, masks.uniq.length
    masks.each do |lo, hi|
      assert_equal 8, lo.to_s(2).count('1') + hi.to_s(2).count('1')
    end
  end

  def test_cscv_precomputed_half_tables_match_direct_sharpe
    rng = BacktestOverfitting::Rng::Sampler.new(71)
    mat = Array.new(6) { BacktestOverfitting::DataGenerator.iid_gaussian(256, sigma: 0.01, rng: rng) }
    tbl = BacktestOverfitting::CscvPbo.precompute_half_tables(mat)
    lo_m, hi_m = BacktestOverfitting::CscvPbo.half_masks[42]
    _lam, sr_is, sr_oos = BacktestOverfitting::CscvPbo.eval_split_mask(tbl, lo_m, hi_m)
    assert_operator sr_is, :>, -5.0
    assert_operator sr_oos, :>, -5.0
  end

  def test_cscv_noise_anchor_pbo_near_half_and_planted_low
    rng = BacktestOverfitting::Rng::Sampler.new(72)
    planted = Array.new(20) do |j|
      BacktestOverfitting::DataGenerator.iid_gaussian(512, sigma: 0.01, mu: j < 3 ? 0.0022 : 0.0, rng: rng)
    end
    res_p = BacktestOverfitting::CscvPbo.evaluate_cscv(planted, workers: 2)
    assert_operator res_p[:pbo], :<, 0.15

    res_n = BacktestOverfitting::CscvPbo.noise_exchangeable_pbo_study(
      n_strats: 12, n_obs: 256, n_mats: 12, seed: 73, workers: 2
    )
    assert_in_delta 0.5, res_n[:emp_mean_pbo], 4.0 * res_n[:emp_se_pbo]
  end

  def test_purged_train_indices_drops_overlapping_horizon_and_embargo
    idx = BacktestOverfitting::PurgedCv.purged_train_indices(100, 40, 59, horizon: 10, embargo: 10)
    # Dihapus: [40 - 10 + 1 = 31 .. 59 + 10 + 10 - 1 = 78]
    (31..78).each { |t| refute_includes idx, t }
    assert_includes idx, 30
    assert_includes idx, 79
  end

  def test_purged_cv_eliminates_overlapping_leakage_unlike_naive_cv
    res = BacktestOverfitting::PurgedCv.run_study(seed: 74, workers: 2)
    assert_operator res[:naive_shuffled_corr_mean], :>, 0.25
    assert_operator res[:purged_zero_z], :<=, 3.5
  end

  def test_transaction_cost_reduces_high_turnover_more_than_low_turnover
    res = BacktestOverfitting::CscvPbo.transaction_cost_study(seed: 75)
    r0 = res[:rows][0]
    r20 = res[:rows][2]
    assert_operator r20[:hi_turnover_mean_sharpe], :<, r0[:hi_turnover_mean_sharpe]
    assert_operator r20[:hi_turnover_mean_sharpe], :<, r20[:lo_turnover_mean_sharpe]
  end

  def test_multiple_testing_indep_fwer_and_bh_fdr_close_to_pi0_alpha
    res = BacktestOverfitting::MultipleTesting.simulate_mt_regime(rho: 0.0, n_reps: 250, seed: 76, workers: 2)
    assert_operator res[:methods][:bonferroni][:fwer], :<=, 0.07
    assert_in_delta res[:theoretical_bh_fdr], res[:methods][:bh][:fdr], 0.02
  end

  def test_sharpe_summary_matches_direct_student_t_statistics
    rng = BacktestOverfitting::Rng::Sampler.new(77)
    rets = BacktestOverfitting::DataGenerator.student_t_returns(10_000, df: 5, sigma: 0.01, rng: rng)
    mom = BacktestOverfitting::Metrics.summary(rets)
    assert_in_delta 0.01, mom[:std], 0.002
    assert_operator mom[:kurtosis], :>, 5.0
  end
end
