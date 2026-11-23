# frozen_string_literal: true

require 'minitest/autorun'
require 'backtest_overfitting'

class TestCore < Minitest::Test
  def test_splitmix64_seed_determinism_and_avalanche
    s1 = BacktestOverfitting::Rng.derive_seed(BacktestOverfitting::MAIN_SEED, 1, 0)
    s2 = BacktestOverfitting::Rng.derive_seed(BacktestOverfitting::MAIN_SEED, 1, 1)
    assert_equal s1, BacktestOverfitting::Rng.derive_seed(BacktestOverfitting::MAIN_SEED, 1, 0)
    refute_equal s1, s2
    assert_operator (s1 ^ s2).to_s(2).count('1'), :>=, 15
  end

  def test_parallel_worker_invariance_1_vs_4
    run_one = lambda do |w|
      BacktestOverfitting::Parallel.map_chunks(n_items: 400, workers: w, chunk_size: 100) do |cidx, len|
        rng = BacktestOverfitting::Rng::Sampler.new(BacktestOverfitting::Rng.derive_seed(99, 1, cidx))
        s = 0.0
        len.times { s += rng.next_normal }
        s
      end.sum(0.0)
    end
    assert_equal [run_one.call(1)].pack('D'), [run_one.call(4)].pack('D')
  end

  def test_norm_cdf_and_norm_inv_roundtrip
    [1.0e-12, 1.0e-6, 0.01, 0.25, 0.5, 0.85, 0.999, 1.0 - 1.0e-12].each do |p|
      x = BacktestOverfitting::Special.norm_inv(p)
      assert_in_delta p, BacktestOverfitting::Special.norm_cdf(x), 1.0e-13
    end
  end

  def test_exact_binomial_combinations_with_big_integers
    assert_equal 12_870, BacktestOverfitting::Special.binom(16, 8)
    assert_equal 155_117_520, BacktestOverfitting::Special.binom(30, 15)
    assert_equal 1, BacktestOverfitting::Special.binom(20, 0)
    assert_equal 0, BacktestOverfitting::Special.binom(10, 12)
  end

  def test_gauss_legendre_integrates_polynomial_degree_9_exactly
    nodes, weights = BacktestOverfitting::Special.gauss_legendre(10, 0.0, 2.0)
    approx = nodes.zip(weights).sum(0.0) { |x, w| w * ((x**7) - (3.0 * (x**3)) + 2.0) }
    exact = ((2.0**8) / 8.0) - (3.0 * (2.0**4) / 4.0) + 4.0
    assert_in_delta exact, approx, 1.0e-12
  end

  def test_student_t_pvalue_matches_quadrature_at_known_quantiles
    assert_in_delta 1.0, BacktestOverfitting::Special.student_t_pvalue(0.0, 10), 1.0e-12
    assert_in_delta 0.05, BacktestOverfitting::Special.student_t_pvalue(2.2281388519, 10), 1.0e-8
    assert_in_delta 0.05, BacktestOverfitting::Special.student_t_pvalue(1.96, 5_000), 1.0e-4
  end

  def test_ranks_handles_ties_properly
    assert_equal [4.0, 1.0, 2.5, 2.5], BacktestOverfitting::Special.ranks([30.0, 10.0, 20.0, 20.0])
  end
end
