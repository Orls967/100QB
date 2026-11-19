# frozen_string_literal: true

require_relative 'test_helper'

class TestMathAndSolver < Minitest::Test
  def test_norm_cdf_and_bs_put
    assert_in_delta 0.5, AmericanPut::Stats.norm_cdf(0.0), 1e-15
    assert_in_delta 0.8413447460685429, AmericanPut::Stats.norm_cdf(1.0), 1e-14

    bs = AmericanPut::Stats.bs_put(AmericanPut::DEFAULT_PARAMS)
    assert_in_delta 5.573526022256971, bs, 1e-12
  end

  def test_box_muller_determinism_and_moments
    rng1 = AmericanPut::Rng::NormalGenerator.new(123)
    rng2 = AmericanPut::Rng::NormalGenerator.new(123)
    seq = Array.new(20) { rng1.next_normal }
    buf = Array.new(20, 0.0)
    rng2.fill_normals!(buf)

    seq.each_with_index { |val, idx| assert_in_delta val, buf[idx], 1e-15 }

    big = Array.new(20_000, 0.0)
    rng1.fill_normals!(big)
    stats = AmericanPut::Stats.summarize_samples(big)
    assert_in_delta 0.0, stats[:price], 0.03
    assert_in_delta 1.0 / Math.sqrt(20_000), stats[:se], 0.002
  end

  def test_linear_solver_dense_and_tridiagonal
    a = [[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]]
    b = [8.0, -11.0, -3.0]
    x = AmericanPut::LinearSolver.solve(a, b)
    assert_in_delta 2.0, x[0], 1e-12
    assert_in_delta 3.0, x[1], 1e-12
    assert_in_delta(-1.0, x[2], 1e-12)

    lower = [0.0, -1.0, -1.0]
    diag = [2.0, 2.0, 2.0]
    upper = [-1.0, -1.0, 0.0]
    rhs = [1.0, 0.0, 1.0]
    xt = AmericanPut::LinearSolver.solve_tridiagonal(lower, diag, upper, rhs)
    xt.each { |v| assert_in_delta 1.0, v, 1e-12 }
  end
end
