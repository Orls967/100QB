# frozen_string_literal: true

require 'minitest/autorun'
require 'copula_tail_dependence'

class TestCore < Minitest::Test
  def test_special_norm_cdf_and_inv_roundtrip
    [1.0e-6, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0 - 1.0e-6].each do |p|
      z = CopulaTailDependence::Special.norm_inv(p)
      assert_in_delta p, CopulaTailDependence::Special.norm_cdf(z), 1.0e-12
    end
  end

  def test_special_t_cdf_and_inv_roundtrip
    [3.0, 4.0, 10.0].each do |df|
      [0.01, 0.1, 0.5, 0.9, 0.99].each do |p|
        x = CopulaTailDependence::Special.t_inv(p, df)
        assert_in_delta p, CopulaTailDependence::Special.t_cdf(x, df), 1.0e-10
      end
    end
  end

  def test_bivariate_norm_cdf_anchors
    sp = CopulaTailDependence::Special
    # rho = 0 -> Phi(a)*Phi(b)
    assert_in_delta sp.norm_cdf(-1.2) * sp.norm_cdf(0.7), sp.bivariate_norm_cdf(-1.2, 0.7, 0.0), 1.0e-12
    # rho -> 1 -> min(Phi(a), Phi(b))
    assert_in_delta sp.norm_cdf(-1.5), sp.bivariate_norm_cdf(-1.5, 0.3, 1.0 - 1.0e-15), 1.0e-12
    # rho -> -1 -> max(0, Phi(a)+Phi(b)-1)
    assert_in_delta sp.norm_cdf(0.5) + sp.norm_cdf(0.5) - 1.0, sp.bivariate_norm_cdf(0.5, 0.5, -1.0 + 1.0e-15), 1.0e-12
    # Simetri P(a, b) = P(b, a) dan P(0, 0, rho) = 0.25 + asin(rho)/(2*pi)
    assert_in_delta sp.bivariate_norm_cdf(-1.0, 0.5, 0.6), sp.bivariate_norm_cdf(0.5, -1.0, 0.6), 1.0e-12
    assert_in_delta 0.25 + (Math.asin(0.5) / (2.0 * Math::PI)), sp.bivariate_norm_cdf(0.0, 0.0, 0.5), 1.0e-10
  end

  def test_positive_stable_laplace_transform
    rng = CopulaTailDependence::Rng::Sampler.new(12_345)
    alpha = 0.5
    n = 40_000
    s_val = 0.8
    sum = 0.0
    n.times do
      x = rng.next_positive_stable(alpha)
      sum += Math.exp(-s_val * x)
    end
    emp = sum / n
    exact = Math.exp(-(s_val**alpha))
    assert_in_delta exact, emp, 0.015
  end

  def test_cholesky_decomposition
    sigma = [[4.0, 2.0, 0.8], [2.0, 5.0, 1.5], [0.8, 1.5, 3.0]]
    l = CopulaTailDependence::Linalg.cholesky(sigma)
    3.times do |i|
      3.times do |j|
        val = (0...3).sum { |k| l[i][k] * l[j][k] }
        assert_in_delta sigma[i][j], val, 1.0e-12
      end
    end
  end

  def test_kendall_tau_and_spearman_on_copulas
    rng = CopulaTailDependence::Rng::Sampler.new(99_999)
    n = 12_000
    [
      [:gaussian, { rho: 0.6 }],
      [:t, { rho: 0.6, df: 4.0 }],
      [:clayton, { theta: 2.0 }],
      [:gumbel, { theta: 2.0 }]
    ].each do |type, params|
      u1 = Array.new(n, 0.0)
      u2 = Array.new(n, 0.0)
      n.times do |i|
        a, b = CopulaTailDependence::Copulas.sample_pair(type, params, rng)
        u1[i] = a
        u2[i] = b
      end
      tau_emp = CopulaTailDependence::Dependence.kendall_tau(u1, u2)
      tau_ex = CopulaTailDependence::Copulas.exact_kendall_tau(type, params)
      assert_in_delta tau_ex, tau_emp, 0.025, "#{type} tau"
    end
  end

  def test_parallel_worker_invariance_1_vs_4
    run_one = lambda do |w|
      chunks = CopulaTailDependence::Parallel.map_chunks(n_items: 500, workers: w, chunk_size: 100) do |cidx, len|
        rng = CopulaTailDependence::Rng::Sampler.new(CopulaTailDependence::Rng.derive_seed(42, 1, cidx))
        s = 0.0
        len.times { s += rng.next_normal }
        s
      end
      chunks.sum(0.0)
    end
    v1 = run_one.call(1)
    v4 = run_one.call(4)
    assert_equal [v1].pack('D'), [v4].pack('D')
  end
end
