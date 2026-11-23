# frozen_string_literal: true

require_relative 'test_helper'

class TestLinalgDist < Minitest::Test
  include FactorPortfolioLab

  def test_rng_deterministic_streams
    s1 = Rng.derive_seed(MAIN_SEED, 1, 0)
    s2 = Rng.derive_seed(MAIN_SEED, 1, 0)
    s3 = Rng.derive_seed(MAIN_SEED, 2, 0)
    assert_equal s1, s2
    refute_equal s1, s3
  end

  def test_rng_std_t_and_wishart_moments
    rng = Rng::Sampler.new(12_345)
    draws = Array.new(6000) { rng.next_std_t(6.0) }
    mean = draws.sum / draws.length
    var = draws.sum { |x| (x - mean)**2 } / (draws.length - 1)
    assert_in_delta 0.0, mean, 0.06
    assert_in_delta 1.0, var, 0.12

    scale = [[1.0, 0.3], [0.3, 1.0]]
    w = rng.next_wishart(scale, 10)
    assert_operator w[0][0], :>, 0.0
    assert_in_delta w[0][1], w[1][0], 1e-12
  end

  def test_norm_cdf_and_inv_roundtrip
    [-3.2, -1.96, -0.5, 0.0, 0.75, 1.96, 3.2].each do |x|
      p = Distributions.norm_cdf(x)
      assert_in_delta x, Distributions.norm_inv(p), 1e-8
    end
  end

  def test_chi2_and_f_distributions
    p_chi = Distributions.chi2_sf(3.84145882, 1)
    assert_in_delta 0.05, p_chi, 1e-4

    f_val = Distributions.f_inv(0.95, 10, 40)
    assert_in_delta 0.05, Distributions.f_sf(f_val, 10, 40), 1e-7
    assert_in_delta f_val, 1.0 / Distributions.f_inv(0.05, 40, 10), 1e-6
  end

  def test_student_t_roundtrip
    [-2.5, -1.0, 0.0, 1.2, 2.8].each do |t|
      p = Distributions.student_t_cdf(t, 7.0)
      assert_in_delta t, Distributions.student_t_inv(p, 7.0), 1e-7
    end
  end

  def test_root_finders
    r1 = RootFinder.brent(1.0, 2.0) { |x| (x * x) - 2.0 }
    assert_in_delta Math.sqrt(2.0), r1, 1e-9

    r2 = RootFinder.newton(0.7) { |x| [Math.cos(x) - x, -Math.sin(x) - 1.0] }
    assert_in_delta 0.7390851332, r2, 1e-9
  end

  def test_solve_invert_and_cholesky
    a = [[4.0, 1.2, -0.4], [1.2, 3.0, 0.5], [-0.4, 0.5, 2.5]]
    inv_a = LinAlg.invert(a)
    prod = LinAlg.mat_mul(a, inv_a)
    3.times do |i|
      3.times { |j| assert_in_delta(i == j ? 1.0 : 0.0, prod[i][j], 1e-10) }
    end

    chol = LinAlg.cholesky(a)
    rec = LinAlg.mat_mul(chol, LinAlg.transpose(chol))
    3.times do |i|
      3.times { |j| assert_in_delta a[i][j], rec[i][j], 1e-10 }
    end

    rhs = [1.0, -2.0, 3.0]
    x_spd = LinAlg.solve_spd(a, rhs)
    x_gen = LinAlg.solve(a, rhs)
    3.times { |i| assert_in_delta x_gen[i], x_spd[i], 1e-10 }
  end

  def test_qr_and_lu_det
    a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]
    q, r = LinAlg.qr_householder(a)
    qtq = LinAlg.mat_mul(LinAlg.transpose(q), q)
    3.times do |i|
      3.times { |j| assert_in_delta(i == j ? 1.0 : 0.0, qtq[i][j], 1e-10) }
    end
    rec = LinAlg.mat_mul(q, r)
    3.times do |i|
      3.times { |j| assert_in_delta a[i][j], rec[i][j], 1e-10 }
    end
    assert_in_delta 8.0, LinAlg.lu_det(a), 1e-10
  end

  def test_jacobi_sym_eigen_decomposition
    a = [[5.0, 1.5, -0.5], [1.5, 3.0, 0.8], [-0.5, 0.8, 2.0]]
    eig = LinAlg.jacobi_sym(a)
    vals = eig[:eigenvalues]
    vecs = eig[:eigenvectors]
    assert_operator vals[0], :>=, vals[1]
    assert_operator vals[1], :>=, vals[2]
    assert_in_delta LinAlg.trace(a), vals.sum, 1e-10

    vtv = LinAlg.mat_mul(LinAlg.transpose(vecs), vecs)
    3.times do |i|
      3.times { |j| assert_in_delta(i == j ? 1.0 : 0.0, vtv[i][j], 1e-10) }
    end
  end

  def test_parallel_map_items_ordered
    res = Parallel.map_items(n_items: 45, workers: 2, chunk_size: 10) { |i| (i * i) + 1 }
    expected = Array.new(45) { |i| (i * i) + 1 }
    assert_equal expected, res
  end
end
