# frozen_string_literal: true

require 'minitest/autorun'
require 'garch_family'

class TestCoreModels < Minitest::Test
  def test_splitmix64_and_parallel_worker_invariance
    s1 = GarchFamily::Rng.derive_seed(GarchFamily::MAIN_SEED, 1, 0)
    s2 = GarchFamily::Rng.derive_seed(GarchFamily::MAIN_SEED, 1, 1)
    refute_equal s1, s2
    run_w = lambda do |w|
      GarchFamily::Parallel.map_chunks(n_items: 200, workers: w, chunk_size: 50) do |cidx, len|
        rng = GarchFamily::Rng::Sampler.new(GarchFamily::Rng.derive_seed(42, 7, cidx))
        Array.new(len) { rng.next_std_t(5.0) }.sum(0.0)
      end.sum(0.0)
    end
    assert_equal [run_w.call(1)].pack('D'), [run_w.call(4)].pack('D')
  end

  def test_special_normal_and_student_t_inv_roundtrip
    [1.0e-6, 0.01, 0.05, 0.5, 0.99].each do |p|
      assert_in_delta p, GarchFamily::Special.norm_cdf(GarchFamily::Special.norm_inv(p)), 1.0e-12
      assert_in_delta p, GarchFamily::Special.student_t_cdf(GarchFamily::Special.student_t_inv(p, 5.0), 5.0), 1.0e-9
    end
  end

  def test_student_t_expected_shortfall_matches_quadrature
    dist = GarchFamily::Innovations::StudentT.new(nu: 5.0)
    scale = Math.sqrt(5.0 / 3.0)
    [0.01, 0.05].each do |alpha|
      es_closed = dist.expected_shortfall(alpha)
      zq = dist.var_quantile(alpha)
      es_quad = [[zq, 25.0], [25.0, 200.0], [200.0, 5000.0]].sum(0.0) do |lo, hi|
        nodes, weights = GarchFamily::Special.gauss_legendre(64, lo, hi)
        nodes.zip(weights).sum(0.0) { |z, w| w * z * scale * GarchFamily::Special.student_t_pdf(z * scale, 5.0) }
      end / alpha
      assert_in_delta es_quad, es_closed, 1.0e-8
    end
  end

  def test_skew_t_density_integrates_to_one_with_unit_variance
    a, b, = GarchFamily::Innovations::SkewT.constants(6.0, -0.3)
    kink = -a / b
    mass = m2 = 0.0
    [[-400.0, -40.0], [-40.0, kink], [kink, 40.0], [40.0, 400.0]].each do |lo, hi|
      nodes, weights = GarchFamily::Special.gauss_legendre(64, lo, hi)
      nodes.zip(weights).each do |z, w|
        pdf = Math.exp(GarchFamily::Innovations::SkewT.log_pdf(z, 6.0, -0.3))
        mass += w * pdf
        m2 += w * z * z * pdf
      end
    end
    assert_in_delta 1.0, mass, 1.0e-7
    assert_in_delta 1.0, m2, 1.0e-6
  end

  def test_linalg_solve_invert_and_newey_west
    a = [[4.0, 1.0], [1.0, 3.0]]
    inv = GarchFamily::Linalg.invert(a)
    prod = GarchFamily::Linalg.mat_mul(a, inv)
    assert_in_delta 1.0, prod[0][0], 1.0e-12
    assert_in_delta 0.0, prod[0][1], 1.0e-12
    ols = GarchFamily::Linalg.ols_with_nw([[1.0, 1.0], [1.0, 2.0], [1.0, 3.0], [1.0, 4.0]], [3.0, 5.0, 7.0, 9.0])
    assert_in_delta 1.0, ols[:beta][0], 1.0e-12
    assert_in_delta 2.0, ols[:beta][1], 1.0e-12
  end

  def test_chi2_mixture_critical_value_is_2_7055
    p_at_crit = GarchFamily::Special.chi2_mixture_01_sf(2.705543454)
    assert_in_delta 0.05, p_at_crit, 1.0e-6
  end
end
