# frozen_string_literal: true

require 'minitest/autorun'
require 'copula_tail_dependence'

class TestStudies < Minitest::Test
  def test_euler_gaussian_sum_and_finite_diff
    ex = CopulaTailDependence::Euler.gaussian_exact
    assert_in_delta ex[:es], ex[:es_contribs].sum(0.0), 1.0e-12
    assert_in_delta ex[:var], ex[:var_contribs].sum(0.0), 1.0e-12

    fd = CopulaTailDependence::Euler.gaussian_es_finite_diff
    ex[:es_contribs].each_with_index do |c, i|
      assert_in_delta c, fd[i], 1.0e-7
    end
  end

  def test_euler_mc_matches_gaussian_exact_smoke
    ex = CopulaTailDependence::Euler.gaussian_exact
    mc = CopulaTailDependence::Euler.mc_euler_es(copula: :gaussian, n_paths: 40_000, seed: 101, workers: 2)
    assert_in_delta mc[:es], mc[:contribs].sum(0.0), 1.0e-10
    ex[:es_contribs].each_with_index do |c, i|
      z = (mc[:contribs][i] - c).abs / mc[:contrib_ses][i]
      assert_operator z, :<=, 3.5, "Euler contrib #{i} z=#{z}"
    end
  end

  def test_fit_gaussian_mle_and_t_profile
    rng = CopulaTailDependence::Rng::Sampler.new(7_777)
    n = 600
    u1 = Array.new(n, 0.0)
    u2 = Array.new(n, 0.0)
    n.times do |i|
      a, b = CopulaTailDependence::Copulas.sample_pair(:t, { rho: 0.5, df: 4.0 }, rng)
      u1[i] = a
      u2[i] = b
    end
    g = CopulaTailDependence::Fit.fit_gaussian_mle(u1, u2)
    t_prof = CopulaTailDependence::Fit.fit_t_profile_mle(u1, u2, df_grid: [3.0, 4.0, 10.0, Float::INFINITY])
    assert_in_delta 0.5, g[:rho], 0.08
    assert_operator t_prof[:best][:loglik], :>, g[:loglik]
  end

  def test_importance_sampling_matches_exact_anchors
    res2 = CopulaTailDependence::Importance.compare_naive_vs_is(
      d: 2, a: -3.0, rho: 0.5, n_paths: 20_000, seed: 55, workers: 2
    )
    assert_operator res2[:is_z].abs, :<=, 3.0
    assert_operator res2[:variance_reduction], :>, 50.0

    res5 = CopulaTailDependence::Importance.compare_naive_vs_is(
      d: 5, a: -3.0, rho: 0.8, n_paths: 20_000, seed: 56, workers: 2
    )
    assert_operator res5[:is_z].abs, :<=, 3.0
    assert_operator res5[:variance_reduction], :>, 100.0
  end

  def test_seed_determinism_importance_sampling
    a = CopulaTailDependence::Importance.compare_naive_vs_is(
      d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: 99, workers: 2
    )
    b = CopulaTailDependence::Importance.compare_naive_vs_is(
      d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: 99, workers: 1
    )
    assert_equal [a[:is_mean]].pack('D'), [b[:is_mean]].pack('D')
  end
end
