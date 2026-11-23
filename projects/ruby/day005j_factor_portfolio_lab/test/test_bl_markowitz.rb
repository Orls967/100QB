# frozen_string_literal: true

require_relative 'test_helper'

class TestBlMarkowitz < Minitest::Test
  include FactorPortfolioLab

  def test_black_litterman_posterior_shrinks_cov
    world = FactorWorld.build_world(n_assets: 6, n_factors: 2, seed: 1201)
    sigma = world[:sigma]
    pi = BlackLitterman.implied_equilibrium(sigma, Array.new(6, 1.0 / 6))
    p_mat = [[1.0, -1.0, 0.0, 0.0, 0.0, 0.0]]
    q_vec = [0.02]
    omega = BlackLitterman.he_litterman_omega(sigma, p_mat, tau: 0.05)
    post = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: 0.05)

    6.times do |i|
      assert_operator post[:sigma_mu_post][i][i], :<=, (0.05 * sigma[i][i]) + 1e-15
      assert_operator post[:sigma_bl][i][i], :>=, sigma[i][i]
    end
  end

  def test_black_litterman_zero_surprise_keeps_prior
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 1202)
    sigma = world[:sigma]
    pi = BlackLitterman.implied_equilibrium(sigma, Array.new(5, 0.2))
    p_mat = [[1.0, -1.0, 0.0, 0.0, 0.0]]
    q_vec = [pi[0] - pi[1]]
    omega = BlackLitterman.he_litterman_omega(sigma, p_mat, tau: 0.05)
    post = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: 0.05)

    5.times { |i| assert_in_delta pi[i], post[:mu_post][i], 1e-12 }
  end

  def test_idzorek_confidence_monotone_interpolation
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 1203)
    sigma = world[:sigma]
    pi = BlackLitterman.implied_equilibrium(sigma, Array.new(5, 0.2))
    p_mat = [[1.0, 0.0, 0.0, 0.0, 0.0]]
    q_vec = [pi[0] + 0.02]

    low = BlackLitterman.posterior(
      sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec,
      omega: BlackLitterman.idzorek_omega(sigma, p_mat, [0.2]), tau: 0.05
    )
    high = BlackLitterman.posterior(
      sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec,
      omega: BlackLitterman.idzorek_omega(sigma, p_mat, [0.8]), tau: 0.05
    )
    assert_operator low[:mu_post][0], :>, pi[0]
    assert_operator high[:mu_post][0], :>, low[:mu_post][0]
    assert_operator high[:mu_post][0], :<, q_vec[0]
  end

  def test_bl_importance_mc_matches_exact_bayes
    study = BlStudy.run(master_seed: MAIN_SEED, n_mc: 80_000, n_assets: 5, tau: 0.05)
    assert_operator study[:max_mu_diff], :<, 2.5e-4
    assert_operator study[:max_cov_diff], :<, 1.0e-5
    assert_operator study[:ess], :>, 5000.0
  end

  def test_capped_simplex_projection_exact
    vec = [0.8, -0.3, 0.4, 0.1, 0.5]
    lo = Array.new(5, 0.0)
    hi = Array.new(5, 0.35)
    proj = Constraints.project_capped_simplex(vec, lower: lo, upper: hi)
    assert_in_delta 1.0, proj.sum, 1e-12
    5.times do |i|
      assert_operator proj[i], :>=, lo[i] - 1e-12
      assert_operator proj[i], :<=, hi[i] + 1e-12
    end
  end

  def test_capped_simplex_raises_infeasible
    assert_raises(ArgumentError) do
      Constraints.project_capped_simplex([0.5, 0.5], lower: [0.0, 0.0], upper: [0.3, 0.3])
    end
  end

  def test_unconstrained_gmv_and_two_fund_spanning
    world = FactorWorld.build_world(n_assets: 8, n_factors: 2, seed: 1204)
    sigma = world[:sigma]
    mu = world[:mu]
    w_gmv = Markowitz.unconstrained_gmv(sigma)
    w_tan = Markowitz.unconstrained_tangency(sigma, mu)
    assert_in_delta 1.0, w_gmv.sum, 1e-12
    assert_in_delta 1.0, w_tan.sum, 1e-12

    target = 0.5 * (LinAlg.dot(w_gmv, mu) + LinAlg.dot(w_tan, mu))
    w_mid = Markowitz.unconstrained_frontier_point(sigma, mu, target)
    assert_in_delta target, LinAlg.dot(w_mid, mu), 1e-12
    8.times { |i| assert_in_delta 0.5 * (w_gmv[i] + w_tan[i]), w_mid[i], 1e-10 }
  end

  def test_constrained_qp_kkt_long_only
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 1205)
    sol = Markowitz.solve_qp(world[:sigma], mu: world[:mu], kappa: 0.4)
    assert_in_delta 1.0, sol[:weights].sum, 1e-10
    assert(sol[:weights].all? { |w| w >= -1e-10 })
    assert_operator sol[:kkt][:max_kkt], :<, 1e-8
  end

  def test_constrained_qp_kkt_with_upper_cap
    world = FactorWorld.build_world(n_assets: 15, n_factors: 3, seed: 1206)
    lo = Array.new(15, 0.0)
    hi = Array.new(15, 0.12)
    sol = Markowitz.solve_qp(world[:sigma], mu: world[:mu], kappa: 0.8, lower: lo, upper: hi)
    assert_in_delta 1.0, sol[:weights].sum, 1e-10
    assert(sol[:weights].all? { |w| w.between?(-1e-10, 0.12 + 1e-10) })
    assert_operator sol[:kkt][:max_kkt], :<, 1e-8
  end

  def test_constrained_gmv_matches_unconstrained_when_all_positive
    sigma = [[0.04, 0.0, 0.0], [0.0, 0.09, 0.0], [0.0, 0.0, 0.16]]
    w_unc = Markowitz.unconstrained_gmv(sigma)
    sol = Markowitz.solve_qp(sigma, kappa: 0.0)
    3.times { |i| assert_in_delta w_unc[i], sol[:weights][i], 1e-9 }
  end

  def test_capped_frontier_has_higher_or_equal_vol_than_long_only
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 1207)
    s_lo = Markowitz.solve_qp(world[:sigma], kappa: 0.0, lower: Array.new(12, 0.0), upper: Array.new(12, 1.0))
    s_cap = Markowitz.solve_qp(world[:sigma], kappa: 0.0, lower: Array.new(12, 0.0), upper: Array.new(12, 0.12))
    assert_operator s_cap[:volatility], :>=, s_lo[:volatility] - 1e-10
  end

  def test_mv_study_45_points_all_kkt_clean
    study = MvStudy.run(master_seed: MAIN_SEED, n_assets: 15, n_points: 45)
    assert_equal 45, study[:long_only].length
    assert_equal 45, study[:capped].length
    assert_operator study[:max_kkt_long_only], :<, 1e-8
    assert_operator study[:max_kkt_capped], :<, 1e-8
  end

  def test_frontier_return_monotone_in_kappa
    study = MvStudy.run(master_seed: MAIN_SEED, n_assets: 12, n_points: 10)
    rets = study[:long_only].map { |p| p[:expected_return] }
    (1...rets.length).each { |i| assert_operator rets[i], :>=, rets[i - 1] - 1e-10 }
  end
end
