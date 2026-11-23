# frozen_string_literal: true

require 'minitest/autorun'
require 'garch_family'

class TestMomentsAndScores < Minitest::Test
  def test_garch11_unconditional_variance_kurtosis_and_acf
    g = GarchFamily::Garch11.new
    p = { mu: 0.0, omega: 1.0e-5, alpha: 0.10, beta: 0.80 }
    assert_in_delta 1.0e-4, g.unconditional_variance(p), 1.0e-12
    assert_in_delta 3.0 * 0.19 / 0.17, g.theoretical_kurtosis(p), 1.0e-12
    rho1 = 0.10 * (1.0 - (0.80 * 0.90)) / (1.0 - (2.0 * 0.80 * 0.90) + 0.64)
    assert_in_delta rho1, g.theoretical_acf(p, 1), 1.0e-12
    assert_in_delta rho1 * 0.90, g.theoretical_acf(p, 2), 1.0e-12
  end

  def test_garch11_multistep_forecast_converges_to_unconditional_variance
    g = GarchFamily::Garch11.new
    p = { mu: 0.0, omega: 2.0e-6, alpha: 0.08, beta: 0.88 }
    assert_in_delta g.unconditional_variance(p), g.forecast_variance(p, 5.0e-4, 400), 1.0e-8
  end

  def test_gjr_gamma_zero_reduces_identically_to_garch11
    g = GarchFamily::Garch11.new
    gjr = GarchFamily::GjrGarch11.new
    pg = { mu: 0.0002, omega: 5.0e-6, alpha: 0.08, beta: 0.85 }
    pj = pg.merge(gamma: 0.0)
    rng = GarchFamily::Rng::Sampler.new(101)
    rets = g.simulate(pg, 500, GarchFamily::Innovations::Normal.new, rng)[:returns]
    assert_in_delta g.unconditional_variance(pg), gjr.unconditional_variance(pj), 1.0e-15
    assert_in_delta g.loglik(pg, rets, GarchFamily::Innovations::Normal.new),
                    gjr.loglik(pj, rets, GarchFamily::Innovations::Normal.new), 1.0e-12
  end

  def test_egarch_gamma_zero_has_symmetric_news_impact_curve
    eg = GarchFamily::Egarch11.new
    p = { mu: 0.0, omega: -0.35, alpha: 0.15, gamma: 0.0, beta: 0.95 }
    nic = eg.news_impact(p, [-0.03, -0.01, 0.01, 0.03])
    assert_in_delta nic[0], nic[3], 1.0e-15
    assert_in_delta nic[1], nic[2], 1.0e-15
  end

  def test_analytic_scores_match_central_finite_differences_for_all_models
    rng = GarchFamily::Rng::Sampler.new(202)
    dist = GarchFamily::Innovations::Normal.new
    g = GarchFamily::Garch11.new
    rets = g.simulate({ mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85 }, 600, dist, rng)[:returns]
    check_score_fd(g, { mu: 0.0001, omega: 6.0e-6, alpha: 0.07, beta: 0.84 }, rets, dist)
    check_score_fd(GarchFamily::GjrGarch11.new,
                   { mu: 0.0001, omega: 6.0e-6, alpha: 0.05, gamma: 0.06, beta: 0.84 }, rets, dist)
    check_score_fd(GarchFamily::Egarch11.new,
                   { mu: 0.0001, omega: -0.4, alpha: 0.14, gamma: -0.07, beta: 0.93 }, rets, dist)
  end

  def test_analytic_score_matches_fd_with_student_t_innovation
    rng = GarchFamily::Rng::Sampler.new(303)
    dist = GarchFamily::Innovations::StudentT.new(nu: 6.0)
    rets = GarchFamily::Garch11.new.simulate({ mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85, nu: 6.0 },
                                             500, dist, rng)[:returns]
    check_score_fd(GarchFamily::Garch11.new, { mu: 0.0, omega: 6.0e-6, alpha: 0.08, beta: 0.84, nu: 5.5 }, rets, dist)
  end

  def test_ewma_filter_matches_recursive_definition
    ewma = GarchFamily::Ewma.new(lambda_decay: 0.94)
    rets = Array.new(30) { |i| (i.even? ? 1 : -1) * 0.01 }
    vars = ewma.filter({}, rets)
    assert_in_delta (0.94 * vars[0]) + (0.06 * rets[0] * rets[0]), vars[1], 1.0e-15
  end

  def test_filter_matches_simulation_true_variances
    g = GarchFamily::Garch11.new
    p = { mu: 0.0003, omega: 5.0e-6, alpha: 0.08, beta: 0.86 }
    rng = GarchFamily::Rng::Sampler.new(404)
    sim = g.simulate(p, 400, GarchFamily::Innovations::Normal.new, rng)
    filt = g.filter(p, sim[:returns])
    assert_in_delta sim[:true_vars].last, filt.last, 1.0e-12
  end

  private

  def check_score_fd(model, params, rets, dist)
    an = model.score(params, rets, dist)
    model.param_names(dist).each_with_index do |k, idx|
      eps = k == :omega ? 1.0e-9 : 1.0e-6
      fd = (model.loglik(params.merge(k => params[k] + eps), rets, dist) -
            model.loglik(params.merge(k => params[k] - eps), rets, dist)) / (2.0 * eps)
      rel = (an[idx] - fd).abs / [fd.abs, 1.0].max
      assert_operator rel, :<, 1.0e-5, "#{model.name} #{k}: an=#{an[idx]} fd=#{fd}"
    end
  end
end
