# frozen_string_literal: true

require 'minitest/autorun'
require 'garch_family'

class TestMleAndInference < Minitest::Test
  def test_bfgs_minimizes_rosenbrock_to_machine_precision
    obj = lambda do |x|
      f = ((1.0 - x[0])**2) + (100.0 * ((x[1] - (x[0] * x[0]))**2))
      g = [(-2.0 * (1.0 - x[0])) - (400.0 * x[0] * (x[1] - (x[0] * x[0]))), 200.0 * (x[1] - (x[0] * x[0]))]
      [f, g]
    end
    res = GarchFamily::Optim.bfgs([-1.2, 1.0], obj, max_iter: 150, gtol: 1.0e-10)
    assert res[:converged]
    assert_in_delta 1.0, res[:x][0], 1.0e-7
    assert_in_delta 1.0, res[:x][1], 1.0e-7
  end

  def test_mle_normal_recovers_garch11_with_zero_score
    rng = GarchFamily::Rng::Sampler.new(1001)
    g = GarchFamily::Garch11.new
    dist = GarchFamily::Innovations::Normal.new
    p_true = { mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85 }
    rets = g.simulate(p_true, 2000, dist, rng)[:returns]
    fit = GarchFamily::Optim.fit(g, rets, dist: dist)
    assert fit[:converged]
    assert_operator fit[:avg_grad_inf], :<, 1.0e-4
    assert_in_delta 0.93, g.persistence(fit[:params]), 0.06
  end

  def test_mle_also_recovers_gjr_and_egarch_params
    rng = GarchFamily::Rng::Sampler.new(1002)
    dist = GarchFamily::Innovations::Normal.new
    gjr = GarchFamily::GjrGarch11.new
    r_gjr = gjr.simulate({ mu: 0.0, omega: 5.0e-6, alpha: 0.04, gamma: 0.08, beta: 0.85 }, 2500, dist, rng)[:returns]
    fit_gjr = GarchFamily::Optim.fit(gjr, r_gjr, dist: dist)
    assert fit_gjr[:converged]
    assert_operator fit_gjr[:params][:gamma], :>, 0.02

    eg = GarchFamily::Egarch11.new
    r_eg = eg.simulate({ mu: 0.0, omega: -0.40, alpha: 0.14, gamma: -0.09, beta: 0.95 }, 2500, dist, rng)[:returns]
    fit_eg = GarchFamily::Optim.fit(eg, r_eg, dist: dist)
    assert fit_eg[:converged]
    assert_operator fit_eg[:params][:gamma], :<, -0.02
  end

  def test_hessian_matches_opg_ratio_near_one_on_correct_normal_dgp
    rng = GarchFamily::Rng::Sampler.new(1003)
    g = GarchFamily::Garch11.new
    dist = GarchFamily::Innovations::Normal.new
    rets = g.simulate(GarchFamily::Mle::TRUE_GARCH, 20_000, dist, rng)[:returns]
    fit = GarchFamily::Optim.fit(g, rets, dist: dist)
    rep = GarchFamily::Inference.covariance_report(g, fit[:params], rets, dist)
    (0..3).each { |i| assert_in_delta 1.0, rep[:se_sandwich][i] / rep[:se_naive][i], 0.15 }
  end

  def test_sandwich_se_exceeds_naive_se_under_student_t_qmle
    rng = GarchFamily::Rng::Sampler.new(1004)
    g = GarchFamily::Garch11.new
    t5 = GarchFamily::Innovations::StudentT.new(nu: 5.0)
    rets = g.simulate(GarchFamily::InnovationStudy::TRUE_T_PARAMS, 4000, t5, rng)[:returns]
    norm = GarchFamily::Innovations::Normal.new
    fit = GarchFamily::Optim.fit(g, rets, dist: norm)
    rep = GarchFamily::Inference.covariance_report(g, fit[:params], rets, norm)
    assert_operator rep[:se_sandwich][2], :>, rep[:se_naive][2] * 1.15
  end

  def test_student_t_mle_recovers_nu_and_beats_normal_loglik
    rng = GarchFamily::Rng::Sampler.new(1005)
    g = GarchFamily::Garch11.new
    t5 = GarchFamily::Innovations::StudentT.new(nu: 5.0)
    rets = g.simulate(GarchFamily::InnovationStudy::TRUE_T_PARAMS, 2500, t5, rng)[:returns]
    fn = GarchFamily::Optim.fit(g, rets, dist: GarchFamily::Innovations::Normal.new)
    ft = GarchFamily::Optim.fit(g, rets, dist: t5, init_params: fn[:params].merge(nu: 6.0))
    assert ft[:converged]
    assert_in_delta 5.0, ft[:params][:nu], 2.0
    assert_operator ft[:loglik] - fn[:loglik], :>, 20.0
  end

  def test_wald_and_lr_test_calculations
    w = GarchFamily::Inference.wald_test([0.08, 0.0], [[0.0004, 0.0], [0.0, 1.0]], [0], [0.0])
    assert_in_delta 16.0, w[:stat], 1.0e-12
    assert_operator w[:p_value], :<, 1.0e-4
    lr = GarchFamily::Inference.lr_test(-100.0, -101.35277, mixture_zero_one: true)
    assert_in_delta 0.05, lr[:p_value], 1.0e-4
  end

  def test_optim_initial_params_covers_all_models_and_distributions
    rng = GarchFamily::Rng::Sampler.new(1008)
    rets = Array.new(200) { 0.01 * rng.next_normal }
    t5 = GarchFamily::Innovations::StudentT.new(nu: 5.0)
    [GarchFamily::Garch11.new, GarchFamily::GjrGarch11.new, GarchFamily::Egarch11.new].each do |m|
      pn = GarchFamily::Optim.initial_params(m, rets, GarchFamily::Innovations::Normal.new)
      pt = GarchFamily::Optim.initial_params(m, rets, t5)
      assert_equal m.param_names(GarchFamily::Innovations::Normal.new).sort, pn.keys.sort
      assert_equal m.param_names(t5).sort, pt.keys.sort
    end
  end

  def test_mle_run_qmle_cell_smoke_summary
    res = GarchFamily::Mle.run_qmle_cell(
      n_obs: 600, n_reps: 8, true_dist: GarchFamily::Innovations::Normal.new,
      stream_id: 900, seed: 42, workers: 2
    )
    assert_operator res[:n_converged], :>=, 6
    assert_operator res[:params][:alpha][:mean_se_sandwich], :>, 0.0
  end
end
