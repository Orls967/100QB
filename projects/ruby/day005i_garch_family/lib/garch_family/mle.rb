# frozen_string_literal: true

module GarchFamily
  module Mle
    TRUE_GARCH = { mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85 }.freeze

    module_function

    def run_study(seed: MAIN_SEED, workers: nil, n_reps: 1000)
      t_grid = [500, 1000, 2500]
      norm_dist = Innovations::Normal.new
      t5_dist = Innovations::StudentT.new(nu: 5.0)
      norm_rows = t_grid.each_with_index.map do |n_obs, idx|
        run_qmle_cell(n_obs: n_obs, n_reps: n_reps, true_dist: norm_dist, stream_id: 200 + idx, seed: seed,
                      workers: workers)
      end
      t5_rows = t_grid.each_with_index.map do |n_obs, idx|
        run_qmle_cell(n_obs: n_obs, n_reps: n_reps, true_dist: t5_dist, stream_id: 220 + idx, seed: seed,
                      workers: workers)
      end
      anchor = single_sample_score_anchor(seed)
      { true_params: TRUE_GARCH, normal_qmle: norm_rows, student_t5_qmle: t5_rows,
        single_sample_anchor: anchor }
    end

    def run_qmle_cell(n_obs:, n_reps:, true_dist:, stream_id:, seed:, workers:)
      model = Garch11.new
      fit_dist = Innovations::Normal.new
      trials = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 25) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, idx))
        rets = model.simulate(TRUE_GARCH, n_obs, true_dist, rng)[:returns]
        evaluate_single_qmle(model, rets, fit_dist)
      end
      summarize_qmle_trials(trials, n_obs, true_dist.name)
    end

    def evaluate_single_qmle(model, rets, fit_dist)
      extra = [{ mu: 0.0, omega: 4.0e-6, alpha: 0.10, beta: 0.82 }]
      fit = Optim.fit(model, rets, dist: fit_dist, init_params: TRUE_GARCH, extra_starts: extra, max_iter: 120)
      return { ok: false, grad_inf: fit[:grad_inf] } unless fit[:converged]

      rep = Inference.covariance_report(model, fit[:params], rets, fit_dist)
      return { ok: false, grad_inf: fit[:grad_inf] } if rep.nil?

      est = model.param_names(fit_dist).map { |k| fit[:params][k] }
      { ok: true, est: est, se_naive: rep[:se_naive], se_sw: rep[:se_sandwich],
        grad_inf: fit[:grad_inf], loglik: fit[:loglik] }
    end
    private_class_method :evaluate_single_qmle

    def summarize_qmle_trials(trials, n_obs, dist_name)
      good = trials.select { |t| t[:ok] }
      failed = trials.length - good.length
      keys = %i[mu omega alpha beta]
      truth = keys.map { |k| TRUE_GARCH[k] }
      stats = keys.each_with_index.with_object({}) do |(k, idx), acc|
        acc[k] = param_summary(good, idx, truth[idx])
      end
      { n_obs: n_obs, dgp_dist: dist_name, n_reps: trials.length,
        n_converged: good.length, n_failed: failed, params: stats }
    end
    private_class_method :summarize_qmle_trials

    def param_summary(good, idx, true_val)
      n = good.length.to_f
      ests = good.map { |r| r[:est][idx] }
      mean = ests.sum(0.0) / n
      bias = mean - true_val
      rmse = Math.sqrt(ests.sum(0.0) { |v| (v - true_val)**2 } / n)
      emp_sd = Math.sqrt(ests.sum(0.0) { |v| (v - mean)**2 } / [n - 1.0, 1.0].max)
      se_nv = good.sum(0.0) { |r| r[:se_naive][idx] } / n
      se_sw = good.sum(0.0) { |r| r[:se_sw][idx] } / n
      cov_nv = coverage_rate(good, idx, true_val, :se_naive)
      cov_sw = coverage_rate(good, idx, true_val, :se_sw)
      { truth: true_val, mean: mean, bias: bias, rmse: rmse, emp_sd: emp_sd,
        mean_se_naive: se_nv, mean_se_sandwich: se_sw,
        coverage_naive_95: cov_nv, coverage_sandwich_95: cov_sw,
        coverage_se: Math.sqrt(cov_sw * (1.0 - cov_sw) / n) }
    end
    private_class_method :param_summary

    def coverage_rate(good, idx, true_val, key)
      hits = good.count do |r|
        se = r[key][idx]
        se.finite? && (r[:est][idx] - true_val).abs <= 1.959963984540054 * se
      end
      hits.to_f / good.length
    end
    private_class_method :coverage_rate

    def single_sample_score_anchor(seed)
      model = Garch11.new
      dist = Innovations::Normal.new
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 240, 0))
      rets = model.simulate(TRUE_GARCH, 5000, dist, rng)[:returns]
      fit = Optim.fit(model, rets, dist: dist, init_params: TRUE_GARCH, max_iter: 200, gtol: 1.0e-8)
      sc = model.score(fit[:params], rets, dist)
      scales = [Math.sqrt(model.unconditional_variance(fit[:params])), fit[:params][:omega], 1.0, 1.0]
      dimless_sc = sc.each_with_index.map { |v, i| (v * scales[i]).abs / 5000.0 }.max
      { n_obs: 5000, converged: fit[:converged], grad_theta_inf: fit[:grad_inf],
        dimless_score_inf: dimless_sc, loglik: fit[:loglik] }
    end
    private_class_method :single_sample_score_anchor
  end
end
