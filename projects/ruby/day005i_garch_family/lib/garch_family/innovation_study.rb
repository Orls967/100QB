# frozen_string_literal: true

module GarchFamily
  module InnovationStudy
    TRUE_T_PARAMS = { mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85, nu: 5.0 }.freeze

    module_function

    def run_study(seed: MAIN_SEED, workers: nil, n_reps: 300)
      eff_rows = [1000, 2500].each_with_index.map do |n_obs, idx|
        run_efficiency_cell(n_obs: n_obs, n_reps: n_reps, stream_id: 300 + idx, seed: seed, workers: workers)
      end
      lr_null = run_lr_boundary_cell(n_obs: 1500, n_reps: n_reps, dgp_dist: Innovations::Normal.new,
                                     stream_id: 320, seed: seed, workers: workers)
      lr_alt = run_lr_boundary_cell(n_obs: 1500, n_reps: n_reps, dgp_dist: Innovations::StudentT.new(nu: 5.0),
                                    stream_id: 321, seed: seed, workers: workers)
      { efficiency_comparison: eff_rows, lr_boundary_null_normal: lr_null, lr_boundary_alt_t5: lr_alt }
    end

    def run_efficiency_cell(n_obs:, n_reps:, stream_id:, seed:, workers:)
      model = Garch11.new
      t5 = Innovations::StudentT.new(nu: 5.0)
      norm = Innovations::Normal.new
      trials = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 20) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, idx))
        rets = model.simulate(TRUE_T_PARAMS, n_obs, t5, rng)[:returns]
        fq = Optim.fit(model, rets, dist: norm, init_params: Mle::TRUE_GARCH, max_iter: 100)
        ft = Optim.fit(model, rets, dist: t5, init_params: fq[:params].merge(nu: 6.0),
                                    extra_starts: [TRUE_T_PARAMS], max_iter: 120)
        { ok: fq[:converged] && ft[:converged], q: fq[:params], t: ft[:params] }
      end
      summarize_efficiency(trials, n_obs)
    end

    def summarize_efficiency(trials, n_obs)
      good = trials.select { |r| r[:ok] }
      nus = good.map { |r| r[:t][:nu] }.sort
      nu_mean = nus.sum(0.0) / nus.length
      nu_med = nus[nus.length / 2]
      inv_nu = good.map { |r| 1.0 / r[:t][:nu] }
      metrics = %i[omega alpha beta].each_with_object({}) do |k, acc|
        vq = mean_sq_err(good.map { |r| r[:q][k] }, TRUE_T_PARAMS[k])
        vt = mean_sq_err(good.map { |r| r[:t][k] }, TRUE_T_PARAMS[k])
        acc[k] = { rmse_qmle: Math.sqrt(vq), rmse_t_mle: Math.sqrt(vt), mse_eff_ratio: vq / [vt, 1.0e-18].max }
      end
      { n_obs: n_obs, n_reps: trials.length, n_converged: good.length, n_failed: trials.length - good.length,
        nu_mean: nu_mean, nu_median: nu_med, inv_nu_mean: inv_nu.sum(0.0) / inv_nu.length,
        inv_nu_bias: (inv_nu.sum(0.0) / inv_nu.length) - 0.2, params: metrics }
    end
    private_class_method :summarize_efficiency

    def run_lr_boundary_cell(n_obs:, n_reps:, dgp_dist:, stream_id:, seed:, workers: nil)
      model = Garch11.new
      norm = Innovations::Normal.new
      t_fit = Innovations::StudentT.new(nu: 12.0)
      p_dgp = dgp_dist.is_a?(Innovations::StudentT) ? TRUE_T_PARAMS : Mle::TRUE_GARCH
      stats = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 20) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, idx))
        rets = model.simulate(p_dgp, n_obs, dgp_dist, rng)[:returns]
        single_lr_boundary(model, rets, norm, t_fit)
      end
      summarize_lr_boundary(stats, dgp_dist.name, n_obs)
    end

    def single_lr_boundary(model, rets, norm, t_fit)
      f0 = Optim.fit(model, rets, dist: norm, init_params: Mle::TRUE_GARCH, max_iter: 100)
      starts = [f0[:params].merge(nu: 15.0), f0[:params].merge(nu: 80.0)]
      f1 = Optim.fit(model, rets, dist: t_fit, init_params: starts[0], extra_starts: [starts[1]], max_iter: 100)
      lr = [2.0 * (f1[:loglik] - f0[:loglik]), 0.0].max
      { lr: lr, hit_zero: lr < 1.0e-5, p_mix: Special.chi2_mixture_01_sf(lr), p_chi1: Special.chi2_sf(lr, 1) }
    end
    private_class_method :single_lr_boundary

    def summarize_lr_boundary(stats, dgp_name, n_obs)
      n = stats.length.to_f
      rej_mix = stats.count { |s| s[:p_mix] < 0.05 } / n
      rej_chi1 = stats.count { |s| s[:p_chi1] < 0.05 } / n
      zero_frac = stats.count { |s| s[:hit_zero] } / n
      { dgp: dgp_name, n_obs: n_obs, n_reps: stats.length,
        zero_mass_fraction: zero_frac, rej_rate_mixture_5pct: rej_mix,
        rej_rate_se: Math.sqrt(rej_mix * (1.0 - rej_mix) / n),
        rej_rate_naive_chi2_5pct: rej_chi1 }
    end
    private_class_method :summarize_lr_boundary

    def mean_sq_err(vals, target) = vals.sum(0.0) { |v| (v - target)**2 } / vals.length
    private_class_method :mean_sq_err
  end
end
