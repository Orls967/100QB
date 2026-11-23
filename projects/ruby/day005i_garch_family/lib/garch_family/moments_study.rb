# frozen_string_literal: true

module GarchFamily
  module MomentsStudy
    module_function

    GARCH_GRID = [
      { mu: 0.0, omega: 2.0e-6, alpha: 0.05, beta: 0.90 },
      { mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85 },
      { mu: 0.0, omega: 1.0e-5, alpha: 0.10, beta: 0.80 },
      { mu: 0.0, omega: 4.0e-6, alpha: 0.06, beta: 0.88 },
      { mu: 0.0, omega: 8.0e-6, alpha: 0.12, beta: 0.78 }
    ].freeze

    GJR_GRID = [
      { mu: 0.0, omega: 5.0e-6, alpha: 0.04, gamma: 0.06, beta: 0.88 },
      { mu: 0.0, omega: 5.0e-6, alpha: 0.05, gamma: 0.08, beta: 0.85 },
      { mu: 0.0, omega: 8.0e-6, alpha: 0.03, gamma: 0.10, beta: 0.85 },
      { mu: 0.0, omega: 4.0e-6, alpha: 0.06, gamma: 0.04, beta: 0.88 },
      { mu: 0.0, omega: 1.0e-5, alpha: 0.05, gamma: 0.12, beta: 0.80 }
    ].freeze

    EGARCH_GRID = [
      { mu: 0.0, omega: -0.28, alpha: 0.14, gamma: -0.08, beta: 0.97 },
      { mu: 0.0, omega: -0.45, alpha: 0.16, gamma: -0.10, beta: 0.95 },
      { mu: 0.0, omega: -0.60, alpha: 0.18, gamma: -0.05, beta: 0.93 },
      { mu: 0.0, omega: -0.20, alpha: 0.12, gamma: 0.0, beta: 0.97 },
      { mu: 0.0, omega: -0.80, alpha: 0.20, gamma: -0.12, beta: 0.90 }
    ].freeze

    def run_study(seed: MAIN_SEED, workers: nil)
      garch_rows = run_grid(Garch11.new, GARCH_GRID, n_obs: 400_000, stream_base: 100, seed: seed, workers: workers)
      gjr_rows = run_grid(GjrGarch11.new, GJR_GRID, n_obs: 200_000, stream_base: 120, seed: seed, workers: workers)
      eg_rows = run_grid(Egarch11.new, EGARCH_GRID, n_obs: 200_000, stream_base: 140, seed: seed, workers: workers)
      base_p = { mu: 0.0, omega: 1.0e-5, alpha: 0.10, beta: 0.80 }
      acf_res = verify_garch_acf_and_forecast(base_p, seed: seed)
      { garch_grid: garch_rows, gjr_grid: gjr_rows, egarch_grid: eg_rows }.merge(acf_res)
    end

    def run_grid(model, grid, n_obs:, stream_base:, seed:, workers:)
      dist = Innovations::Normal.new
      Parallel.map_chunks(n_items: grid.length, workers: workers, chunk_size: 1) do |idx, _len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_base + idx, 0))
        p = grid[idx]
        rets = model.simulate(p, n_obs, dist, rng)[:returns]
        single_moment_check(model, p, rets)
      end
    end

    def single_moment_check(model, params, rets)
      v_emp, se_v = batch_means_se(rets.map { |r| r * r }, 40)
      v_theo = model.unconditional_variance(params)
      res = { params: params, persistence: model.persistence(params), var_theo: v_theo,
              var_emp: v_emp, var_se: se_v, var_z: (v_emp - v_theo).abs / [se_v, 1.0e-18].max }
      if model.is_a?(Garch11)
        res[:kurt_theo] = model.theoretical_kurtosis(params)
        res[:kurt_emp], res[:kurt_se] = batch_kurtosis_se(rets, 40)
        res[:kurt_z] = (res[:kurt_emp] - res[:kurt_theo]).abs / [res[:kurt_se], 1.0e-12].max
      end
      res
    end

    def batch_means_se(vals, n_batches)
      sz = vals.length / n_batches
      bmeans = Array.new(n_batches) { |b| vals[b * sz, sz].sum(0.0) / sz }
      m = bmeans.sum(0.0) / n_batches
      sd = Math.sqrt(bmeans.sum(0.0) { |x| (x - m)**2 } / (n_batches - 1))
      [m, sd / Math.sqrt(n_batches)]
    end

    def batch_kurtosis_se(rets, n_batches)
      sz = rets.length / n_batches
      bks = Array.new(n_batches) do |b|
        sl = rets[b * sz, sz]
        m2 = sl.sum(0.0) { |x| x * x } / sz
        m4 = sl.sum(0.0) { |x| x**4 } / sz
        m4 / (m2 * m2)
      end
      m = bks.sum(0.0) / n_batches
      sd = Math.sqrt(bks.sum(0.0) { |x| (x - m)**2 } / (n_batches - 1))
      [m, sd / Math.sqrt(n_batches)]
    end

    def verify_garch_acf_and_forecast(params, seed:)
      g = Garch11.new
      dist = Innovations::Normal.new
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 170, 0))
      rets = g.simulate(params, 400_000, dist, rng)[:returns]
      sq = rets.map { |r| r * r }
      acf_rows = (1..5).map do |k|
        emp, se = sample_lag_corr_se(sq, k, 40)
        theo = g.theoretical_acf(params, k)
        { lag: k, theo: theo, emp: emp, se: se, z: (emp - theo).abs / [se, 1.0e-12].max }
      end
      fc = verify_forecast_term_structure(g, params, dist, seed)
      { acf_sq_lags: acf_rows, forecast_mc: fc }
    end

    def sample_lag_corr_se(sq, lag, n_batches)
      sz = sq.length / n_batches
      corrs = Array.new(n_batches) do |b|
        sl = sq[b * sz, sz]
        Special.pearson_corr(sl[lag..], sl[0...-lag])
      end
      m = corrs.sum(0.0) / n_batches
      sd = Math.sqrt(corrs.sum(0.0) { |x| (x - m)**2 } / (n_batches - 1))
      [m, sd / Math.sqrt(n_batches)]
    end

    def verify_forecast_term_structure(garch, params, dist, seed)
      horizons = [1, 5, 10, 20, 60]
      h1 = 2.5 * garch.unconditional_variance(params)
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 190, 0))
      n_paths = 15_000
      sums = Array.new(horizons.length, 0.0)
      n_paths.times do
        cur = Array.new(horizons.last, h1)
        simulate_forward_vars!(params, cur, dist, rng)
        horizons.each_with_index { |hz, idx| sums[idx] += cur[hz - 1] }
      end
      horizons.each_with_index.map do |hz, idx|
        { horizon: hz, mc_var: sums[idx] / n_paths, theo_var: garch.forecast_variance(params, h1, hz) }
      end
    end

    def simulate_forward_vars!(params, cur, dist, rng)
      w = params[:omega]
      a = params[:alpha]
      b = params[:beta]
      (1...cur.length).each do |i|
        z = dist.sample(rng)
        ep2 = cur[i - 1] * z * z
        cur[i] = w + (a * ep2) + (b * cur[i - 1])
      end
    end
    private_class_method :simulate_forward_vars!
  end
end
