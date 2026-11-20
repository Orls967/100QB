# frozen_string_literal: true

module AsianOptions
  module Estimators
    GbmSpec = Data.define(:s0, :k, :n_steps, :drift, :vol, :disc, :log_s0, :inv_n)

    module_function

    def build_gbm_spec(params)
      dt = params.t.to_f / params.n_steps
      GbmSpec.new(
        s0: params.s0.to_f,
        k: params.k.to_f,
        n_steps: params.n_steps.to_i,
        drift: (params.r - (0.5 * params.sigma * params.sigma)) * dt,
        vol: params.sigma * Math.sqrt(dt),
        disc: Math.exp(-params.r * params.t),
        log_s0: Math.log(params.s0),
        inv_n: 1.0 / params.n_steps
      )
    end

    def validate_mc!(params, n_paths)
      params.validate!
      raise ArgumentError, 'n_paths harus >= 2' if n_paths < 2
    end

    def price_european_mc(params, n_paths:, seed: DEFAULT_SEED, workers: nil)
      validate_mc!(params, n_paths)
      spec = build_gbm_spec(params)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers) do |chunk_idx, len|
        ChunkSimulator.run_european(spec, chunk_idx, len, seed)
      end
      Stats.summarize_univariate(Stats.merge_univariate(chunks))
    end

    def price_asian_naive(params, n_paths:, seed: DEFAULT_SEED, workers: nil)
      validate_mc!(params, n_paths)
      spec = build_gbm_spec(params)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers) do |chunk_idx, len|
        ChunkSimulator.run_naive(spec, chunk_idx, len, seed)
      end
      Stats.summarize_univariate(Stats.merge_univariate(chunks))
    end

    def price_asian_antithetic(params, n_paths:, seed: DEFAULT_SEED, workers: nil)
      params.validate!
      raise ArgumentError, 'n_paths harus genap dan >= 4 untuk antithetic' if n_paths < 4 || n_paths.odd?

      spec = build_gbm_spec(params)
      chunks = Parallel.map_chunks(n_items: n_paths / 2, workers: workers) do |chunk_idx, len|
        ChunkSimulator.run_antithetic(spec, chunk_idx, len, seed)
      end
      Stats.summarize_univariate(Stats.merge_univariate(chunks))
    end

    def price_asian_control_variate(params, n_paths:, seed: DEFAULT_SEED, workers: nil)
      validate_mc!(params, n_paths)
      y_exact = BlackScholes.geometric_asian_call(params)
      spec = build_gbm_spec(params)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers) do |chunk_idx, len|
        ChunkSimulator.run_control_variate(spec, chunk_idx, len, seed)
      end
      finalize_control_variate(Stats.merge_bivariate(chunks), y_exact)
    end

    def finalize_control_variate(merged, y_exact)
      n = merged.sample_count.to_f
      denom = n - 1.0
      var_x = merged.m2_x / denom
      var_y = merged.m2_y / denom
      raise ArgumentError, 'varians control variate bernilai nol' unless var_y.positive?

      cov_xy = merged.c_xy / denom
      beta = cov_xy / var_y
      z_mean = merged.mean_x - (beta * (merged.mean_y - y_exact))
      var_z = [var_x + (beta * beta * var_y) - (2.0 * beta * cov_xy), 0.0].max
      se = Math.sqrt(var_z / n)

      ControlVariateResult.new(
        price: z_mean, se: se, ci_low: z_mean - (Stats::CI_Z95 * se),
        ci_high: z_mean + (Stats::CI_Z95 * se), beta: beta, y_exact: y_exact
      )
    end
    private_class_method :finalize_control_variate
  end
end
