# frozen_string_literal: true

module AsianOptions
  module Greeks
    BumpSpec = Data.define(:scale_up, :scale_dn, :factor, :strike)

    module_function

    def delta_pathwise(params, n_paths:, seed: DEFAULT_SEED, workers: nil)
      Estimators.validate_mc!(params, n_paths)
      spec = Estimators.build_gbm_spec(params)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers) do |chunk_idx, len|
        run_pathwise_chunk(spec, chunk_idx, len, seed)
      end
      to_delta_result(Stats.summarize_univariate(Stats.merge_univariate(chunks)))
    end

    def delta_bump(params, n_paths:, h_rel: 0.01, seed: DEFAULT_SEED, workers: nil)
      Estimators.validate_mc!(params, n_paths)
      raise ArgumentError, 'h_rel harus > 0' unless h_rel.positive?

      spec = Estimators.build_gbm_spec(params)
      bump = build_bump_spec(spec, h_rel.to_f)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers) do |chunk_idx, len|
        run_bump_chunk(spec, bump, chunk_idx, len, seed)
      end
      to_delta_result(Stats.summarize_univariate(Stats.merge_univariate(chunks)))
    end

    def build_bump_spec(spec, h_rel)
      h = h_rel * spec.s0
      BumpSpec.new(
        scale_up: (spec.s0 + h) / spec.s0,
        scale_dn: (spec.s0 - h) / spec.s0,
        factor: spec.disc / (2.0 * h),
        strike: spec.k
      )
    end
    private_class_method :build_bump_spec

    def to_delta_result(est)
      DeltaResult.new(delta: est.price, se: est.se, ci_low: est.ci_low, ci_high: est.ci_high)
    end
    private_class_method :to_delta_result

    def run_pathwise_chunk(spec, chunk_idx, len, seed)
      rng = Rng::NormalGenerator.new(Rng.derive_chunk_seed(seed, chunk_idx))
      normals = Array.new(spec.n_steps, 0.0)
      inv_s0 = 1.0 / spec.s0
      mean = 0.0
      m2 = 0.0
      i = 0
      while i < len
        rng.fill_normals!(normals)
        avg = PathKernels.arithmetic_average(normals, spec)
        val = avg > spec.k ? spec.disc * avg * inv_s0 : 0.0
        i += 1
        delta = val - mean
        mean += delta / i
        m2 += delta * (val - mean)
      end
      Stats::UnivariateChunk.new(sample_count: len, mean: mean, m2: m2)
    end
    private_class_method :run_pathwise_chunk

    def run_bump_chunk(spec, bump, chunk_idx, len, seed)
      rng = Rng::NormalGenerator.new(Rng.derive_chunk_seed(seed, chunk_idx))
      normals = Array.new(spec.n_steps, 0.0)
      mean = 0.0
      m2 = 0.0
      i = 0
      while i < len
        rng.fill_normals!(normals)
        avg = PathKernels.arithmetic_average(normals, spec)
        val = bump_sample(avg, bump)
        i += 1
        delta = val - mean
        mean += delta / i
        m2 += delta * (val - mean)
      end
      Stats::UnivariateChunk.new(sample_count: len, mean: mean, m2: m2)
    end
    private_class_method :run_bump_chunk

    def bump_sample(avg, bump)
      up = avg * bump.scale_up
      dn = avg * bump.scale_dn
      pay_up = up > bump.strike ? up - bump.strike : 0.0
      pay_dn = dn > bump.strike ? dn - bump.strike : 0.0
      (pay_up - pay_dn) * bump.factor
    end
    private_class_method :bump_sample
  end
end
