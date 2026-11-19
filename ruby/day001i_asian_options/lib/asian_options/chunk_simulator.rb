# frozen_string_literal: true

module AsianOptions
  module ChunkSimulator
    module_function

    def run_european(spec, chunk_idx, len, seed)
      run_univariate(spec, chunk_idx, len, seed) do |normals|
        s_t = Math.exp(PathKernels.terminal_log_price(normals, spec))
        s_t > spec.k ? (s_t - spec.k) * spec.disc : 0.0
      end
    end

    def run_naive(spec, chunk_idx, len, seed)
      run_univariate(spec, chunk_idx, len, seed) do |normals|
        avg = PathKernels.arithmetic_average(normals, spec)
        avg > spec.k ? (avg - spec.k) * spec.disc : 0.0
      end
    end

    def run_antithetic(spec, chunk_idx, len, seed)
      run_univariate(spec, chunk_idx, len, seed) do |normals|
        PathKernels.antithetic_pair_payoff(normals, spec)
      end
    end

    def run_control_variate(spec, chunk_idx, len, seed)
      rng = Rng::NormalGenerator.new(Rng.derive_chunk_seed(seed, chunk_idx))
      normals = Array.new(spec.n_steps, 0.0)
      mean_x = 0.0
      mean_y = 0.0
      m2_x = 0.0
      m2_y = 0.0
      c_xy = 0.0
      i = 0
      while i < len
        rng.fill_normals!(normals)
        x, y = PathKernels.cv_path_payoffs(normals, spec)
        i += 1
        dx = x - mean_x
        mean_x += dx / i
        dy = y - mean_y
        mean_y += dy / i
        m2_x += dx * (x - mean_x)
        m2_y += dy * (y - mean_y)
        c_xy += dx * (y - mean_y)
      end
      Stats::BivariateChunk.new(
        sample_count: len, mean_x: mean_x, mean_y: mean_y, m2_x: m2_x, m2_y: m2_y, c_xy: c_xy
      )
    end

    def run_univariate(spec, chunk_idx, len, seed)
      rng = Rng::NormalGenerator.new(Rng.derive_chunk_seed(seed, chunk_idx))
      normals = Array.new(spec.n_steps, 0.0)
      mean = 0.0
      m2 = 0.0
      i = 0
      while i < len
        rng.fill_normals!(normals)
        val = yield(normals)
        i += 1
        delta = val - mean
        mean += delta / i
        m2 += delta * (val - mean)
      end
      Stats::UnivariateChunk.new(sample_count: len, mean: mean, m2: m2)
    end
  end
end
