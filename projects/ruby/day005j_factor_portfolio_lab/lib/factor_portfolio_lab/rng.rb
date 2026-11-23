# frozen_string_literal: true

module FactorPortfolioLab
  module Rng
    MASK64 = (1 << 64) - 1
    GOLDEN_GAMMA = 0x9E3779B97F4A7C15
    MIX_MUL1 = 0xBF58476D1CE4E5B9
    MIX_MUL2 = 0x94D049BB133111EB

    module_function

    def derive_seed(master_seed, stream_id, chunk_index)
      step = ((stream_id.to_i + 1) * 0x517CC1B727220A95) + ((chunk_index.to_i + 1) * GOLDEN_GAMMA)
      z = (master_seed.to_i + step) & MASK64
      z = ((z ^ (z >> 30)) * MIX_MUL1) & MASK64
      z = ((z ^ (z >> 27)) * MIX_MUL2) & MASK64
      (z ^ (z >> 31)) & MASK64
    end

    class Sampler
      def initialize(seed)
        @rng = Random.new(seed)
        @has_spare = false
        @spare = 0.0
      end

      def next_float = @rng.rand

      def next_int(limit) = @rng.rand(limit)

      def next_normal
        if @has_spare
          @has_spare = false
          return @spare
        end

        rng = @rng
        loop do
          u = (2.0 * rng.rand) - 1.0
          v = (2.0 * rng.rand) - 1.0
          s = (u * u) + (v * v)
          next if s >= 1.0 || s <= 0.0

          mul = Math.sqrt(-2.0 * Math.log(s) / s)
          @spare = v * mul
          @has_spare = true
          return u * mul
        end
      end

      def next_gamma(shape)
        if shape < 1.0
          u = [next_float, 1.0e-15].max
          return next_gamma(shape + 1.0) * (u**(1.0 / shape))
        end

        d = shape - (1.0 / 3.0)
        c = 1.0 / Math.sqrt(9.0 * d)
        loop do
          x, v = draw_gamma_candidate(c)
          u = [next_float, 1.0e-15].max
          return d * v if u < 1.0 - (0.0331 * (x**4)) || Math.log(u) < (0.5 * x * x) + (d * (1.0 - v + Math.log(v)))
        end
      end

      def next_chi2(df) = 2.0 * next_gamma(0.5 * df)

      def next_std_t(df)
        z = next_normal
        g = next_gamma(0.5 * df)
        z / Math.sqrt((2.0 * g) / (df - 2.0))
      end

      def next_wishart(scale_mat, df)
        n = scale_mat.length
        chol = LinAlg.cholesky(scale_mat)
        a_mat = Array.new(n) { Array.new(n, 0.0) }
        n.times do |i|
          a_mat[i][i] = Math.sqrt(next_chi2(df - i))
          (0...i).each { |j| a_mat[i][j] = next_normal }
        end
        la = LinAlg.mat_mul(chol, a_mat)
        LinAlg.mat_mul(la, LinAlg.transpose(la))
      end

      private

      def draw_gamma_candidate(c)
        loop do
          x = next_normal
          v = 1.0 + (c * x)
          return [x, v * v * v] if v.positive?
        end
      end
    end
  end
end
