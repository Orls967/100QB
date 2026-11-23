# frozen_string_literal: true

module CopulaTailDependence
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

      def next_uniform
        @rng.rand
      end

      def next_exp(rate = 1.0)
        u = @rng.rand
        u = @rng.rand while u <= 0.0
        -Math.log(u) / rate
      end

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

      # Marsaglia & Tsang (2000) untuk alpha >= 1, dengan boost untuk alpha < 1.
      def next_gamma(alpha, scale = 1.0)
        raise ArgumentError, 'alpha harus > 0' if alpha <= 0.0

        if alpha < 1.0
          u = @rng.rand
          u = @rng.rand while u <= 0.0
          return next_gamma(alpha + 1.0, scale) * (u**(1.0 / alpha))
        end

        d = alpha - (1.0 / 3.0)
        c = 1.0 / Math.sqrt(9.0 * d)
        loop do
          x = next_normal
          v = 1.0 + (c * x)
          next if v <= 0.0

          v3 = v * v * v
          u = @rng.rand
          return d * v3 * scale if u < 1.0 - (0.0331 * x * x * x * x)
          return d * v3 * scale if Math.log(u) < (0.5 * x * x) + (d * (1.0 - v3 + Math.log(v3)))
        end
      end

      def next_chi2(df)
        2.0 * next_gamma(0.5 * df, 1.0)
      end

      def next_t(df)
        z = next_normal
        chi2 = next_chi2(df)
        z / Math.sqrt(chi2 / df)
      end

      # Chambers-Mallows-Stuck untuk distribusi stable positif dengan transformasi Laplace exp(-s^alpha).
      def next_positive_stable(alpha)
        return 1.0 if (alpha - 1.0).abs < 1e-12

        u = Math::PI * @rng.rand
        u = Math::PI * @rng.rand while u <= 0.0 || u >= Math::PI
        w = next_exp(1.0)

        a = alpha
        sin_au = Math.sin(a * u)
        sin_u = Math.sin(u)
        sin_1au = Math.sin((1.0 - a) * u)
        (sin_au / (sin_u**(1.0 / a))) * ((sin_1au / w)**((1.0 - a) / a))
      end
    end
  end
end
