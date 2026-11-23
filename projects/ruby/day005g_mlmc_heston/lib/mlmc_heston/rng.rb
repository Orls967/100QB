# frozen_string_literal: true

module MlmcHeston
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
    end
  end
end
