# frozen_string_literal: true

module AmericanPut
  module Rng
    MASK64 = (1 << 64) - 1
    GOLDEN_GAMMA = 0x9E3779B97F4A7C15
    MIX_MUL1 = 0xBF58476D1CE4E5B9
    MIX_MUL2 = 0x94D049BB133111EB

    def self.derive_chunk_seed(master_seed, chunk_index)
      z = (master_seed.to_i + ((chunk_index + 1) * GOLDEN_GAMMA)) & MASK64
      z = ((z ^ (z >> 30)) * MIX_MUL1) & MASK64
      z = ((z ^ (z >> 27)) * MIX_MUL2) & MASK64
      (z ^ (z >> 31)) & MASK64
    end

    class NormalGenerator
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

        z0, z1 = draw_polar_pair
        @spare = z1
        @has_spare = true
        z0
      end

      # Mengisi buffer Float pra-alokasi dengan sampel N(0, 1) memakai Box-Muller polar in-place.
      def fill_normals!(buffer)
        rng = @rng
        len = buffer.length
        i = 0

        if @has_spare && i < len
          buffer[0] = @spare
          @has_spare = false
          i = 1
        end

        while i < len
          u = (2.0 * rng.rand) - 1.0
          v = (2.0 * rng.rand) - 1.0
          s = (u * u) + (v * v)
          next if s >= 1.0 || s.zero?

          mul = Math.sqrt(-2.0 * Math.log(s) / s)
          buffer[i] = u * mul
          i += 1
          if i < len
            buffer[i] = v * mul
            i += 1
          else
            @spare = v * mul
            @has_spare = true
          end
        end

        buffer
      end

      private

      def draw_polar_pair
        rng = @rng
        loop do
          u = (2.0 * rng.rand) - 1.0
          v = (2.0 * rng.rand) - 1.0
          s = (u * u) + (v * v)
          next if s >= 1.0 || s.zero?

          mul = Math.sqrt(-2.0 * Math.log(s) / s)
          return [u * mul, v * mul]
        end
      end
    end
  end
end
