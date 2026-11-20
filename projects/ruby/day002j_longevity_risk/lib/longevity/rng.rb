# frozen_string_literal: true

module Longevity
  class Rng
    MASK64 = (1 << 64) - 1

    def self.derive_seed(base_seed, index)
      z = (base_seed + 0x9e3779b97f4a7c15 + (index * 0xbf58476d1ce4e5b9)) & MASK64
      z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
      z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & MASK64
      z ^ (z >> 31)
    end

    def initialize(seed)
      @prng = Random.new(seed)
      @has_spare = false
      @spare = 0.0
    end

    def next_uniform
      @prng.rand
    end

    def next_normal
      if @has_spare
        @has_spare = false
        return @spare
      end

      u1 = 0.0
      u1 = @prng.rand while u1 <= 1e-16
      u2 = @prng.rand
      mag = Math.sqrt(-2.0 * Math.log(u1))
      theta = 2.0 * Math::PI * u2
      @spare = mag * Math.sin(theta)
      @has_spare = true
      mag * Math.cos(theta)
    end
  end
end
