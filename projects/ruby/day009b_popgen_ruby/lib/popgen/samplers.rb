# frozen_string_literal: true

require_relative 'binomial'

module PopGen
  # Exact discrete and continuous samplers using stdlib Random only.
  module Samplers
    module_function

    def multinomial(n, probs, rng)
      k = probs.length
      counts = Array.new(k, 0)
      rem_n = n
      rem_p = 1.0
      (0...(k - 1)).each do |i|
        break if rem_n <= 0

        pi = probs[i].to_f
        p_cond = rem_p > 0.0 ? (pi / rem_p) : 0.0
        p_cond = 1.0 if p_cond > 1.0
        p_cond = 0.0 if p_cond < 0.0
        c = Binomial.sample(rem_n, p_cond, rng)
        counts[i] = c
        rem_n -= c
        rem_p -= pi
      end
      counts[k - 1] = rem_n if k.positive?
      counts
    end

    def exponential(rate, rng)
      return Float::INFINITY if rate <= 0.0

      u = rng.rand
      u = 1e-16 if u <= 0.0
      -Math.log(u) / rate
    end

    # Geometric(p) returning number of trials >= 1 before first success
    def geometric_trials(p, rng)
      return 1 if p >= 1.0

      u = rng.rand
      u = 1e-16 if u <= 0.0
      (Math.log(u) / Math.log(1.0 - p)).floor + 1
    end

    # Poisson via Knuth inversion (lambda < 10) or PTRS (Hormann 1993) exact rejection
    def poisson_knuth(lam, rng)
      l = Math.exp(-lam)
      k = 0
      p = 1.0
      loop do
        k += 1
        p *= rng.rand
        break if p <= l
      end
      k - 1
    end

    def poisson_ptrs(lam, rng)
      b = 0.931 + 2.53 * Math.sqrt(lam)
      a = -0.059 + 0.02483 * b
      invalpha = 1.1239 + 1.1328 / (b - 3.4)
      vr = 0.9277 - 3.6224 / (b - 2.0)
      log_lam = Math.log(lam)

      loop do
        u = rng.rand - 0.5
        v = rng.rand
        us = 0.5 - u.abs
        k = ((2.0 * a / us + b) * u + lam + 0.43).floor
        next if k.negative?

        return k if us >= 0.07 && v <= vr
        next if us < 0.013 && v > us

        lhs = Math.log(v * invalpha / (a / (us * us) + b))
        rhs = -lam + k * log_lam - Math.lgamma(k + 1).first
        return k if lhs <= rhs
      end
    end

    def poisson(lam, rng)
      return 0 if lam <= 0.0

      lam < 10.0 ? poisson_knuth(lam, rng) : poisson_ptrs(lam, rng)
    end
  end
end
