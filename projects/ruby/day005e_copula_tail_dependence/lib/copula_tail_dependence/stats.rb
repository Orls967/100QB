# frozen_string_literal: true

module CopulaTailDependence
  module Stats
    module_function

    def mean(xs)
      return 0.0 if xs.empty?

      xs.sum(0.0) / xs.length
    end

    def variance(xs)
      n = xs.length
      return 0.0 if n < 2

      m = mean(xs)
      ss = xs.sum(0.0) { |x| (x - m) * (x - m) }
      ss / (n - 1)
    end

    def sd(xs)
      Math.sqrt(variance(xs))
    end

    def se(xs)
      return 0.0 if xs.empty?

      sd(xs) / Math.sqrt(xs.length)
    end

    def quantile(xs, q)
      return 0.0 if xs.empty?

      sorted = xs.sort
      pos = q.clamp(0.0, 1.0) * (sorted.length - 1)
      lo = pos.floor
      hi = pos.ceil
      frac = pos - lo
      (sorted[lo] * (1.0 - frac)) + (sorted[hi] * frac)
    end

    # Jackknife SE dari estimator blok (mis. K blok independen atau delete-1 blok).
    def jackknife_se(block_estimates)
      k = block_estimates.length
      return 0.0 if k < 2

      m = mean(block_estimates)
      ss = block_estimates.sum(0.0) { |v| (v - m) * (v - m) }
      Math.sqrt(ss / (k * (k - 1)))
    end
  end
end
