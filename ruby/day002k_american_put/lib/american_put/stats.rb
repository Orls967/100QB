# frozen_string_literal: true

module AmericanPut
  module Stats
    CI_Z95 = 1.96
    INV_SQRT2 = 1.0 / Math.sqrt(2.0)

    module_function

    def norm_cdf(value)
      0.5 * Math.erfc(-value * INV_SQRT2)
    end

    def bs_put(params)
      params.validate!
      spot = params.s0
      strike = params.k
      rate = params.r
      maturity = params.t
      return [(strike * Math.exp(-rate * maturity)) - spot, 0.0].max if maturity <= 0.0 || params.sigma <= 0.0

      vol_sqrt_t = params.sigma * Math.sqrt(maturity)
      d1 = (Math.log(spot / strike) + ((rate + (0.5 * params.sigma * params.sigma)) * maturity)) / vol_sqrt_t
      d2 = d1 - vol_sqrt_t
      (strike * Math.exp(-rate * maturity) * norm_cdf(-d2)) - (spot * norm_cdf(-d1))
    end

    def summarize_samples(samples)
      n = samples.length
      raise ArgumentError, 'jumlah sampel harus >= 2' if n < 2

      mean = 0.0
      m2 = 0.0
      i = 0
      while i < n
        x = samples[i]
        delta = x - mean
        mean += delta / (i + 1)
        m2 += delta * (x - mean)
        i += 1
      end

      variance = [m2 / (n - 1).to_f, 0.0].max
      se = Math.sqrt(variance / n.to_f)
      {
        price: mean,
        se: se,
        ci_low: mean - (CI_Z95 * se),
        ci_high: mean + (CI_Z95 * se)
      }
    end

    def z_score(diff, se_a, se_b = 0.0)
      combined = Math.hypot(se_a, se_b)
      return 0.0 if combined.zero?

      diff.abs / combined
    end

    def fit_loglog_slope(counts, errors)
      raise ArgumentError, 'panjang counts dan errors harus >= 2' if counts.length != errors.length || counts.length < 2

      xs = counts.map { |n| Math.log10(n.to_f) }
      ys = errors.map { |e| Math.log10([e.to_f.abs, 1e-16].max) }
      mx = xs.sum / xs.length
      my = ys.sum / ys.length

      num = 0.0
      den = 0.0
      xs.each_index do |i|
        dx = xs[i] - mx
        num += dx * (ys[i] - my)
        den += dx * dx
      end
      num / den
    end
  end
end
