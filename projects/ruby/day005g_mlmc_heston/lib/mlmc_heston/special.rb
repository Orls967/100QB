# frozen_string_literal: true

module MlmcHeston
  module Special
    SQRT_2 = Math.sqrt(2.0)

    module_function

    def norm_cdf(x) = 0.5 * (1.0 + Math.erf(x / SQRT_2))

    def cexp(z)
      er = Math.exp(z.real)
      Complex(er * Math.cos(z.imag), er * Math.sin(z.imag))
    end

    def clog(z)
      Complex(Math.log(z.abs), z.arg)
    end

    def csqrt(z)
      r = z.abs
      return Complex(0.0, 0.0) if r.zero?

      ang = 0.5 * z.arg
      sr = Math.sqrt(r)
      Complex(sr * Math.cos(ang), sr * Math.sin(ang))
    end

    def bs_call(spot, strike, rate, sigma, maturity)
      return [spot - strike, 0.0].max if maturity <= 0.0 || sigma <= 0.0

      vol_t = sigma * Math.sqrt(maturity)
      d1 = (Math.log(spot / strike) + ((rate + (0.5 * sigma * sigma)) * maturity)) / vol_t
      d2 = d1 - vol_t
      (spot * norm_cdf(d1)) - (strike * Math.exp(-rate * maturity) * norm_cdf(d2))
    end

    def gauss_legendre(n)
      @gl_cache ||= {}
      return @gl_cache[n] if @gl_cache.key?(n)

      nodes = Array.new(n, 0.0)
      weights = Array.new(n, 0.0)
      ((n + 1) / 2).times do |i|
        z = Math.cos(Math::PI * (i + 0.75) / (n + 0.5))
        pp = 0.0
        loop do
          p1 = 1.0
          p2 = 0.0
          (1..n).each do |j|
            p3 = p2
            p2 = p1
            p1 = ((((2.0 * j) - 1.0) * z * p2) - ((j - 1.0) * p3)) / j
          end
          pp = n * ((z * p1) - p2) / ((z * z) - 1.0)
          z1 = z
          z -= p1 / pp
          break if (z - z1).abs < 1.0e-15
        end
        nodes[i] = -z
        nodes[n - 1 - i] = z
        w = 2.0 / ((1.0 - (z * z)) * pp * pp)
        weights[i] = w
        weights[n - 1 - i] = w
      end
      @gl_cache[n] = [nodes, weights]
    end
  end
end
