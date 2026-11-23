# frozen_string_literal: true

module CdoTranche
  module NormInvCoeffs
    A0 = [2.5090809287301226727e3, 3.3430575583588128105e4, 6.7265770927008700853e4,
          4.5921953931549871457e4, 1.3731693765509461125e4, 1.9715909503025514427e3,
          1.3314166789178437745e2, 3.3871328727963666080e0].freeze
    B0 = [5.2264952788528545610e3, 2.8729085735721942674e4, 3.9307895800092710610e4,
          2.1213794301586595867e4, 5.3941960214247511077e3, 6.8718700749205790830e2,
          4.2313330701600911252e1, 1.0].freeze
    A1 = [7.74545014278341407640e-4, 2.27238449892691845833e-2, 2.41780725177450611770e-1,
          1.27045825245236838258e0, 3.64784832476320460504e0, 5.76949722146069140550e0,
          4.63033784615654529590e0, 1.42343711074968357734e0].freeze
    B1 = [1.05075007164441684324e-9, 5.47593808499534494600e-4, 1.51986665636164571966e-2,
          1.48103976427480074590e-1, 6.89767334985100004550e-1, 1.67638483018380384930e0,
          2.05319162663775882187e0, 1.0].freeze
    A2 = [2.01033439929228813265e-7, 2.71155556874348757815e-5, 1.24266094738807843860e-3,
          2.65321895265761230930e-2, 2.96560571828504891230e-1, 1.78482653991729133580e0,
          5.46378491116411436990e0, 6.65790464350110377720e0].freeze
    B2 = [2.04426310338993978564e-15, 1.42151175831644588870e-7, 1.84631831751005468180e-5,
          7.86869131145613259100e-4, 1.48753612908506148525e-2, 1.36929880922735805310e-1,
          5.99832206555887937690e-1, 1.0].freeze
  end

  module Special
    INV_SQRT_2PI = 1.0 / Math.sqrt(2.0 * Math::PI)
    SQRT_2 = Math.sqrt(2.0)

    module_function

    def norm_pdf(x) = INV_SQRT_2PI * Math.exp(-0.5 * x * x)
    def norm_cdf(x) = 0.5 * (1.0 + Math.erf(x / SQRT_2))

    def horner(x, coeffs)
      acc = 0.0
      coeffs.each { |c| acc = (acc * x) + c }
      acc
    end
    private_class_method :horner

    def norm_inv(p)
      raise ArgumentError, "p harus di (0, 1), dapat #{p}" if p <= 0.0 || p >= 1.0

      q = p - 0.5
      if q.abs <= 0.425
        r = 0.180625 - (q * q)
        return q * horner(r, NormInvCoeffs::A0) / horner(r, NormInvCoeffs::B0)
      end
      r = Math.sqrt(-Math.log(q < 0.0 ? p : 1.0 - p))
      val = if r <= 5.0
              horner(r - 1.6,
                     NormInvCoeffs::A1) / horner(r - 1.6,
                                                 NormInvCoeffs::B1)
            else
              horner(r - 5.0,
                     NormInvCoeffs::A2) / horner(
                       r - 5.0, NormInvCoeffs::B2
                     )
            end
      q < 0.0 ? -val : val
    end

    def log_gamma(x) = Math.lgamma(x).first

    def beta_inc(a, b, x)
      return 0.0 if x <= 0.0
      return 1.0 if x >= 1.0

      bt = Math.exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + (a * Math.log(x)) + (b * Math.log(1.0 - x)))
      x < (a + 1.0) / (a + b + 2.0) ? bt * beta_cf(a, b, x) / a : 1.0 - (bt * beta_cf(b, a, 1.0 - x) / b)
    end

    def beta_cf(a, b, x)
      tiny = 1.0e-30
      qab = a + b
      qap = a + 1.0
      qam = a - 1.0
      c = 1.0
      d = [1.0 - (qab * x / qap), tiny].max
      d = 1.0 / d
      h = d
      (1..200).each do |m|
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 / [1.0 + (aa * d), tiny].max_by(&:abs)
        c = [1.0 + (aa / c), tiny].max_by(&:abs)
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 / [1.0 + (aa * d), tiny].max_by(&:abs)
        c = [1.0 + (aa / c), tiny].max_by(&:abs)
        del = d * c
        h *= del
        break if (del - 1.0).abs < 3.0e-14
      end
      h
    end
    private_class_method :beta_cf

    def t_cdf(x, df)
      z = df / (df + (x * x))
      ib = 0.5 * beta_inc(0.5 * df, 0.5, z)
      x >= 0.0 ? 1.0 - ib : ib
    end

    def t_inv(p, df)
      x = norm_inv(p)
      lo = -1.0e4
      hi = 1.0e4
      40.times do
        err = t_cdf(x, df) - p
        return x if err.abs < 1.0e-12

        err.negative? ? (lo = [lo, x].max) : (hi = [hi, x].min)
        x = 0.5 * (lo + hi)
      end
      x
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

    def bivariate_norm_cdf(a, b, rho)
      return norm_cdf(a) * norm_cdf(b) if rho.zero?
      return norm_cdf([a, b].min) if rho >= 1.0 - 1.0e-15
      return [0.0, norm_cdf(a) + norm_cdf(b) - 1.0].max if rho <= -1.0 + 1.0e-15

      nodes, weights = gauss_legendre(32)
      half = 0.5 * Math.asin(rho)
      sum = 0.0
      nodes.each_with_index do |xi, idx|
        u = half * (xi + 1.0)
        c2 = Math.cos(u)**2
        next if c2 < 1.0e-15

        expo = -(((a * a) - (2.0 * Math.sin(u) * a * b) + (b * b)) / (2.0 * c2))
        sum += weights[idx] * Math.exp(expo) if expo > -700.0
      end
      (norm_cdf(a) * norm_cdf(b)) + (sum * half / (2.0 * Math::PI))
    end
  end
end
