# frozen_string_literal: true

module CopulaTailDependence
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

    def norm_pdf(x)
      INV_SQRT_2PI * Math.exp(-0.5 * x * x)
    end

    def norm_cdf(x)
      0.5 * (1.0 + Math.erf(x / SQRT_2))
    end

    def horner(x, coeffs)
      acc = 0.0
      coeffs.each { |c| acc = (acc * x) + c }
      acc
    end
    private_class_method :horner

    # Wichura AS241 rational approximation untuk kuantil normal standar (< 1e-15).
    def norm_inv(p)
      raise ArgumentError, "p harus di (0, 1), dapat #{p}" if p <= 0.0 || p >= 1.0

      q = p - 0.5
      if q.abs <= 0.425
        r = 0.180625 - (q * q)
        return q * horner(r, NormInvCoeffs::A0) / horner(r, NormInvCoeffs::B0)
      end

      r = Math.sqrt(-Math.log(q < 0.0 ? p : 1.0 - p))
      val = if r <= 5.0
              s = r - 1.6
              horner(s, NormInvCoeffs::A1) / horner(s, NormInvCoeffs::B1)
            else
              s = r - 5.0
              horner(s, NormInvCoeffs::A2) / horner(s, NormInvCoeffs::B2)
            end
      q < 0.0 ? -val : val
    end

    def log_gamma(x)
      Math.lgamma(x).first
    end

    # Pecahan lanjut Lentz untuk regularized incomplete beta I_x(a, b).
    def beta_inc(a, b, x)
      return 0.0 if x <= 0.0
      return 1.0 if x >= 1.0

      bt = Math.exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + (a * Math.log(x)) + (b * Math.log(1.0 - x)))
      if x < (a + 1.0) / (a + b + 2.0)
        bt * beta_cf(a, b, x) / a
      else
        1.0 - (bt * beta_cf(b, a, 1.0 - x) / b)
      end
    end

    def beta_cf(a, b, x)
      tiny = 1.0e-30
      eps = 3.0e-14
      qab = a + b
      qap = a + 1.0
      qam = a - 1.0
      c = 1.0
      d = 1.0 - (qab * x / qap)
      d = tiny if d.abs < tiny
      d = 1.0 / d
      h = d

      (1..200).each do |m|
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + (aa * d)
        d = tiny if d.abs < tiny
        c = 1.0 + (aa / c)
        c = tiny if c.abs < tiny
        d = 1.0 / d
        h *= d * c

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + (aa * d)
        d = tiny if d.abs < tiny
        c = 1.0 + (aa / c)
        c = tiny if c.abs < tiny
        d = 1.0 / d
        del = d * c
        h *= del
        break if (del - 1.0).abs < eps
      end
      h
    end
    private_class_method :beta_cf

    def t_cdf(x, df)
      return norm_cdf(x) if df.infinite?

      z = df / (df + (x * x))
      ib = 0.5 * beta_inc(0.5 * df, 0.5, z)
      x >= 0.0 ? 1.0 - ib : ib
    end

    def t_pdf(x, df)
      return norm_pdf(x) if df.infinite?

      c = Math.exp(log_gamma(0.5 * (df + 1.0)) - log_gamma(0.5 * df)) / Math.sqrt(df * Math::PI)
      c * ((1.0 + ((x * x) / df))**(-0.5 * (df + 1.0)))
    end

    def t_inv(p, df)
      raise ArgumentError, "p harus di (0, 1), dapat #{p}" if p <= 0.0 || p >= 1.0
      return norm_inv(p) if df.infinite?

      x = norm_inv(p)
      lo = -1.0e4
      hi = 1.0e4
      40.times do
        err = t_cdf(x, df) - p
        return x if err.abs < 1.0e-13

        if err.negative?
          lo = x if x > lo
        elsif x < hi
          hi = x
        end
        pdf = t_pdf(x, df)
        nxt = pdf > 1.0e-18 ? x - (err / pdf) : 0.5 * (lo + hi)
        x = nxt > lo && nxt < hi ? nxt : 0.5 * (lo + hi)
      end
      x
    end

    # Node dan bobot Gauss-Legendre pada [-1, 1].
    def gauss_legendre(n)
      @gl_cache ||= {}
      return @gl_cache[n] if @gl_cache.key?(n)

      m = (n + 1) / 2
      nodes = Array.new(n, 0.0)
      weights = Array.new(n, 0.0)
      m.times do |i|
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
          z = z1 - (p1 / pp)
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

    # Probabilitas ortan bivariat Gaussian P(Z1 <= a, Z2 <= b) via integral Plackett
    # dengan substitusi r = sin(u) agar singularitas 1/sqrt(1-r^2) hilang di |rho| -> 1.
    def bivariate_norm_cdf(a, b, rho)
      return norm_cdf(a) * norm_cdf(b) if rho.zero?
      return norm_cdf([a, b].min) if rho >= 1.0 - 1.0e-15
      return [0.0, norm_cdf(a) + norm_cdf(b) - 1.0].max if rho <= -1.0 + 1.0e-15

      nodes, weights = gauss_legendre(32)
      u_max = Math.asin(rho)
      half = 0.5 * u_max
      sum = 0.0
      two_pi_inv = 1.0 / (2.0 * Math::PI)
      a2 = a * a
      b2 = b * b
      ab2 = 2.0 * a * b

      nodes.each_with_index do |xi, idx|
        u = half * (xi + 1.0)
        r = Math.sin(u)
        c2 = Math.cos(u)**2
        next if c2 < 1.0e-15

        expo = -((a2 - (r * ab2) + b2) / (2.0 * c2))
        sum += weights[idx] * Math.exp(expo) if expo > -700.0
      end

      (norm_cdf(a) * norm_cdf(b)) + (two_pi_inv * half * sum)
    end
  end
end
