# frozen_string_literal: true

module FactorPortfolioLab
  module Distributions
    SQRT_2 = Math.sqrt(2.0)
    SQRT_2PI = Math.sqrt(2.0 * Math::PI)

    module_function

    def lgamma(x) = Math.lgamma(x)[0]

    def norm_pdf(x) = Math.exp(-0.5 * x * x) / SQRT_2PI

    def norm_cdf(x)
      return 0.5 * Math.erfc(-x / SQRT_2) if x.negative?

      0.5 * (1.0 + Math.erf(x / SQRT_2))
    end

    def norm_inv(prob)
      raise ArgumentError, 'prob harus di (0, 1)' if prob <= 0.0 || prob >= 1.0
      return 0.0 if (prob - 0.5).abs < 1.0e-16
      return -norm_inv(1.0 - prob) if prob > 0.5

      t = Math.sqrt(-2.0 * Math.log(prob))
      x = -(t - ((2.515517 + (0.802853 * t) + (0.010328 * t * t)) /
                 (1.0 + (1.432788 * t) + (0.189269 * t * t) + (0.001308 * t * t * t))))
      4.times do
        err = norm_cdf(x) - prob
        u = err / [norm_pdf(x), 1.0e-300].max
        denom = 1.0 + (0.5 * x * u)
        x -= denom.abs > 0.1 ? (u / denom) : u
      end
      x
    end

    def beta_inc(a, b, x)
      return 0.0 if x <= 0.0
      return 1.0 if x >= 1.0

      bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + (a * Math.log(x)) + (b * Math.log(1.0 - x)))
      if x < (a + 1.0) / (a + b + 2.0)
        bt * beta_cf(a, b, x) / a
      else
        1.0 - (bt * beta_cf(b, a, 1.0 - x) / b)
      end
    end

    def beta_cf(a, b, x)
      qab = a + b
      qap = a + 1.0
      qam = a - 1.0
      c = 1.0
      d = 1.0 / tiny_guard(1.0 - (qab * x / qap))
      h = d
      (1..160).each { |m| c, d, h = beta_cf_step(m, a, b, x, qab, qap, qam, c, d, h) }
      h
    end
    private_class_method :beta_cf

    def beta_cf_step(m, a, b, x, qab, qap, qam, c, d, h)
      m2 = 2 * m
      aa = m * (b - m) * x / ((qam + m2) * (a + m2))
      d = 1.0 / tiny_guard(1.0 + (aa * d))
      c = tiny_guard(1.0 + (aa / c))
      h *= d * c
      aa2 = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
      d = 1.0 / tiny_guard(1.0 + (aa2 * d))
      c = tiny_guard(1.0 + (aa2 / c))
      [c, d, h * d * c]
    end
    private_class_method :beta_cf_step

    def tiny_guard(val) = val.abs < 1.0e-30 ? 1.0e-30 : val
    private_class_method :tiny_guard

    def gamma_inc_q(a, x)
      return 1.0 if x <= 0.0
      return 1.0 - gamma_series(a, x) if x < a + 1.0

      gamma_cf(a, x)
    end

    def gamma_series(a, x)
      gln = lgamma(a)
      ap = a
      sum = 1.0 / a
      del = sum
      180.times do
        ap += 1.0
        del *= x / ap
        sum += del
        break if del.abs < sum.abs * 1.0e-14
      end
      sum * Math.exp(-x + (a * Math.log(x)) - gln)
    end
    private_class_method :gamma_series

    def gamma_cf(a, x)
      gln = lgamma(a)
      b = x + 1.0 - a
      c = 1.0e30
      d = 1.0 / tiny_guard(b)
      h = d
      (1..180).each do |i|
        an = -i * (i - a)
        b += 2.0
        d = 1.0 / tiny_guard((an * d) + b)
        c = tiny_guard(b + (an / c))
        del = d * c
        h *= del
        break if (del - 1.0).abs < 1.0e-14
      end
      Math.exp(-x + (a * Math.log(x)) - gln) * h
    end
    private_class_method :gamma_cf

    def chi2_sf(x, df) = x <= 0.0 ? 1.0 : gamma_inc_q(0.5 * df, 0.5 * x)

    def chi2_cdf(x, df) = 1.0 - chi2_sf(x, df)

    def f_cdf(x, d1, d2)
      return 0.0 if x <= 0.0

      z = (d1 * x) / ((d1 * x) + d2)
      beta_inc(0.5 * d1, 0.5 * d2, z)
    end

    def f_sf(x, d1, d2) = 1.0 - f_cdf(x, d1, d2)

    def f_inv(prob, d1, d2)
      raise ArgumentError, 'prob harus di (0, 1)' if prob <= 0.0 || prob >= 1.0

      hi = 2.0
      hi *= 2.0 while f_cdf(hi, d1, d2) < prob && hi < 1.0e8
      RootFinder.brent(0.0, hi, tol: 1.0e-10) { |x| f_cdf(x, d1, d2) - prob }
    end

    def student_t_cdf(t, df)
      return 0.5 if t.zero?

      x = df.to_f / (df + (t * t))
      ib = beta_inc(0.5 * df, 0.5, x)
      t.positive? ? 1.0 - (0.5 * ib) : 0.5 * ib
    end

    def student_t_inv(prob, df)
      return 0.0 if (prob - 0.5).abs < 1.0e-16
      return -student_t_inv(1.0 - prob, df) if prob > 0.5

      RootFinder.brent(-200.0, 0.0, tol: 1.0e-10) { |t| student_t_cdf(t, df) - prob }
    end

    def pearson_corr(xs, ys)
      n = xs.length
      mx = xs.sum(0.0) / n
      my = ys.sum(0.0) / n
      num = sx = sy = 0.0
      n.times do |i|
        dx = xs[i] - mx
        dy = ys[i] - my
        num += dx * dy
        sx += dx * dx
        sy += dy * dy
      end
      denom = Math.sqrt(sx * sy)
      denom.positive? ? num / denom : 0.0
    end
  end
end
