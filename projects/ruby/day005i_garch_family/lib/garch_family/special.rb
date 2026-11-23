# frozen_string_literal: true

module GarchFamily
  module Special
    SQRT_2 = Math.sqrt(2.0)
    SQRT_2PI = Math.sqrt(2.0 * Math::PI)
    SQRT_2_OVER_PI = Math.sqrt(2.0 / Math::PI)

    module_function

    def lgamma(x) = Math.lgamma(x)[0]

    def digamma(x)
      res = 0.0
      xx = x
      while xx < 7.0
        res -= 1.0 / xx
        xx += 1.0
      end
      inv = 1.0 / xx
      inv2 = inv * inv
      res + Math.log(xx) - (0.5 * inv) - (inv2 * ((1.0 / 12.0) - (inv2 * ((1.0 / 120.0) - (inv2 / 252.0)))))
    end

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
      (1..140).each { |m| c, d, h = beta_cf_step(m, a, b, x, qab, qap, qam, c, d, h) }
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
      160.times do
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
      (1..160).each do |i|
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

    def chi2_mixture_01_sf(x) = x <= 0.0 ? 1.0 : 0.5 * chi2_sf(x, 1.0)

    def student_t_pdf(t, df)
      c = Math.exp(lgamma(0.5 * (df + 1.0)) - lgamma(0.5 * df)) / Math.sqrt(Math::PI * df)
      c * ((1.0 + ((t * t) / df))**(-0.5 * (df + 1.0)))
    end

    def student_t_cdf(t, df)
      return 0.5 if t.zero?

      x = df.to_f / (df + (t * t))
      ib = beta_inc(0.5 * df, 0.5, x)
      t.positive? ? 1.0 - (0.5 * ib) : 0.5 * ib
    end

    def student_t_inv(prob, df)
      raise ArgumentError, 'prob harus di (0, 1)' if prob <= 0.0 || prob >= 1.0
      return 0.0 if (prob - 0.5).abs < 1.0e-16
      return -student_t_inv(1.0 - prob, df) if prob > 0.5

      lo = -200.0
      hi = 0.0
      80.times do
        mid = 0.5 * (lo + hi)
        if student_t_cdf(mid, df) < prob
          lo = mid
        else
          hi = mid
        end
      end
      0.5 * (lo + hi)
    end

    def std_t_inv(prob, df) = student_t_inv(prob, df) * Math.sqrt((df - 2.0) / df)

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

    def gauss_legendre(n, lo = -1.0, hi = 1.0)
      nodes = Array.new(n, 0.0)
      weights = Array.new(n, 0.0)
      ((n + 1) / 2).times { |i| gl_fill_pair(i, n, nodes, weights) }
      half = 0.5 * (hi - lo)
      mid = 0.5 * (hi + lo)
      [nodes.map { |x| mid + (half * x) }, weights.map { |w| half * w }]
    end

    def gl_fill_pair(idx, n, nodes, weights)
      z = Math.cos(Math::PI * (idx + 0.75) / (n + 0.5))
      pp = 0.0
      loop do
        p1, p2 = gl_eval(n, z)
        pp = n * ((z * p1) - p2) / ((z * z) - 1.0)
        dz = p1 / pp
        z -= dz
        break if dz.abs < 1.0e-15
      end
      nodes[idx] = -z
      nodes[n - 1 - idx] = z
      w = 2.0 / ((1.0 - (z * z)) * pp * pp)
      weights[idx] = w
      weights[n - 1 - idx] = w
    end
    private_class_method :gl_fill_pair

    def gl_eval(n, z)
      p1 = 1.0
      p2 = 0.0
      (1..n).each do |j|
        p3 = p2
        p2 = p1
        p1 = ((((2.0 * j) - 1.0) * z * p2) - ((j - 1.0) * p3)) / j
      end
      [p1, p2]
    end
    private_class_method :gl_eval
  end
end
