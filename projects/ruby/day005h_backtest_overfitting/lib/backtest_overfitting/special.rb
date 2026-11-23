# frozen_string_literal: true

module BacktestOverfitting
  module Special
    SQRT_2 = Math.sqrt(2.0)
    SQRT_2PI = Math.sqrt(2.0 * Math::PI)
    EULER_GAMMA = 0.5772156649015328606

    module_function

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

    def binom(n, k)
      return 0 if k.negative? || k > n

      m = [k, n - k].min
      return 1 if m.zero?

      num = 1
      den = 1
      m.times do |i|
        num *= (n - i)
        den *= (i + 1)
      end
      num / den
    end

    def beta_inc(a, b, x)
      return 0.0 if x <= 0.0
      return 1.0 if x >= 1.0

      bt = Math.exp(Math.lgamma(a + b)[0] - Math.lgamma(a)[0] - Math.lgamma(b)[0] +
                    (a * Math.log(x)) + (b * Math.log(1.0 - x)))
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
      d = [1.0 - (qab * x / qap), 1.0e-30].max_by(&:abs)
      d = 1.0e-30 if d.abs < 1.0e-30
      d = 1.0 / d
      h = d
      (1..120).each do |m|
        c, d, h = beta_cf_step(m, a, b, x, qab, qap, qam, c, d, h)
      end
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

    def student_t_pvalue(t_stat, df)
      x = df.to_f / (df + (t_stat * t_stat))
      beta_inc(0.5 * df, 0.5, x)
    end

    def gauss_legendre(n, lo = -1.0, hi = 1.0)
      nodes, weights = gauss_legendre_unit(n)
      half = 0.5 * (hi - lo)
      mid = 0.5 * (hi + lo)
      [nodes.map { |x| mid + (half * x) }, weights.map { |w| half * w }]
    end

    def gauss_legendre_unit(n)
      @gl_cache ||= {}
      return @gl_cache[n] if @gl_cache.key?(n)

      nodes = Array.new(n, 0.0)
      weights = Array.new(n, 0.0)
      ((n + 1) / 2).times { |i| gl_fill_pair(i, n, nodes, weights) }
      @gl_cache[n] = [nodes, weights]
    end
    private_class_method :gauss_legendre_unit

    def gl_fill_pair(idx, n, nodes, weights)
      z = Math.cos(Math::PI * (idx + 0.75) / (n + 0.5))
      pp = 0.0
      loop do
        p1, p2 = gl_legendre_eval(n, z)
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

    def gl_legendre_eval(n, z)
      p1 = 1.0
      p2 = 0.0
      (1..n).each do |j|
        p3 = p2
        p2 = p1
        p1 = ((((2.0 * j) - 1.0) * z * p2) - ((j - 1.0) * p3)) / j
      end
      [p1, p2]
    end
    private_class_method :gl_legendre_eval

    def pearson_corr(xs, ys)
      n = xs.length
      mx = xs.sum(0.0) / n
      my = ys.sum(0.0) / n
      num = 0.0
      sx = 0.0
      sy = 0.0
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

    def ols_fit(xs, ys)
      n = xs.length
      mx = xs.sum(0.0) / n
      my = ys.sum(0.0) / n
      num = 0.0
      den = 0.0
      n.times do |i|
        dx = xs[i] - mx
        num += dx * (ys[i] - my)
        den += dx * dx
      end
      slope = den.positive? ? num / den : 0.0
      { slope: slope, intercept: my - (slope * mx) }
    end

    def ranks(values)
      n = values.length
      order = (0...n).sort_by { |i| values[i] }
      out = Array.new(n, 0.0)
      i = 0
      while i < n
        j = i + 1
        j += 1 while j < n && values[order[j]] == values[order[i]]
        avg_rank = 0.5 * (i + 1 + j)
        (i...j).each { |k| out[order[k]] = avg_rank }
        i = j
      end
      out
    end
  end
end
