# frozen_string_literal: true

module AsianOptions
  module Stats
    CI_Z95 = 1.96
    INV_SQRT2 = 1.0 / Math.sqrt(2.0)

    UnivariateChunk = Data.define(:sample_count, :mean, :m2)
    BivariateChunk = Data.define(:sample_count, :mean_x, :mean_y, :m2_x, :m2_y, :c_xy)

    module_function

    def norm_cdf(value)
      0.5 * Math.erfc(-value * INV_SQRT2)
    end

    def merge_univariate(chunks)
      acc = UnivariateChunk.new(sample_count: 0, mean: 0.0, m2: 0.0)
      chunks.each do |chunk|
        next if chunk.sample_count.zero?

        acc = acc.sample_count.zero? ? chunk : combine_univariate(acc, chunk)
      end
      acc
    end

    def combine_univariate(left, right)
      total = left.sample_count + right.sample_count
      delta = right.mean - left.mean
      scale = (left.sample_count * right.sample_count).to_f / total
      UnivariateChunk.new(
        sample_count: total,
        mean: left.mean + (delta * (right.sample_count.to_f / total)),
        m2: left.m2 + right.m2 + (delta * delta * scale)
      )
    end
    private_class_method :combine_univariate

    def summarize_univariate(merged)
      raise ArgumentError, 'jumlah sampel harus >= 2' if merged.sample_count < 2

      n = merged.sample_count.to_f
      variance = [merged.m2 / (n - 1.0), 0.0].max
      se = Math.sqrt(variance / n)
      EstimatorResult.new(
        price: merged.mean,
        se: se,
        ci_low: merged.mean - (CI_Z95 * se),
        ci_high: merged.mean + (CI_Z95 * se)
      )
    end

    def merge_bivariate(chunks)
      acc = BivariateChunk.new(sample_count: 0, mean_x: 0.0, mean_y: 0.0, m2_x: 0.0, m2_y: 0.0, c_xy: 0.0)
      chunks.each do |chunk|
        next if chunk.sample_count.zero?

        acc = acc.sample_count.zero? ? chunk : combine_bivariate(acc, chunk)
      end
      acc
    end

    def combine_bivariate(left, right)
      total = left.sample_count + right.sample_count
      scale = (left.sample_count * right.sample_count).to_f / total
      dx = right.mean_x - left.mean_x
      dy = right.mean_y - left.mean_y
      ratio = right.sample_count.to_f / total
      build_combined_bivariate(left, right, total, ratio, scale, dx, dy)
    end
    private_class_method :combine_bivariate

    def build_combined_bivariate(left, right, total, ratio, scale, delta_x, delta_y)
      BivariateChunk.new(
        sample_count: total, mean_x: left.mean_x + (delta_x * ratio), mean_y: left.mean_y + (delta_y * ratio),
        m2_x: left.m2_x + right.m2_x + (delta_x * delta_x * scale),
        m2_y: left.m2_y + right.m2_y + (delta_y * delta_y * scale),
        c_xy: left.c_xy + right.c_xy + (delta_x * delta_y * scale)
      )
    end
    private_class_method :build_combined_bivariate

    def fit_loglog_slope(counts, ses)
      if counts.length != ses.length || counts.length < 2
        raise ArgumentError, 'panjang counts dan ses harus sama dan >= 2'
      end

      xs = counts.map { |n| Math.log10(n.to_f) }
      ys = ses.map { |se| Math.log10(se.to_f) }
      mx = xs.sum / xs.length
      my = ys.sum / ys.length

      num = 0.0
      den = 0.0
      xs.each_index do |i|
        dx = xs[i] - mx
        num += dx * (ys[i] - my)
        den += dx * dx
      end

      raise ArgumentError, 'varians x nol pada regresi log-log' if den <= 0.0

      num / den
    end

    def variance_reduction_factor(se_naive, se_method)
      raise ArgumentError, 'se_method harus > 0' unless se_method.positive?

      (se_naive / se_method)**2
    end

    def z_score(diff, se_a, se_b = 0.0)
      denom = Math.hypot(se_a, se_b)
      raise ArgumentError, 'SE gabungan harus > 0' unless denom.positive?

      diff.abs / denom
    end
  end
end
