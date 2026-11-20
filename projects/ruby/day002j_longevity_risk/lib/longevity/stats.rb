# frozen_string_literal: true

module Longevity
  module Stats
    module_function

    def mean(arr)
      return 0.0 if arr.empty?

      arr.sum(0.0) / arr.size.to_f
    end

    def variance(arr, ddof: 1)
      n = arr.size
      return 0.0 if n <= ddof

      m = mean(arr)
      arr.sum(0.0) { |x| (x - m)**2 } / (n - ddof).to_f
    end

    def std_dev(arr, ddof: 1)
      Math.sqrt([0.0, variance(arr, ddof: ddof)].max)
    end

    def correlation(xs_vals, ys_vals)
      mx = mean(xs_vals)
      my = mean(ys_vals)
      cov = 0.0
      vx = 0.0
      vy = 0.0
      xs_vals.each_index do |i|
        dx = xs_vals[i] - mx
        dy = ys_vals[i] - my
        cov += dx * dy
        vx += dx * dx
        vy += dy * dy
      end
      denom = Math.sqrt(vx * vy)
      denom.zero? ? 0.0 : cov / denom
    end

    def quantile_sorted(sorted, prob)
      return 0.0 if sorted.empty?
      return sorted.first if prob <= 0.0
      return sorted.last if prob >= 1.0

      pos = prob * (sorted.size - 1)
      lo = pos.floor
      hi = pos.ceil
      frac = pos - lo
      (sorted[lo] * (1.0 - frac)) + (sorted[hi] * frac)
    end

    def var_and_es(values, var_prob: 0.995, es_prob: 0.99)
      sorted = values.sort
      var_val = quantile_sorted(sorted, var_prob)
      tail_start = (es_prob * sorted.size).floor.clamp(0, sorted.size - 1)
      tail = sorted[tail_start..]
      es_val = mean(tail)
      { var: var_val, es: es_val }
    end

    def batch_means_se(values, num_batches: 20)
      n = values.size
      b_count = [num_batches, n / 10].min.clamp(2, n)
      b_size = n / b_count
      b_means = []
      b_vars = []
      b_count.times do |b|
        slice = values[(b * b_size)...((b + 1) * b_size)]
        b_means << mean(slice)
        b_vars << variance(slice, ddof: 1)
      end
      {
        mean_se: std_dev(b_means, ddof: 1) / Math.sqrt(b_count.to_f),
        var_se: std_dev(b_vars, ddof: 1) / Math.sqrt(b_count.to_f)
      }
    end
  end
end
