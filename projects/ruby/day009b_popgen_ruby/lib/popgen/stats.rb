# frozen_string_literal: true

module PopGen
  # Summary statistics, block averaging, chi-square p-values, KS, OLS, FNV-1a, and compression.
  module Stats
    FNV_OFFSET = 14_695_981_039_346_656_037
    FNV_PRIME = 1_099_511_628_211
    MASK64 = 0xFFFFFFFFFFFFFFFF

    module_function

    def mean(xs)
      return 0.0 if xs.empty?

      xs.sum(0.0) / xs.length
    end

    def variance(xs)
      n = xs.length
      return 0.0 if n <= 1

      m = mean(xs)
      ss = xs.sum(0.0) { |x| (x - m) * (x - m) }
      ss / (n - 1)
    end

    def sd(xs)
      Math.sqrt(variance(xs))
    end

    def se_iid(xs)
      n = xs.length
      return 0.0 if n <= 1

      Math.sqrt(variance(xs) / n)
    end

    def se_proportion(p, n)
      return 0.0 if n <= 0

      Math.sqrt([p * (1.0 - p), 0.0].max / n)
    end

    def block_average(xs, n_blocks = 10)
      n = xs.length
      nb = [n_blocks, n].min
      nb = 2 if nb < 2 && n >= 2
      bsz = n / nb
      bsz = 1 if bsz < 1
      b_means = Array.new(nb) do |b|
        slice = xs[(b * bsz)...((b + 1) * bsz)]
        mean(slice)
      end
      m = mean(b_means)
      se = Math.sqrt(variance(b_means) / nb)
      { mean: m, se: se, sd: sd(xs), n_blocks: nb, block_size: bsz, block_means: b_means }
    end

    def median_mad(xs)
      return [0.0, 0.0] if xs.empty?

      s = xs.sort
      n = s.length
      med = n.odd? ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2])
      devs = s.map { |x| (x - med).abs }.sort
      mad = n.odd? ? devs[n / 2] : 0.5 * (devs[n / 2 - 1] + devs[n / 2])
      [med, mad]
    end

    def ols_fit(xs, ys)
      n = xs.length
      return { slope: 0.0, intercept: 0.0, slope_se: 0.0, r2: 0.0 } if n < 2

      mx = mean(xs)
      my = mean(ys)
      sxx = 0.0
      sxy = 0.0
      syy = 0.0
      n.times do |i|
        dx = xs[i] - mx
        dy = ys[i] - my
        sxx += dx * dx
        sxy += dx * dy
        syy += dy * dy
      end
      return { slope: 0.0, intercept: my, slope_se: 0.0, r2: 0.0 } if sxx <= 1e-30

      slope = sxy / sxx
      intercept = my - slope * mx
      sse = 0.0
      n.times do |i|
        r = ys[i] - (intercept + slope * xs[i])
        sse += r * r
      end
      r2 = syy > 1e-30 ? (1.0 - sse / syy) : 1.0
      slope_se = n > 2 ? Math.sqrt((sse / (n - 2)) / sxx) : 0.0
      { slope: slope, intercept: intercept, slope_se: slope_se, r2: r2 }
    end

    def pearson_corr(xs, ys)
      n = xs.length
      return 0.0 if n < 2

      mx = mean(xs)
      my = mean(ys)
      sxx = 0.0
      syy = 0.0
      sxy = 0.0
      n.times do |i|
        dx = xs[i] - mx
        dy = ys[i] - my
        sxx += dx * dx
        syy += dy * dy
        sxy += dx * dy
      end
      den = Math.sqrt(sxx * syy)
      den > 1e-30 ? (sxy / den) : 0.0
    end

    # Upper regularized incomplete gamma Q(a, x) = Gamma(a, x) / Gamma(a)
    def upper_gamma_q(a, x)
      return 1.0 if x <= 0.0

      gln = Math.lgamma(a).first
      if x < a + 1.0
        # Series representation for P(a, x), then Q = 1 - P
        ap = a
        sum = 1.0 / a
        del = sum
        200.times do
          ap += 1.0
          del *= x / ap
          sum += del
          break if del.abs < sum.abs * 1e-14
        end
        p_val = sum * Math.exp(-x + a * Math.log(x) - gln)
        [[1.0 - p_val, 0.0].max, 1.0].min
      else
        # Lentz continued fraction for Q(a, x)
        b = x + 1.0 - a
        c = 1.0 / 1e-30
        d = 1.0 / b
        h = d
        (1..200).each do |i|
          an = -i * (i - a)
          b += 2.0
          d = an * d + b
          d = 1e-30 if d.abs < 1e-30
          c = b + an / c
          c = 1e-30 if c.abs < 1e-30
          d = 1.0 / d
          del = d * c
          h *= del
          break if (del - 1.0).abs < 1e-14
        end
        q_val = Math.exp(-x + a * Math.log(x) - gln) * h
        [[q_val, 0.0].max, 1.0].min
      end
    end

    def chi2_p_value(chi2, df)
      return 1.0 if df <= 0 || chi2 <= 0.0

      upper_gamma_q(0.5 * df, 0.5 * chi2)
    end

    # Chi-square goodness-of-fit pooling tail bins whose expected count < min_exp
    def chi2_gof_pooled(observed, expected, min_exp: 5.0)
      obs_p = []
      exp_p = []
      acc_o = 0.0
      acc_e = 0.0
      observed.each_with_index do |o, i|
        acc_o += o.to_f
        acc_e += expected[i].to_f
        if acc_e >= min_exp
          obs_p << acc_o
          exp_p << acc_e
          acc_o = 0.0
          acc_e = 0.0
        end
      end
      if acc_e.positive?
        if exp_p.empty?
          obs_p << acc_o
          exp_p << acc_e
        else
          obs_p[-1] += acc_o
          exp_p[-1] += acc_e
        end
      end
      chi2 = 0.0
      obs_p.each_with_index do |o, i|
        e = exp_p[i]
        chi2 += (o - e) * (o - e) / e if e > 0.0
      end
      df = [obs_p.length - 1, 1].max
      { chi2: chi2, df: df, p_value: chi2_p_value(chi2, df), bins_used: obs_p.length }
    end

    def ks_against_cdf(samples, cdf_proc)
      s = samples.sort
      n = s.length
      return { d: 0.0, p_value: 1.0 } if n.zero?

      d_max = 0.0
      n.times do |i|
        f = cdf_proc.call(s[i])
        d1 = (f - i.to_f / n).abs
        d2 = ((i + 1).to_f / n - f).abs
        d_max = d1 if d1 > d_max
        d_max = d2 if d2 > d_max
      end
      en = Math.sqrt(n)
      lam = (en + 0.12 + 0.11 / en) * d_max
      p = if lam < 1e-8
            1.0
          else
            sum = 0.0
            (1..100).each do |k|
              term = 2.0 * Math.exp(-2.0 * k * k * lam * lam)
              sum += (k.odd? ? term : -term)
              break if term < 1e-12
            end
            [[sum, 0.0].max, 1.0].min
          end
      { d: d_max, p_value: p }
    end

    def fnv1a_bytes(str, init = FNV_OFFSET)
      h = init
      str.each_byte do |b|
        h = ((h ^ b) * FNV_PRIME) & MASK64
      end
      h
    end

    def fnv1a_hex(str)
      format('%016x', fnv1a_bytes(str))
    end

    # Quantized delta + RLE + ZigZag varint lossless roundtrip on trajectory matrix
    def compress_trajectory_roundtrip(frames, quant_step = 1e-4)
      n_frames = frames.length
      n_vals = frames[0].length
      raw_bytes = n_frames * n_vals * 8

      q = frames.map { |row| row.map { |v| (v / quant_step).round } }
      deltas = Array.new(n_frames) { Array.new(n_vals, 0) }
      n_vals.times { |j| deltas[0][j] = q[0][j] }
      (1...n_frames).each do |f|
        n_vals.times { |j| deltas[f][j] = q[f][j] - q[f - 1][j] }
      end

      flat = deltas.flatten
      encoded = +''
      encoded.force_encoding(Encoding::ASCII_8BIT)
      i = 0
      tot = flat.length
      while i < tot
        if flat[i].zero?
          run = 0
          while i + run < tot && flat[i + run].zero? && run < 65_535
            run += 1
          end
          encoded << encode_zigzag_varint(0)
          encoded << encode_zigzag_varint(run)
          i += run
        else
          encoded << encode_zigzag_varint(flat[i])
          i += 1
        end
      end

      # Decode and reconstruct
      dec = []
      pos = 0
      bytes = encoded.bytes
      while pos < bytes.length
        val, pos = decode_zigzag_varint(bytes, pos)
        if val.zero?
          run, pos = decode_zigzag_varint(bytes, pos)
          run.times { dec << 0 }
        else
          dec << val
        end
      end

      rec = Array.new(n_frames) { Array.new(n_vals, 0) }
      n_vals.times { |j| rec[0][j] = dec[j] }
      (1...n_frames).each do |f|
        n_vals.times { |j| rec[f][j] = rec[f - 1][j] + dec[f * n_vals + j] }
      end

      max_err_q = 0
      n_frames.times do |f|
        n_vals.times do |j|
          err = (rec[f][j] - q[f][j]).abs
          max_err_q = err if err > max_err_q
        end
      end

      comp_bytes = encoded.bytesize
      {
        raw_bytes: raw_bytes,
        comp_bytes: comp_bytes,
        ratio: comp_bytes.positive? ? (raw_bytes.to_f / comp_bytes) : 0.0,
        exact_roundtrip: (dec.length == tot && max_err_q.zero?),
        max_err_q: max_err_q
      }
    end

    def encode_zigzag_varint(v)
      zz = v >= 0 ? (v << 1) : ((-v << 1) - 1)
      out = +''
      out.force_encoding(Encoding::ASCII_8BIT)
      while zz >= 0x80
        out << ((zz & 0x7F) | 0x80)
        zz >>= 7
      end
      out << zz
      out
    end

    def decode_zigzag_varint(bytes, pos)
      zz = 0
      shift = 0
      loop do
        b = bytes[pos]
        pos += 1
        zz |= ((b & 0x7F) << shift)
        break if (b & 0x80).zero?

        shift += 7
      end
      val = zz.even? ? (zz >> 1) : -((zz + 1) >> 1)
      [val, pos]
    end
  end
end
