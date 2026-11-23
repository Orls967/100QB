# frozen_string_literal: true

module CopulaTailDependence
  module Dependence
    module_function

    # Peringkat rata-rata (1..n) dalam O(n log n).
    def average_ranks(xs)
      n = xs.length
      order = (0...n).sort_by { |i| xs[i] }
      ranks = Array.new(n, 0.0)
      i = 0
      while i < n
        j = i + 1
        j += 1 while j < n && xs[order[j]] == xs[order[i]]
        avg = 0.5 * (i + j + 1)
        k = i
        while k < j
          ranks[order[k]] = avg
          k += 1
        end
        i = j
      end
      ranks
    end

    def spearman_rho(xs, ys)
      n = xs.length
      rx = average_ranks(xs)
      ry = average_ranks(ys)
      mx = 0.5 * (n + 1)
      num = 0.0
      dx2 = 0.0
      dy2 = 0.0
      i = 0
      while i < n
        dx = rx[i] - mx
        dy = ry[i] - mx
        num += dx * dy
        dx2 += dx * dx
        dy2 += dy * dy
        i += 1
      end
      den = Math.sqrt(dx2 * dy2)
      den.positive? ? num / den : 0.0
    end

    # Kendall's tau-a/b dalam O(n log n) via Fenwick tree setelah pengurutan berdasarkan xs.
    def kendall_tau(xs, ys)
      n = xs.length
      return 0.0 if n < 2

      order = (0...n).sort_by { |i| [xs[i], ys[i]] }
      ry = average_ranks(ys)
      # Kompresi peringkat ys menjadi integer 1..n untuk BIT.
      y_sorted = order.map { |idx| ry[idx] }
      bit = Array.new(n + 1, 0)
      concordant = 0
      discordant = 0
      i = 0
      while i < n
        # Hitung elemen sebelumnya yang < y_sorted[i] dan > y_sorted[i].
        j = i
        j += 1 while j < n && xs[order[j]] == xs[order[i]]
        # Untuk blok x yang sama, jangan hitung pasangan sesama x-tie.
        k = i
        while k < j
          r = y_sorted[k].ceil
          less = bit_query(bit, r - 1)
          concordant += less
          discordant += (i - bit_query(bit, r))
          k += 1
        end
        k = i
        while k < j
          r = y_sorted[k].ceil
          bit_add(bit, r, 1)
          k += 1
        end
        i = j
      end
      total = concordant + discordant
      total.positive? ? (concordant - discordant).to_f / total : 0.0
    end

    def bit_add(bit, idx, val)
      n = bit.length - 1
      while idx <= n
        bit[idx] += val
        idx += idx & -idx
      end
    end
    private_class_method :bit_add

    def bit_query(bit, idx)
      s = 0
      while idx.positive?
        s += bit[idx]
        idx -= idx & -idx
      end
      s
    end
    private_class_method :bit_query

    def empirical_tail(u1, u2, u_thresh, side: :upper)
      n = u1.length
      count = 0
      i = 0
      if side == :upper
        while i < n
          count += 1 if u1[i] > u_thresh && u2[i] > u_thresh
          i += 1
        end
      else
        low = 1.0 - u_thresh
        while i < n
          count += 1 if u1[i] <= low && u2[i] <= low
          i += 1
        end
      end
      count.to_f / (n * (1.0 - u_thresh))
    end
  end
end
