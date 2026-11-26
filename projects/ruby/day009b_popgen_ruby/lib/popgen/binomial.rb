# frozen_string_literal: true

module PopGen
  # Exact Binomial sampler with two branches:
  # - Inversion (inv) for n * min(p, 1-p) < 10
  # - BTPE (Kachitvichyanukul & Schmeiser 1988) exact rejection for n * min(p, 1-p) >= 10
  module Binomial
    module_function

    def log_choose(n, k)
      return 0.0 if k <= 0 || k >= n

      Math.lgamma(n + 1).first - Math.lgamma(k + 1).first - Math.lgamma(n - k + 1).first
    end

    def pmf(k, n, p)
      return 0.0 if k < 0 || k > n
      return (k.zero? ? 1.0 : 0.0) if p <= 0.0
      return (k == n ? 1.0 : 0.0) if p >= 1.0

      lp = log_choose(n, k) + k * Math.log(p) + (n - k) * Math.log(1.0 - p)
      Math.exp(lp)
    end

    def pmf_array(n, p)
      Array.new(n + 1) { |k| pmf(k, n, p) }
    end

    def sample_inv(n, r, rng)
      q = 1.0 - r
      s = r / q
      a = (n + 1) * s
      p0 = Math.exp(n * Math.log(q))
      u = rng.rand
      k = 0
      pk = p0
      while u > pk && k < n
        u -= pk
        k += 1
        pk = pk * (a / k - s)
        break if pk <= 0.0
      end
      k
    end

    def sample_btpe(n, r, rng)
      q = 1.0 - r
      fm = n * r + r
      m = fm.floor
      p1 = (2.195 * Math.sqrt(n * r * q) - 4.6 * q).floor + 0.5
      xm = m + 0.5
      xl = xm - p1
      xr = xm + p1
      c = 0.134 + 20.5 / (15.3 + m)
      al = (fm - xl) / (fm - xl * r)
      lam_l = al * (1.0 + 0.5 * al)
      ar = (xr - fm) / (xr * q)
      lam_r = ar * (1.0 + 0.5 * ar)
      p2 = p1 * (1.0 + 2.0 * c)
      p3 = p2 + c / lam_l
      p4 = p3 + c / lam_r

      s = r / q
      a = s * (n + 1)
      npq = n * r * q
      half_npq_m1 = 0.5 * npq - 1.0
      lg_m = Math.lgamma(m + 1).first
      lg_nm = Math.lgamma(n - m + 1).first
      log_rq = Math.log(r / q)

      loop do
        u = rng.rand * p4
        v = rng.rand

        if u <= p1
          return (xm - p1 * v + u).floor
        elsif u <= p2
          x = xl + (u - p1) / c
          v = v * c + 1.0 - (m - x + 0.5).abs / p1
          next if v > 1.0 || v <= 0.0

          ix = x.floor
        elsif u <= p3
          ix = (xl + Math.log(v) / lam_l).floor
          next if ix < 0

          v *= (u - p2) * lam_l
        else
          ix = (xr - Math.log(v) / lam_r).floor
          next if ix > n

          v *= (u - p3) * lam_r
        end

        k = (ix - m).abs
        if k <= 20
          f = 1.0
          if m < ix
            i = m + 1
            while i <= ix
              f *= (a / i - s)
              i += 1
            end
          elsif m > ix
            i = ix + 1
            while i <= m
              f /= (a / i - s)
              i += 1
            end
          end
          return ix if v <= f
        else
          if k < half_npq_m1
            amaxp = (k / npq) * ((k * (k / 3.0 + 0.625) + 1.0 / 6.0) / npq + 0.5)
            ynorm = -(k * k) / (2.0 * npq)
            alv = Math.log(v)
            return ix if alv < ynorm - amaxp
            next if alv > ynorm + amaxp
          end
          log_f = (lg_m + lg_nm) -
                  (Math.lgamma(ix + 1).first + Math.lgamma(n - ix + 1).first) +
                  (ix - m) * log_rq
          return ix if Math.log(v) <= log_f
        end
      end
    end

    def sample(n, p, rng, force_branch: nil)
      return 0 if n <= 0 || p <= 0.0
      return n if p >= 1.0

      flipped = p > 0.5
      r = flipped ? (1.0 - p) : p
      k = case force_branch
          when :inv
            sample_inv(n, r, rng)
          when :btpe
            sample_btpe(n, r, rng)
          else
            (n * r < 10.0) ? sample_inv(n, r, rng) : sample_btpe(n, r, rng)
          end
      flipped ? (n - k) : k
    end
  end
end
