# frozen_string_literal: true

module FactorPortfolioLab
  module RootFinder
    module_function

    def brent(lo, hi, tol: 1.0e-10, max_iter: 120, &block)
      a = lo.to_f
      b = hi.to_f
      fa = yield(a)
      fb = yield(b)
      raise ArgumentError, "Akar tidak diapit: f(#{a})=#{fa}, f(#{b})=#{fb}" if fa * fb > 0.0
      return a if fa.abs < tol
      return b if fb.abs < tol

      c = a
      fc = fa
      d = b - a
      e = d
      max_iter.times do
        a, b, c, fa, fb, fc, d, e = brent_step(a, b, c, fa, fb, fc, d, e, tol, &block)
        return b if fb.zero? || (0.5 * (c - b)).abs <= tol
      end
      b
    end

    def brent_step(a, b, c, fa, fb, fc, d, e, tol)
      if (fb.positive? && fc.positive?) || (fb.negative? && fc.negative?)
        c = a
        fc = fa
        d = b - a
        e = d
      end
      if fc.abs < fb.abs
        a, b, c = b, c, b
        fa, fb, fc = fb, fc, fb
      end
      tol1 = (2.0 * Float::EPSILON * b.abs) + (0.5 * tol)
      xm = 0.5 * (c - b)
      return [a, b, c, fa, fb, fc, d, e] if xm.abs <= tol1 || fb.zero?

      d, e = next_brent_delta(a, b, fa, fb, fc, d, e, xm, tol1)
      a = b
      fa = fb
      b += if d.abs > tol1
             d
           else
             (xm.positive? ? tol1 : -tol1)
           end
      fb = yield(b)
      [a, b, c, fa, fb, fc, d, e]
    end
    private_class_method :brent_step

    def next_brent_delta(a, b, fa, fb, fc, d, e, xm, tol1)
      if e.abs >= tol1 && fa.abs > fb.abs
        p_val, q_val = brent_interp(a, b, fa, fb, fc, xm)
        if p_val.positive?
          q_val = -q_val
        else
          p_val = -p_val
        end
        return [p_val / q_val, d] if (2.0 * p_val) < [((3.0 * xm * q_val) - (tol1 * q_val).abs), (e * q_val).abs].min
      end
      [xm, xm]
    end
    private_class_method :next_brent_delta

    def brent_interp(a, b, fa, fb, fc, xm)
      s = fb / fa
      if a == b
        p_val = 2.0 * xm * s
        q_val = 1.0 - s
      else
        q_r = fa / fc
        r = fb / fc
        p_val = s * ((2.0 * xm * q_r * (q_r - r)) - ((b - a) * (r - 1.0)))
        q_val = (q_r - 1.0) * (r - 1.0) * (s - 1.0)
      end
      p_val.positive? ? [p_val, -q_val] : [-p_val, q_val]
    end
    private_class_method :brent_interp

    def newton(x0, tol: 1.0e-10, max_iter: 60)
      x = x0.to_f
      max_iter.times do
        fx, dfx = yield(x)
        return x if fx.abs < tol
        break if dfx.abs < 1.0e-15

        x -= fx / dfx
      end
      x
    end
  end
end
