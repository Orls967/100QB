# frozen_string_literal: true

module GarchFamily
  module Linalg
    module_function

    def solve(mat, rhs)
      n = rhs.length
      a = mat.map(&:dup)
      b = rhs.dup
      n.times do |col|
        piv = (col...n).max_by { |r| a[r][col].abs }
        a[col], a[piv] = a[piv], a[col]
        b[col], b[piv] = b[piv], b[col]
        raise 'Matriks singular' if a[col][col].abs < 1.0e-18

        eliminate_col!(a, b, col, n)
      end
      back_sub(a, b, n)
    end

    def invert(mat)
      n = mat.length
      cols = n.times.map do |j|
        rhs = Array.new(n, 0.0)
        rhs[j] = 1.0
        solve(mat, rhs)
      end
      Array.new(n) { |i| Array.new(n) { |j| cols[j][i] } }
    end

    def mat_mul(a, b)
      rows = a.length
      cols = b.first.length
      inner = b.length
      Array.new(rows) do |i|
        Array.new(cols) { |j| inner.times.sum(0.0) { |k| a[i][k] * b[k][j] } }
      end
    end

    def dot(a, b)
      s = 0.0
      a.length.times { |i| s += a[i] * b[i] }
      s
    end

    def quad_form(vec, mat)
      n = vec.length
      n.times.sum(0.0) { |i| n.times.sum(0.0) { |j| vec[i] * mat[i][j] * vec[j] } }
    end

    def ols_with_nw(x_mat, y_vec, lags: nil)
      n = y_vec.length
      p = x_mat.first.length
      xtx, xty = accumulate_normal_eq(x_mat, y_vec, n, p)
      beta = solve(xtx, xty)
      resid = Array.new(n) { |t| y_vec[t] - dot(x_mat[t], beta) }
      bread = invert(xtx)
      l_max = lags || (4.0 * ((n / 100.0)**(2.0 / 9.0))).floor.clamp(1, 12)
      meat = newey_west_meat(x_mat, resid, n, p, l_max)
      cov = mat_mul(mat_mul(bread, meat), bread)
      se = Array.new(p) { |i| Math.sqrt([cov[i][i], 0.0].max) }
      { beta: beta, residuals: resid, se_nw: se, cov_nw: cov }
    end

    def accumulate_normal_eq(x_mat, y_vec, n, p)
      xtx = Array.new(p) { Array.new(p, 0.0) }
      xty = Array.new(p, 0.0)
      n.times do |t|
        row = x_mat[t]
        yv = y_vec[t]
        p.times do |i|
          xi = row[i]
          xty[i] += xi * yv
          p.times { |j| xtx[i][j] += xi * row[j] }
        end
      end
      [xtx, xty]
    end
    private_class_method :accumulate_normal_eq

    def newey_west_meat(x_mat, resid, n, p, l_max)
      meat = Array.new(p) { Array.new(p, 0.0) }
      add_lag_cov!(meat, x_mat, resid, n, p, 0, 1.0)
      (1..l_max).each { |l| add_lag_cov!(meat, x_mat, resid, n, p, l, 1.0 - (l.to_f / (l_max + 1.0))) }
      meat
    end
    private_class_method :newey_west_meat

    def add_lag_cov!(meat, x_mat, resid, n, p, lag, weight)
      (lag...n).each do |t|
        gt = x_mat[t]
        gl = x_mat[t - lag]
        sc = weight * resid[t] * resid[t - lag]
        p.times do |i|
          p.times do |j|
            meat[i][j] += sc * gt[i] * gl[j]
            meat[i][j] += sc * gl[i] * gt[j] if lag.positive?
          end
        end
      end
    end
    private_class_method :add_lag_cov!

    def eliminate_col!(a, b, col, n)
      diag = a[col][col]
      ((col + 1)...n).each do |r|
        fac = a[r][col] / diag
        (col...n).each { |c| a[r][c] -= fac * a[col][c] }
        b[r] -= fac * b[col]
      end
    end
    private_class_method :eliminate_col!

    def back_sub(a, b, n)
      x = Array.new(n, 0.0)
      (n - 1).downto(0) do |i|
        sum = b[i]
        ((i + 1)...n).each { |j| sum -= a[i][j] * x[j] }
        x[i] = sum / a[i][i]
      end
      x
    end
    private_class_method :back_sub
  end
end
