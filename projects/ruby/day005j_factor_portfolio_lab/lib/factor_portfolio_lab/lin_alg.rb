# frozen_string_literal: true

module FactorPortfolioLab
  module LinAlg
    module_function

    def dot(a, b)
      s = 0.0
      a.length.times { |i| s += a[i] * b[i] }
      s
    end

    def mat_vec(mat, vec) = mat.map { |row| dot(row, vec) }

    def quad_form(vec, mat)
      n = vec.length
      n.times.sum(0.0) { |i| vec[i] * dot(mat[i], vec) }
    end

    def transpose(mat)
      rows = mat.length
      cols = mat.first.length
      Array.new(cols) { |j| Array.new(rows) { |i| mat[i][j] } }
    end

    def mat_mul(a, b)
      rows = a.length
      cols = b.first.length
      bt = transpose(b)
      Array.new(rows) { |i| Array.new(cols) { |j| dot(a[i], bt[j]) } }
    end

    def trace(mat) = mat.each_with_index.sum(0.0) { |row, i| row[i] }

    def eye(n, scale = 1.0) = Array.new(n) { |i| Array.new(n) { |j| i == j ? scale : 0.0 } }

    def solve(mat, rhs)
      n = rhs.length
      a = mat.map(&:dup)
      b = rhs.dup
      n.times do |col|
        piv = (col...n).max_by { |r| a[r][col].abs }
        a[col], a[piv] = a[piv], a[col]
        b[col], b[piv] = b[piv], b[col]
        raise 'Matriks singular' if a[col][col].abs < 1.0e-16

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
      transpose(cols)
    end

    def cholesky(mat)
      n = mat.length
      l = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        (0..i).each do |j|
          s = (0...j).sum(0.0) { |k| l[i][k] * l[j][k] }
          val = i == j ? mat[i][i] - s : (mat[i][j] - s) / l[j][j]
          raise 'Matriks tidak definit positif' if i == j && val <= 1.0e-15

          l[i][j] = i == j ? Math.sqrt(val) : val
        end
      end
      l
    end

    def solve_spd(mat, rhs)
      l = cholesky(mat)
      n = rhs.length
      y = Array.new(n, 0.0)
      n.times { |i| y[i] = (rhs[i] - (0...i).sum(0.0) { |k| l[i][k] * y[k] }) / l[i][i] }
      x = Array.new(n, 0.0)
      (n - 1).downto(0) { |i| x[i] = (y[i] - ((i + 1)...n).sum(0.0) { |k| l[k][i] * x[k] }) / l[i][i] }
      x
    end

    def lu_det(mat)
      n = mat.length
      a = mat.map(&:dup)
      sign = 1.0
      n.times do |col|
        piv = (col...n).max_by { |r| a[r][col].abs }
        if piv != col
          a[col], a[piv] = a[piv], a[col]
          sign = -sign
        end
        return 0.0 if a[col][col].abs < 1.0e-18

        ((col + 1)...n).each do |r|
          fac = a[r][col] / a[col][col]
          (col...n).each { |c| a[r][c] -= fac * a[col][c] }
        end
      end
      sign * n.times.reduce(1.0) { |acc, i| acc * a[i][i] }
    end

    def qr_householder(mat)
      m = mat.length
      n = mat.first.length
      r = mat.map(&:dup)
      q = eye(m)
      [m - 1, n].min.times { |k| qr_step!(q, r, m, n, k) }
      [q, r]
    end

    def qr_step!(q, r, m, n, k)
      norm_x = Math.sqrt((k...m).sum(0.0) { |i| r[i][k]**2 })
      return if norm_x < 1.0e-16

      v = Array.new(m - k) { |i| r[k + i][k] }
      v[0] += v[0] >= 0.0 ? norm_x : -norm_x
      v_norm2 = dot(v, v)
      return if v_norm2 < 1.0e-24

      apply_hh_r!(r, v, v_norm2, m, n, k)
      apply_hh_q!(q, v, v_norm2, m, k)
    end
    private_class_method :qr_step!

    def apply_hh_r!(r, v, v_norm2, m, n, k)
      (k...n).each do |j|
        s = 2.0 * (0...(m - k)).sum(0.0) { |i| v[i] * r[k + i][j] } / v_norm2
        (0...(m - k)).each { |i| r[k + i][j] -= s * v[i] }
      end
    end
    private_class_method :apply_hh_r!

    def apply_hh_q!(q, v, v_norm2, m, k)
      m.times do |i|
        s = 2.0 * (0...(m - k)).sum(0.0) { |j| q[i][k + j] * v[j] } / v_norm2
        (0...(m - k)).each { |j| q[i][k + j] -= s * v[j] }
      end
    end
    private_class_method :apply_hh_q!

    def jacobi_sym(mat, max_sweeps: 60, tol: 1.0e-12)
      n = mat.length
      a = Array.new(n) { |i| Array.new(n) { |j| 0.5 * (mat[i][j] + mat[j][i]) } }
      v = eye(n)
      max_sweeps.times do
        max_off = jacobi_sweep!(a, v, n, tol)
        break if max_off < tol
      end
      order = (0...n).sort_by { |i| -a[i][i] }
      { eigenvalues: order.map { |i| a[i][i] }, eigenvectors: Array.new(n) { |r| order.map { |c| v[r][c] } } }
    end

    def jacobi_sweep!(a, v, n, tol)
      max_off = 0.0
      n.times do |p|
        ((p + 1)...n).each do |q|
          apq = a[p][q].abs
          max_off = apq if apq > max_off
          jacobi_rotate!(a, v, n, p, q) if apq > tol
        end
      end
      max_off
    end
    private_class_method :jacobi_sweep!

    def jacobi_rotate!(a, v, n, p, q)
      tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
      t = tau >= 0.0 ? 1.0 / (tau + Math.sqrt(1.0 + (tau * tau))) : -1.0 / (-tau + Math.sqrt(1.0 + (tau * tau)))
      c = 1.0 / Math.sqrt(1.0 + (t * t))
      s = t * c
      jacobi_update_entries!(a, v, n, p, q, c, s, t)
    end
    private_class_method :jacobi_rotate!

    def jacobi_update_entries!(a, v, n, p, q, c, s, t)
      apq = a[p][q]
      a[p][p] -= t * apq
      a[q][q] += t * apq
      a[p][q] = a[q][p] = 0.0
      n.times do |r|
        next if r == p || r == q

        arp = a[r][p]
        arq = a[r][q]
        a[r][p] = a[p][r] = (c * arp) - (s * arq)
        a[r][q] = a[q][r] = (s * arp) + (c * arq)
      end
      n.times do |r|
        vrp = v[r][p]
        vrq = v[r][q]
        v[r][p] = (c * vrp) - (s * vrq)
        v[r][q] = (s * vrp) + (c * vrq)
      end
    end
    private_class_method :jacobi_update_entries!

    def cov_matrix(rows, ddof: 1)
      t = rows.length
      n = rows.first.length
      means = Array.new(n) { |j| rows.sum(0.0) { |r| r[j] } / t }
      cov = Array.new(n) { Array.new(n, 0.0) }
      rows.each do |r|
        n.times do |i|
          di = r[i] - means[i]
          (0..i).each { |j| cov[i][j] += di * (r[j] - means[j]) }
        end
      end
      denom = (t - ddof).to_f
      n.times { |i| (0..i).each { |j| cov[j][i] = cov[i][j] = cov[i][j] / denom } }
      cov
    end

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
