# frozen_string_literal: true

module Longevity
  module Linalg
    module_function

    def zeros(rows, cols)
      Array.new(rows) { Array.new(cols, 0.0) }
    end

    def eye(size)
      Array.new(size) { |i| Array.new(size) { |j| i == j ? 1.0 : 0.0 } }
    end

    def transpose(mat)
      rows = mat.size
      cols = mat[0].size
      out = zeros(cols, rows)
      rows.times do |i|
        row = mat[i]
        cols.times { |j| out[j][i] = row[j] }
      end
      out
    end

    def matmul(lhs, rhs)
      m = lhs.size
      k = lhs[0].size
      n = rhs[0].size
      rhs_t = transpose(rhs)
      out = zeros(m, n)
      m.times do |i|
        li = lhs[i]
        oi = out[i]
        n.times do |j|
          rj = rhs_t[j]
          s = 0.0
          k.times { |p| s += li[p] * rj[p] }
          oi[j] = s
        end
      end
      out
    end

    def frobenius_norm(mat)
      s = 0.0
      mat.each do |row|
        row.each { |v| s += v * v }
      end
      Math.sqrt(s)
    end

    def reconstruct_from_svd(u_mat, s_vals, v_mat)
      m = u_mat.size
      k = s_vals.size
      n = v_mat.size
      out = zeros(m, n)
      m.times do |i|
        ui = u_mat[i]
        oi = out[i]
        n.times do |j|
          vj = v_mat[j]
          sum = 0.0
          k.times { |p| sum += ui[p] * s_vals[p] * vj[p] }
          oi[j] = sum
        end
      end
      out
    end

    # Thin SVD via one-sided Hestenes-Jacobi rotations: A = U * diag(S) * V^T.
    # For m x n, returns { u: m x k, s: Array(k), v: n x k } where k = min(m, n).
    def svd(mat, max_sweeps: 30, tol: 1e-14)
      m = mat.size
      n = mat[0].size
      if m < n
        res = svd(transpose(mat), max_sweeps: max_sweeps, tol: tol)
        return { u: res[:v], s: res[:s], v: res[:u] }
      end

      b = mat.map(&:dup)
      v = eye(n)

      max_sweeps.times do
        converged = jacobi_sweep!(b, v, m, n, tol)
        break if converged
      end

      finalize_svd(b, v, m, n)
    end

    def jacobi_sweep!(mat_b, mat_v, rows, cols, tol)
      converged = true
      (0...(cols - 1)).each do |col_p|
        ((col_p + 1)...cols).each do |col_q|
          alpha, beta, gamma = pair_inner_products(mat_b, rows, col_p, col_q)
          denom = Math.sqrt(alpha * beta)
          next if denom <= 1e-30 || gamma.abs <= tol * denom

          converged = false
          zeta = (beta - alpha) / (2.0 * gamma)
          tan_t = (zeta >= 0.0 ? 1.0 : -1.0) / (zeta.abs + Math.sqrt(1.0 + (zeta * zeta)))
          cos_t = 1.0 / Math.sqrt(1.0 + (tan_t * tan_t))
          sin_t = cos_t * tan_t

          rotate_cols!(mat_b, rows, col_p, col_q, cos_t, sin_t)
          rotate_cols!(mat_v, cols, col_p, col_q, cos_t, sin_t)
        end
      end
      converged
    end

    def pair_inner_products(mat_b, rows, col_p, col_q)
      alpha = 0.0
      beta = 0.0
      gamma = 0.0
      rows.times do |i|
        bip = mat_b[i][col_p]
        biq = mat_b[i][col_q]
        alpha += bip * bip
        beta += biq * biq
        gamma += bip * biq
      end
      [alpha, beta, gamma]
    end

    def rotate_cols!(mat, rows, col_p, col_q, cos_t, sin_t)
      rows.times do |i|
        row = mat[i]
        xp = row[col_p]
        xq = row[col_q]
        row[col_p] = (cos_t * xp) - (sin_t * xq)
        row[col_q] = (sin_t * xp) + (cos_t * xq)
      end
    end

    def finalize_svd(mat_b, mat_v, rows, cols)
      norms = Array.new(cols) do |j|
        s = 0.0
        rows.times { |i| s += mat_b[i][j] * mat_b[i][j] }
        Math.sqrt(s)
      end

      order = (0...cols).sort_by { |j| -norms[j] }
      s_vals = order.map { |j| norms[j] }

      u_out = zeros(rows, cols)
      v_out = zeros(cols, cols)
      order.each_with_index do |orig_j, new_j|
        sig = norms[orig_j]
        if sig > 1e-14
          inv_sig = 1.0 / sig
          rows.times { |i| u_out[i][new_j] = mat_b[i][orig_j] * inv_sig }
        end
        cols.times { |i| v_out[i][new_j] = mat_v[i][orig_j] }
      end

      orthonormalize_zero_cols!(u_out, s_vals, rows, cols)
      { u: u_out, s: s_vals, v: v_out }
    end

    def orthonormalize_zero_cols!(u_out, s_vals, rows, cols)
      cols.times do |j|
        next if s_vals[j] > 1e-14

        vec = Array.new(rows) { |i| i == j ? 1.0 : 0.0 }
        j.times do |k|
          dot = 0.0
          rows.times { |i| dot += vec[i] * u_out[i][k] }
          rows.times { |i| vec[i] -= dot * u_out[i][k] }
        end
        nrm = Math.sqrt(vec.sum(0.0) { |x| x * x })
        rows.times { |i| u_out[i][j] = vec[i] / nrm } if nrm > 1e-14
      end
    end
  end
end
