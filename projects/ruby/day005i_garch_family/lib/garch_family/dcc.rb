# frozen_string_literal: true

module GarchFamily
  module Dcc
    TRUE_DCC = { a: 0.05, b: 0.90 }.freeze
    M1_PARAMS = { mu: 0.0, omega: 4.0e-6, alpha: 0.07, beta: 0.87 }.freeze
    M2_PARAMS = { mu: 0.0, omega: 6.0e-6, alpha: 0.09, beta: 0.83 }.freeze
    RHO_BAR = 0.45

    module_function

    def sample_corr_matrix(z_rows)
      t = z_rows.length.to_f
      k = z_rows.first.length
      q = Array.new(k) { Array.new(k, 0.0) }
      z_rows.each { |z| k.times { |i| k.times { |j| q[i][j] += z[i] * z[j] / t } } }
      normalize_corr(q)
    end

    def normalize_corr(q)
      k = q.length
      inv_sd = Array.new(k) { |i| 1.0 / Math.sqrt([q[i][i], 1.0e-15].max) }
      Array.new(k) { |i| Array.new(k) { |j| i == j ? 1.0 : (q[i][j] * inv_sd[i] * inv_sd[j]).clamp(-0.9999, 0.9999) } }
    end

    def cholesky_lower(mat)
      n = mat.length
      l = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        (0..i).each do |j|
          s = (0...j).sum(0.0) { |k| l[i][k] * l[j][k] }
          val = i == j ? mat[i][i] - s : (mat[i][j] - s) / l[j][j]
          return nil if i == j && val <= 1.0e-14

          l[i][j] = i == j ? Math.sqrt(val) : val
        end
      end
      l
    end

    def log_det_and_quad(mat, vec)
      l = cholesky_lower(mat)
      return nil if l.nil?

      n = mat.length
      y = Array.new(n, 0.0)
      logdet = 0.0
      n.times do |i|
        s = (0...i).sum(0.0) { |k| l[i][k] * y[k] }
        y[i] = (vec[i] - s) / l[i][i]
        logdet += Math.log(l[i][i])
      end
      [(2.0 * logdet), y.sum(0.0) { |v| v * v }]
    end

    def filter_corr(a_par, b_par, z_rows, q_bar: nil)
      qb = q_bar || sample_corr_matrix(z_rows)
      k = qb.length
      c = 1.0 - a_par - b_par
      q = qb.map(&:dup)
      r_seq = Array.new(z_rows.length)
      z_rows.each_with_index do |_z, t|
        if t.positive?
          zp = z_rows[t - 1]
          q = Array.new(k) { |i| Array.new(k) { |j| (c * qb[i][j]) + (a_par * zp[i] * zp[j]) + (b_par * q[i][j]) } }
        end
        r_seq[t] = normalize_corr(q)
      end
      { q_bar: qb, corr_seq: r_seq }
    end

    def corr_loglik(a_par, b_par, z_rows, q_bar: nil)
      r_seq = filter_corr(a_par, b_par, z_rows, q_bar: q_bar)[:corr_seq]
      ll = 0.0
      z_rows.each_with_index do |z, t|
        res = log_det_and_quad(r_seq[t], z)
        return -Float::INFINITY if res.nil?

        ll -= 0.5 * (res[0] + res[1] - Linalg.dot(z, z))
      end
      ll
    end

    def untransform_ab(theta)
      lim = 0.998
      ea = Math.exp(theta[0].clamp(-25.0, 25.0))
      eb = Math.exp(theta[1].clamp(-25.0, 25.0))
      norm = 1.0 + ea + eb
      [lim * ea / norm, lim * eb / norm]
    end

    def transform_ab(a_par, b_par)
      lim = 0.998
      s = [1.0 - ((a_par + b_par) / lim), 1.0e-5].max
      [Math.log([a_par / (lim * s), 1.0e-5].max), Math.log([b_par / (lim * s), 1.0e-5].max)]
    end

    def fit_stage2(z_rows, init_a: 0.04, init_b: 0.90)
      qb = sample_corr_matrix(z_rows)
      th0 = transform_ab(init_a, init_b)
      inv_n = 1.0 / z_rows.length
      obj = ->(th) { stage2_obj(th, z_rows, qb, inv_n) }
      opt = Optim.bfgs(th0, obj, max_iter: 80, gtol: 1.0e-6)
      a_hat, b_hat = untransform_ab(opt[:x])
      { a: a_hat, b: b_hat, persistence: a_hat + b_hat, loglik: -opt[:f] * z_rows.length,
        q_bar: qb, converged: opt[:converged] }
    end

    def stage2_obj(th, z_rows, qb, inv_n)
      a, b = untransform_ab(th)
      ll = corr_loglik(a, b, z_rows, q_bar: qb)
      return [1.0e30, [0.0, 0.0]] unless ll.finite?

      h = 1.0e-5
      g = Array.new(2) do |j|
        tp = th.dup
        tm = th.dup
        tp[j] += h
        tm[j] -= h
        ap, bp = untransform_ab(tp)
        am, bm = untransform_ab(tm)
        -(corr_loglik(ap, bp, z_rows, q_bar: qb) - corr_loglik(am, bm, z_rows, q_bar: qb)) / (2.0 * h) * inv_n
      end
      [-ll * inv_n, g]
    end
    private_class_method :stage2_obj
  end
end
