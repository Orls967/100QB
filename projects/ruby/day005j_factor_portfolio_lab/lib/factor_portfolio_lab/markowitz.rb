# frozen_string_literal: true

module FactorPortfolioLab
  # Mean-Variance (Markowitz) optimization: analytical unconstrained and custom KKT QP solver.
  module Markowitz
    module_function

    def unconstrained_gmv(sigma)
      n = sigma.length
      ones = Array.new(n, 1.0)
      inv_ones = LinAlg.solve_spd(sigma, ones)
      denom = inv_ones.sum
      inv_ones.map { |v| v / denom }
    end

    def unconstrained_tangency(sigma, mu)
      inv_mu = LinAlg.solve_spd(sigma, mu)
      denom = inv_mu.sum
      inv_mu.map { |v| v / denom }
    end

    def unconstrained_frontier_point(sigma, mu, target_return)
      n = sigma.length
      ones = Array.new(n, 1.0)
      s_inv_1 = LinAlg.solve_spd(sigma, ones)
      s_inv_mu = LinAlg.solve_spd(sigma, mu)
      a = LinAlg.dot(ones, s_inv_mu)
      b = LinAlg.dot(mu, s_inv_mu)
      c = LinAlg.dot(ones, s_inv_1)
      d = (b * c) - (a * a)
      lam1 = ((b - (a * target_return)) / d)
      lam2 = (((c * target_return) - a) / d)
      Array.new(n) { |i| (lam1 * s_inv_1[i]) + (lam2 * s_inv_mu[i]) }
    end

    # Solves min 0.5 * w' Sigma w - kappa * mu' w s.t. sum(w) = 1, lower <= w <= upper.
    def solve_qp(sigma, mu: nil, kappa: 0.0, lower: nil, upper: nil, tol: 1.0e-10, max_iter: 2000)
      n = sigma.length
      lo = lower || Array.new(n, 0.0)
      hi = upper || Array.new(n, 1.0)
      mu_vec = mu || Array.new(n, 0.0)

      w = initial_feasible_weights(lo, hi, n)
      lip = power_spectral_radius(sigma, n)
      step = 1.0 / [lip, 1.0e-6].max

      w = projected_bb_loop(sigma, mu_vec, kappa, w, lo, hi, step, tol, max_iter)
      w = active_set_polish(sigma, mu_vec, kappa, w, lo, hi, tol)

      grad = qp_grad(sigma, mu_vec, kappa, w)
      kkt = Constraints.kkt_residual(w, grad, lower: lo, upper: hi)
      port_var = LinAlg.quad_form(w, sigma)

      {
        weights: w,
        expected_return: LinAlg.dot(w, mu_vec),
        volatility: Math.sqrt([port_var, 0.0].max),
        kkt: kkt
      }
    end

    def projected_bb_loop(sigma, mu_vec, kappa, w, lo, hi, base_step, tol, max_iter)
      step = base_step
      g = qp_grad(sigma, mu_vec, kappa, w)
      max_iter.times do |iter|
        cand = Array.new(w.length) { |i| w[i] - (step * g[i]) }
        w_next = Constraints.project_capped_simplex(cand, lower: lo, upper: hi)
        s = Array.new(w.length) { |i| w_next[i] - w[i] }
        g_next = qp_grad(sigma, mu_vec, kappa, w_next)
        y = Array.new(w.length) { |i| g_next[i] - g[i] }

        if (iter % 25).zero?
          kkt = Constraints.kkt_residual(w_next, g_next, lower: lo, upper: hi)
          return w_next if kkt[:max_kkt] < tol
        end

        s_dot_y = LinAlg.dot(s, y)
        if s_dot_y > 1.0e-16
          step = iter.even? ? LinAlg.dot(s, s) / s_dot_y : s_dot_y / [LinAlg.dot(y, y), 1.0e-20].max
          step = Constraints.clamp(step, base_step * 0.05, base_step * 50.0)
        else
          step = base_step
        end
        w = w_next
        g = g_next
        break if s.map(&:abs).max < 1.0e-14
      end
      w
    end

    def active_set_polish(sigma, mu_vec, kappa, w, lo, hi, tol)
      n = w.length
      15.times do
        g = qp_grad(sigma, mu_vec, kappa, w)
        kkt = Constraints.kkt_residual(w, g, lower: lo, upper: hi)
        return w if kkt[:max_kkt] < [tol * 0.1, 1.0e-11].min

        free = (0...n).select { |i| w[i] > lo[i] + 1.0e-9 && w[i] < hi[i] - 1.0e-9 }
        break if free.length < 2

        w_try = solve_free_subproblem(sigma, mu_vec, kappa, w, free)
        break if w_try.nil?

        feasible = free.all? { |i| w_try[i].between?(lo[i] - 1.0e-12, hi[i] + 1.0e-12) }
        if feasible
          free.each { |i| w[i] = Constraints.clamp(w_try[i], lo[i], hi[i]) }
        else
          alpha = max_feasible_step(w, w_try, free, lo, hi)
          free.each { |i| w[i] = Constraints.clamp(w[i] + (alpha * (w_try[i] - w[i])), lo[i], hi[i]) }
        end
      end
      w
    end

    def solve_free_subproblem(sigma, mu_vec, kappa, w, free)
      m = free.length
      fixed = (0...w.length).reject { |i| free.include?(i) }
      target_sum = 1.0 - fixed.sum { |j| w[j] }

      s_ff = Array.new(m) { |a| Array.new(m) { |b| sigma[free[a]][free[b]] } }
      rhs = Array.new(m) do |a|
        i = free[a]
        fixed_term = fixed.sum { |j| sigma[i][j] * w[j] }
        (kappa * mu_vec[i]) - fixed_term
      end

      ones = Array.new(m, 1.0)
      inv_rhs = LinAlg.solve_spd(s_ff, rhs)
      inv_ones = LinAlg.solve_spd(s_ff, ones)
      nu = (target_sum - inv_rhs.sum) / [inv_ones.sum, 1.0e-15].max

      w_out = w.dup
      m.times { |a| w_out[free[a]] = inv_rhs[a] + (nu * inv_ones[a]) }
      w_out
    rescue StandardError
      nil
    end

    def max_feasible_step(w, w_try, free, lo, hi)
      alpha = 1.0
      free.each do |i|
        d = w_try[i] - w[i]
        if d < -1.0e-14
          alpha = [alpha, (lo[i] - w[i]) / d].min
        elsif d > 1.0e-14
          alpha = [alpha, (hi[i] - w[i]) / d].min
        end
      end
      [alpha, 0.0].max
    end

    def qp_grad(sigma, mu_vec, kappa, w)
      sw = LinAlg.mat_vec(sigma, w)
      sw.each_with_index.map { |v, i| v - (kappa * mu_vec[i]) }
    end

    def initial_feasible_weights(lo, hi, n)
      eq = Array.new(n, 1.0 / n)
      Constraints.project_capped_simplex(eq, lower: lo, upper: hi)
    end

    def power_spectral_radius(sigma, n)
      v = Array.new(n, 1.0 / Math.sqrt(n))
      12.times do
        sv = LinAlg.mat_vec(sigma, v)
        norm = Math.sqrt(LinAlg.dot(sv, sv))
        break if norm < 1.0e-15

        v = sv.map { |x| x / norm }
      end
      LinAlg.quad_form(v, sigma)
    end
  end
end
