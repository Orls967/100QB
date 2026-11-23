# frozen_string_literal: true

module FactorPortfolioLab
  # Equal Risk Contribution (Risk Parity) via Spinu (2013) Newton and Cyclic Coordinate Descent.
  module RiskParity
    module_function

    def solve_newton(sigma, budgets: nil, tol: 1.0e-12, max_iter: 80)
      n = sigma.length
      b = normalize_budgets(budgets, n)
      x = initial_guess(sigma, b, n)
      iters = 0
      converged = false

      max_iter.times do |iter|
        iters = iter + 1
        sx = LinAlg.mat_vec(sigma, x)
        grad = Array.new(n) { |i| sx[i] - (b[i] / x[i]) }
        if grad.map(&:abs).max < tol
          converged = true
          break
        end
        x = newton_line_search(sigma, b, x, sx, grad, n)
      end

      summarize_solution(sigma, x, iters, converged)
    end

    def solve_ccd(sigma, budgets: nil, tol: 1.0e-12, max_iter: 400)
      n = sigma.length
      b = normalize_budgets(budgets, n)
      x = initial_guess(sigma, b, n)
      sx = LinAlg.mat_vec(sigma, x)
      iters = 0
      converged = false

      max_iter.times do |iter|
        iters = iter + 1
        max_diff = ccd_sweep!(sigma, b, x, sx, n)
        if max_diff < tol
          converged = true
          break
        end
      end

      summarize_solution(sigma, x, iters, converged)
    end

    def ccd_sweep!(sigma, budgets, x, sx, n)
      max_diff = 0.0
      n.times do |i|
        aii = sigma[i][i]
        ci = sx[i] - (aii * x[i])
        disc = (ci * ci) + (4.0 * aii * budgets[i])
        x_new = (-ci + Math.sqrt([disc, 0.0].max)) / (2.0 * aii)
        delta = x_new - x[i]
        max_diff = delta.abs if delta.abs > max_diff
        x[i] = x_new
        n.times { |j| sx[j] += sigma[j][i] * delta }
      end
      max_diff
    end

    def newton_line_search(sigma, budgets, x, sx, grad, n)
      hess = Array.new(n) { |i| Array.new(n) { |j| sigma[i][j] } }
      n.times { |i| hess[i][i] += budgets[i] / (x[i] * x[i]) }
      neg_grad = grad.map(&:-@)
      step = LinAlg.solve_spd(hess, neg_grad)

      alpha = 1.0
      alpha *= 0.5 while n.times.any? { |i| x[i] + (alpha * step[i]) <= 1.0e-15 } && alpha > 1.0e-12
      f0 = spinu_obj(sigma, budgets, x, sx)
      g_dot_d = LinAlg.dot(grad, step)

      25.times do
        cand = Array.new(n) { |i| x[i] + (alpha * step[i]) }
        s_cand = LinAlg.mat_vec(sigma, cand)
        return cand if spinu_obj(sigma, budgets, cand, s_cand) <= f0 + (1.0e-4 * alpha * g_dot_d) + 1.0e-15

        alpha *= 0.5
      end
      Array.new(n) { |i| [x[i] + (alpha * step[i]), 1.0e-12].max }
    end

    def spinu_obj(sigma, budgets, x, sx = nil)
      sx ||= LinAlg.mat_vec(sigma, x)
      (0.5 * LinAlg.dot(x, sx)) - budgets.each_with_index.sum { |bi, i| bi * Math.log(x[i]) }
    end

    def normalize_budgets(budgets, n)
      return Array.new(n, 1.0 / n) if budgets.nil?

      total = budgets.sum.to_f
      budgets.map { |v| v / total }
    end

    def initial_guess(sigma, budgets, n)
      raw = Array.new(n) { |i| Math.sqrt(budgets[i] / [sigma[i][i], 1.0e-12].max) }
      total_var = LinAlg.quad_form(raw, sigma)
      scale = Math.sqrt(budgets.sum / [total_var, 1.0e-15].max)
      raw.map { |v| v * scale }
    end

    def summarize_solution(sigma, x, iters, converged)
      sum_x = x.sum
      w = x.map { |v| v / sum_x }
      sw = LinAlg.mat_vec(sigma, w)
      port_var = LinAlg.dot(w, sw)
      port_vol = Math.sqrt([port_var, 0.0].max)
      rc = w.each_with_index.map { |wi, i| (wi * sw[i]) / [port_vol, 1.0e-15].max }
      rel_rc = rc.map { |r| r / [port_vol, 1.0e-15].max }
      target = 1.0 / w.length
      rc_max_dev = rel_rc.map { |r| (r - target).abs }.max

      {
        weights: w,
        risk_contributions: rc,
        rel_risk_contributions: rel_rc,
        rc_max_dev: rc_max_dev,
        portfolio_vol: port_vol,
        iterations: iters,
        converged: converged
      }
    end
  end
end
