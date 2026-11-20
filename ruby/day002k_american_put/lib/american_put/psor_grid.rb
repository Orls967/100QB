# frozen_string_literal: true

module AmericanPut
  module PsorGrid
    module_function

    def build(params, n_x, n_t)
      x_max = [5.0 * params.sigma * Math.sqrt(params.t), 2.0].max
      x_min = -x_max
      step_dx = (x_max - x_min) / n_x.to_f
      step_dt = params.t / n_t.to_f
      xs = Array.new(n_x + 1) { |i| x_min + (i * step_dx) }
      spots = xs.map { |x_val| params.k * Math.exp(x_val) }
      payoff = spots.map { |spot| params.k > spot ? params.k - spot : 0.0 }
      {
        n_x: n_x, n_t: n_t, dx: step_dx, dt: step_dt, xs: xs, spots: spots, payoff: payoff,
        cn: theta_coefficients(params, step_dx, step_dt, 0.5),
        be: theta_coefficients(params, step_dx, 0.5 * step_dt, 1.0)
      }
    end

    def theta_coefficients(params, step_dx, step_dt, theta)
      diff = 0.5 * params.sigma * params.sigma
      inv_dx2 = diff / (step_dx * step_dx)
      adv = (params.r - diff) / (2.0 * step_dx)
      alpha = inv_dx2 - adv
      beta = (-2.0 * inv_dx2) - params.r
      gamma = inv_dx2 + adv
      rhs_w = 1.0 - theta
      {
        lhs_l: -theta * step_dt * alpha, lhs_d: 1.0 - (theta * step_dt * beta), lhs_u: -theta * step_dt * gamma,
        rhs_l: rhs_w * step_dt * alpha, rhs_d: 1.0 + (rhs_w * step_dt * beta), rhs_u: rhs_w * step_dt * gamma
      }
    end

    def build_rhs!(rhs, vals, coeffs, n_x)
      rhs_l = coeffs[:rhs_l]
      rhs_d = coeffs[:rhs_d]
      rhs_u = coeffs[:rhs_u]
      i = 1
      while i < n_x
        rhs[i] = (rhs_l * vals[i - 1]) + (rhs_d * vals[i]) + (rhs_u * vals[i + 1])
        i += 1
      end
    end

    def interpolate_cubic(x_grid, vals, x_target)
      idx = locate_interval(x_grid, x_target)
      i0 = (idx - 1).clamp(0, x_grid.length - 4)
      lagrange4(x_grid[i0, 4], vals[i0, 4], x_target)
    end

    def extract_greeks(x_grid, vals, x_target, spot0)
      idx = locate_interval(x_grid, x_target)
      i0 = (idx - 1).clamp(0, x_grid.length - 4)
      dv_dx, d2v_dx2 = lagrange4_derivatives(x_grid[i0, 4], vals[i0, 4], x_target)
      [dv_dx / spot0, (d2v_dx2 - dv_dx) / (spot0 * spot0)]
    end

    def locate_interval(x_grid, x_target)
      idx = x_grid.bsearch_index { |x_val| x_val >= x_target } || (x_grid.length - 1)
      (idx - 1).clamp(0, x_grid.length - 2)
    end

    def lagrange4(x_pts, y_pts, x_target)
      sum = 0.0
      4.times do |i|
        basis = 1.0
        4.times { |j| basis *= (x_target - x_pts[j]) / (x_pts[i] - x_pts[j]) if i != j }
        sum += y_pts[i] * basis
      end
      sum
    end

    def lagrange4_derivatives(x_pts, y_pts, x_target)
      eps = 1e-4 * (x_pts[1] - x_pts[0])
      vp = lagrange4(x_pts, y_pts, x_target + eps)
      v0 = lagrange4(x_pts, y_pts, x_target)
      vm = lagrange4(x_pts, y_pts, x_target - eps)
      [(vp - vm) / (2.0 * eps), (vp - (2.0 * v0) + vm) / (eps * eps)]
    end
  end
end
