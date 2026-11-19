# frozen_string_literal: true

module AmericanPut
  module Psor
    DEFAULT_OMEGA = 1.25
    DEFAULT_TOL = 1e-9
    DEFAULT_MAX_ITERS = 5_000

    module_function

    def price(params, n_x: DEFAULT_PSOR_NX, n_t: DEFAULT_PSOR_NT,
              omega: DEFAULT_OMEGA, tol: DEFAULT_TOL, max_iters: DEFAULT_MAX_ITERS)
      params.validate!
      raise ArgumentError, 'n_x harus >= 10 dan n_t harus >= 2' if n_x < 10 || n_t < 2

      grid = PsorGrid.build(params, n_x, n_t)
      am_v, avg_sweeps = solve_american_rannacher(grid, omega, tol, max_iters)
      eu_v = solve_european_rannacher(grid, params)

      x_target = Math.log(params.s0 / params.k)
      am_price = PsorGrid.interpolate_cubic(grid[:xs], am_v, x_target)
      eu_price = PsorGrid.interpolate_cubic(grid[:xs], eu_v, x_target)
      delta, gamma = PsorGrid.extract_greeks(grid[:xs], am_v, x_target, params.s0)

      PsorResult.new(
        n_x: n_x, n_t: n_t, american_price: am_price, european_price: eu_price,
        early_exercise_premium: am_price - eu_price, delta: delta, gamma: gamma, avg_sweeps: avg_sweeps
      )
    end

    def convergence_sweep(params, grid_sizes: PSOR_GRID, workers: nil)
      Parallel.map_tasks(grid_sizes, workers: workers) do |(n_x, n_t), _idx|
        price(params, n_x: n_x, n_t: n_t)
      end
    end

    def solve_american_rannacher(grid, omega, tol, max_iters)
      vals = grid[:payoff].dup
      v_next = vals.dup
      rhs = Array.new(grid[:n_x] + 1, 0.0)
      total_sweeps = 0

      # Rannacher smoothing: 4 setengah-langkah Backward Euler untuk meredam osilasi di strike.
      4.times do
        vals, v_next, sweeps = advance_american_step(vals, v_next, rhs, grid, grid[:be], omega, tol, max_iters)
        total_sweeps += sweeps
      end
      (grid[:n_t] - 2).times do
        vals, v_next, sweeps = advance_american_step(vals, v_next, rhs, grid, grid[:cn], omega, tol, max_iters)
        total_sweeps += sweeps
      end
      [vals, total_sweeps.to_f / (grid[:n_t] + 2)]
    end
    private_class_method :solve_american_rannacher

    def advance_american_step(vals, v_next, rhs, grid, coeffs, omega, tol, max_iters)
      n_x = grid[:n_x]
      PsorGrid.build_rhs!(rhs, vals, coeffs, n_x)
      v_next[0] = grid[:payoff][0]
      v_next[n_x] = 0.0
      rhs[1] -= coeffs[:lhs_l] * v_next[0]
      rhs[n_x - 1] -= coeffs[:lhs_u] * v_next[n_x]
      sweeps = psor_step!(v_next, rhs, grid[:payoff], coeffs, n_x, omega, tol, max_iters)
      [v_next, vals, sweeps]
    end
    private_class_method :advance_american_step

    def psor_step!(v_next, rhs, payoff, coeffs, n_x, omega, tol, max_iters)
      sweeps = 0
      while sweeps < max_iters
        sweeps += 1
        max_diff = psor_sweep!(v_next, rhs, payoff, coeffs, n_x, omega)
        break if max_diff < tol
      end
      sweeps
    end
    private_class_method :psor_step!

    def psor_sweep!(v_next, rhs, payoff, coeffs, n_x, omega)
      lhs_l = coeffs[:lhs_l]
      lhs_d = coeffs[:lhs_d]
      lhs_u = coeffs[:lhs_u]
      max_diff = 0.0
      i = 1
      while i < n_x
        r_i = rhs[i] - (lhs_l * v_next[i - 1]) - (lhs_u * v_next[i + 1])
        cand = [v_next[i] + (omega * ((r_i / lhs_d) - v_next[i])), payoff[i]].max
        diff = (cand - v_next[i]).abs
        max_diff = diff if diff > max_diff
        v_next[i] = cand
        i += 1
      end
      max_diff
    end
    private_class_method :psor_sweep!

    def solve_european_rannacher(grid, params)
      vals = grid[:payoff].dup
      rhs_full = Array.new(grid[:n_x] + 1, 0.0)
      tau = 0.0
      4.times do
        tau += 0.5 * grid[:dt]
        vals = advance_european_step(vals, rhs_full, grid, grid[:be], params, tau)
      end
      (grid[:n_t] - 2).times do
        tau += grid[:dt]
        vals = advance_european_step(vals, rhs_full, grid, grid[:cn], params, tau)
      end
      vals
    end
    private_class_method :solve_european_rannacher

    def advance_european_step(vals, rhs_full, grid, coeffs, params, tau)
      n_x = grid[:n_x]
      PsorGrid.build_rhs!(rhs_full, vals, coeffs, n_x)
      left_bc = [(params.k * Math.exp(-params.r * tau)) - grid[:spots][0], 0.0].max
      rhs_int = rhs_full[1...n_x]
      rhs_int[0] -= coeffs[:lhs_l] * left_bc
      m = n_x - 1
      sol = LinearSolver.solve_tridiagonal(
        Array.new(m, coeffs[:lhs_l]), Array.new(m, coeffs[:lhs_d]), Array.new(m, coeffs[:lhs_u]), rhs_int
      )
      [left_bc] + sol + [0.0]
    end
    private_class_method :advance_european_step
  end
end
