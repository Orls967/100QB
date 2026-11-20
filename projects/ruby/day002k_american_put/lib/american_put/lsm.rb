# frozen_string_literal: true

module AmericanPut
  module Lsm
    module_function

    def price(params, n_paths: DEFAULT_LSM_PATHS, n_steps: DEFAULT_LSM_STEPS,
              degree: 3, seed: DEFAULT_SEED, antithetic: true, workers: nil)
      params.validate!
      raise ArgumentError, 'n_paths harus >= 4 dan genap bila antithetic' if n_paths < 4 || (antithetic && n_paths.odd?)
      raise ArgumentError, 'n_steps harus >= 2' if n_steps < 2

      step_dt = params.t / n_steps.to_f
      disc = Math.exp(-params.r * step_dt)
      spots_by_step = LsmSimulator.simulate_all_paths(params, n_paths, n_steps, step_dt, seed, antithetic, workers)
      cashflows = spots_by_step[n_steps - 1].map { |spot| params.k > spot ? params.k - spot : 0.0 }

      backward_regress!(spots_by_step, cashflows, params.k, disc, n_steps, n_paths, degree)
      summarize_lsm(cashflows, params, n_paths, n_steps, degree, antithetic, disc)
    end

    def convergence_sweep(params, paths_grid: LSM_PATHS_GRID, n_steps: DEFAULT_LSM_STEPS,
                          degree: 3, seed: DEFAULT_SEED, antithetic: true, workers: nil)
      paths_grid.map do |n_paths|
        price(params, n_paths: n_paths, n_steps: n_steps, degree: degree,
                      seed: seed, antithetic: antithetic, workers: workers)
      end
    end

    def backward_regress!(spots_by_step, cashflows, strike, disc, n_steps, n_paths, degree)
      step = n_steps - 2
      while step >= 0
        discount_in_place!(cashflows, disc, n_paths)
        spots = spots_by_step[step]
        beta = fit_continuation_poly(spots, cashflows, strike, n_paths, degree)
        apply_exercise_decision!(spots, cashflows, strike, beta, n_paths, degree) if beta
        step -= 1
      end
    end
    private_class_method :backward_regress!

    def discount_in_place!(cashflows, disc, n_paths)
      i = 0
      while i < n_paths
        cashflows[i] *= disc
        i += 1
      end
    end
    private_class_method :discount_in_place!

    def fit_continuation_poly(spots, cashflows, strike, n_paths, degree)
      basis_count = degree + 1
      power_sums = Array.new((2 * degree) + 1, 0.0)
      rhs = Array.new(basis_count, 0.0)
      itm_count = accumulate_normal_eqs!(power_sums, rhs, spots, cashflows, strike, n_paths, degree)
      return nil if itm_count <= basis_count

      gram = Array.new(basis_count) { |row| Array.new(basis_count) { |col| power_sums[row + col] } }
      LinearSolver.solve(gram, rhs)
    rescue ZeroDivisionError
      nil
    end
    private_class_method :fit_continuation_poly

    def accumulate_normal_eqs!(power_sums, rhs, spots, cashflows, strike, n_paths, degree)
      itm_count = 0
      max_pow = 2 * degree
      i = 0
      while i < n_paths
        spot = spots[i]
        if spot < strike
          itm_count += 1
          add_itm_observation!(power_sums, rhs, spot / strike, cashflows[i], degree, max_pow)
        end
        i += 1
      end
      itm_count
    end
    private_class_method :accumulate_normal_eqs!

    def add_itm_observation!(power_sums, rhs, scaled_x, cashflow, degree, max_pow)
      x_pow = 1.0
      pow = 0
      while pow <= max_pow
        power_sums[pow] += x_pow
        rhs[pow] += x_pow * cashflow if pow <= degree
        x_pow *= scaled_x
        pow += 1
      end
    end
    private_class_method :add_itm_observation!

    def apply_exercise_decision!(spots, cashflows, strike, beta, n_paths, degree)
      i = 0
      while i < n_paths
        spot = spots[i]
        if spot < strike
          cont = eval_horner(beta, spot / strike, degree)
          exercise = strike - spot
          cashflows[i] = exercise if exercise > cont
        end
        i += 1
      end
    end
    private_class_method :apply_exercise_decision!

    def eval_horner(beta, scaled_x, degree)
      val = beta[degree]
      deg = degree - 1
      while deg >= 0
        val = (val * scaled_x) + beta[deg]
        deg -= 1
      end
      val
    end
    private_class_method :eval_horner

    def summarize_lsm(cashflows, params, n_paths, n_steps, degree, antithetic, disc)
      samples = if antithetic
                  Array.new(n_paths / 2) { |j| 0.5 * disc * (cashflows[2 * j] + cashflows[(2 * j) + 1]) }
                else
                  cashflows.map { |cf| cf * disc }
                end

      stats = Stats.summarize_samples(samples)
      price_val = [stats[:price], [params.k - params.s0, 0.0].max].max
      LsmResult.new(
        n_paths: n_paths, n_steps: n_steps, degree: degree, antithetic: antithetic,
        price: price_val, se: stats[:se], ci_low: price_val - (Stats::CI_Z95 * stats[:se]),
        ci_high: price_val + (Stats::CI_Z95 * stats[:se]),
        early_exercise_premium: price_val - Stats.bs_put(params)
      )
    end
    private_class_method :summarize_lsm
  end
end
