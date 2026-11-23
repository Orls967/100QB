# frozen_string_literal: true

module FactorPortfolioLab
  # Portfolio constraints: capped-simplex Euclidean projection and KKT optimality residuals.
  module Constraints
    module_function

    def project_capped_simplex(vec, lower:, upper:, target_sum: 1.0, tol: 1.0e-13)
      n = vec.length
      sum_l = lower.sum
      sum_u = upper.sum
      if sum_l > target_sum + 1.0e-9 || sum_u < target_sum - 1.0e-9
        raise ArgumentError,
              'infeasible bounds for target_sum'
      end

      lo = (0...n).map { |i| vec[i] - upper[i] }.min - 1.0
      hi = (0...n).map { |i| vec[i] - lower[i] }.max + 1.0

      80.times do
        mid = 0.5 * (lo + hi)
        s = 0.0
        n.times { |i| s += clamp(vec[i] - mid, lower[i], upper[i]) }
        break if (s - target_sum).abs <= tol

        if s > target_sum
          lo = mid
        else
          hi = mid
        end
      end

      nu = 0.5 * (lo + hi)
      Array.new(n) { |i| clamp(vec[i] - nu, lower[i], upper[i]) }
    end

    def clamp(val, lo, hi)
      return lo if val < lo
      return hi if val > hi

      val
    end

    def kkt_residual(weights, grad, lower:, upper:, tol_bound: 1.0e-10)
      n = weights.length
      free_grads = []
      n.times do |i|
        free_grads << grad[i] if weights[i] > lower[i] + tol_bound && weights[i] < upper[i] - tol_bound
      end

      nu = if free_grads.empty?
             estimate_multiplier(weights, grad, lower, upper, tol_bound)
           else
             free_grads.sum / free_grads.length
           end

      max_stat = 0.0
      n.times do |i|
        res = if weights[i] <= lower[i] + tol_bound
                [nu - grad[i], 0.0].max
              elsif weights[i] >= upper[i] - tol_bound
                [grad[i] - nu, 0.0].max
              else
                (grad[i] - nu).abs
              end
        max_stat = res if res > max_stat
      end

      sum_err = (weights.sum - 1.0).abs
      bound_err = n.times.map do |i|
        [lower[i] - weights[i], weights[i] - upper[i], 0.0].max
      end.max

      {
        stationarity: max_stat,
        primal_sum_err: sum_err,
        primal_bound_err: bound_err,
        max_kkt: [max_stat, sum_err, bound_err].max,
        nu: nu
      }
    end

    def estimate_multiplier(weights, grad, lower, upper, tol_bound)
      lo = -1.0e12
      hi = 1.0e12
      weights.each_with_index do |w, i|
        lo = [lo, grad[i]].max if w >= upper[i] - tol_bound
        hi = [hi, grad[i]].min if w <= lower[i] + tol_bound
      end
      return 0.5 * (lo + hi) if lo <= hi && lo > -1.0e11 && hi < 1.0e11

      grad.sum / grad.length
    end
  end
end
