# frozen_string_literal: true

module Longevity
  module GompertzMakeham
    module_function

    def log_mu(age, param_a, param_b, param_c)
      Math.log(param_a + (param_b * (param_c**age)))
    end

    def fit(ages, a_hat, initial: [0.0005, 0.00003, 1.10])
      u0 = [
        Math.log(initial[0]),
        Math.log(initial[1]),
        Math.log(initial[2] - 1.0)
      ]

      res = Optimizers.nelder_mead(u0, step: [0.25, 0.25, 0.15], tol: 1e-13, max_iter: 4000) do |u|
        pa = Math.exp(u[0])
        pb = Math.exp(u[1])
        pc = 1.0 + Math.exp(u[2])
        sse = 0.0
        ages.each_with_index do |x, i|
          diff = Math.log(pa + (pb * (pc**x))) - a_hat[i]
          sse += diff * diff
        end
        sse
      end

      u_opt = res[:x]
      {
        a: Math.exp(u_opt[0]),
        b: Math.exp(u_opt[1]),
        c: 1.0 + Math.exp(u_opt[2]),
        sse: res[:f]
      }
    end
  end
end
