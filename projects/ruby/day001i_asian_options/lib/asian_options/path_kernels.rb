# frozen_string_literal: true

module AsianOptions
  module PathKernels
    module_function

    def terminal_log_price(normals, spec)
      log_s = spec.log_s0
      drift = spec.drift
      vol = spec.vol
      n_steps = spec.n_steps
      step = 0
      while step < n_steps
        log_s += drift + (vol * normals[step])
        step += 1
      end
      log_s
    end

    def arithmetic_average(normals, spec)
      log_s = spec.log_s0
      drift = spec.drift
      vol = spec.vol
      n_steps = spec.n_steps
      sum_s = 0.0
      step = 0
      while step < n_steps
        log_s += drift + (vol * normals[step])
        sum_s += Math.exp(log_s)
        step += 1
      end
      sum_s * spec.inv_n
    end

    def antithetic_pair_payoff(normals, spec)
      log_pos = spec.log_s0
      log_neg = spec.log_s0
      drift = spec.drift
      vol = spec.vol
      n_steps = spec.n_steps
      sum_pos = 0.0
      sum_neg = 0.0
      step = 0
      while step < n_steps
        inc = vol * normals[step]
        log_pos += drift + inc
        log_neg += drift - inc
        sum_pos += Math.exp(log_pos)
        sum_neg += Math.exp(log_neg)
        step += 1
      end
      avg_pos = sum_pos * spec.inv_n
      avg_neg = sum_neg * spec.inv_n
      pay_pos = avg_pos > spec.k ? avg_pos - spec.k : 0.0
      pay_neg = avg_neg > spec.k ? avg_neg - spec.k : 0.0
      0.5 * (pay_pos + pay_neg) * spec.disc
    end

    def cv_path_payoffs(normals, spec)
      log_s = spec.log_s0
      drift = spec.drift
      vol = spec.vol
      n_steps = spec.n_steps
      sum_s = 0.0
      sum_log_s = 0.0
      step = 0
      while step < n_steps
        log_s += drift + (vol * normals[step])
        sum_s += Math.exp(log_s)
        sum_log_s += log_s
        step += 1
      end
      arith_avg = sum_s * spec.inv_n
      geom_avg = Math.exp(sum_log_s * spec.inv_n)
      x = arith_avg > spec.k ? (arith_avg - spec.k) * spec.disc : 0.0
      y = geom_avg > spec.k ? (geom_avg - spec.k) * spec.disc : 0.0
      [x, y]
    end
  end
end
