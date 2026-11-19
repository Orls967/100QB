# frozen_string_literal: true

module AmericanPut
  module Crr
    module_function

    def price(params, n_steps: DEFAULT_CRR_STEPS)
      params.validate!
      raise ArgumentError, 'n_steps harus >= 2' if n_steps < 2

      tree = build_tree_constants(params, n_steps)
      am_vals = terminal_put_payoffs(params.s0, params.k, tree[:ln_u], n_steps)
      eu_vals = am_vals.dup

      greeks = backward_induct!(am_vals, eu_vals, params, n_steps, tree)
      am_price = am_vals[0]
      eu_price = eu_vals[0]

      CrrResult.new(
        n_steps: n_steps, american_price: am_price, european_price: eu_price,
        early_exercise_premium: am_price - eu_price,
        delta: greeks[:delta], gamma: greeks[:gamma], theta: greeks[:theta]
      )
    end

    def convergence_sweep(params, steps_grid: CRR_STEPS_GRID, workers: nil)
      Parallel.map_tasks(steps_grid, workers: workers) do |n_steps, _idx|
        price(params, n_steps: n_steps)
      end
    end

    def build_tree_constants(params, n_steps)
      step_dt = params.t / n_steps.to_f
      ln_u = params.sigma * Math.sqrt(step_dt)
      up = Math.exp(ln_u)
      down = 1.0 / up
      disc = Math.exp(-params.r * step_dt)
      prob_up = (Math.exp(params.r * step_dt) - down) / (up - down)
      {
        dt: step_dt, ln_u: ln_u, u_sq: Math.exp(2.0 * ln_u),
        disc_p: disc * prob_up, disc_q: disc * (1.0 - prob_up)
      }
    end
    private_class_method :build_tree_constants

    def terminal_put_payoffs(spot0, strike, ln_u, n_steps)
      u_sq = Math.exp(2.0 * ln_u)
      spot = spot0 * Math.exp(-n_steps * ln_u)
      Array.new(n_steps + 1) do
        payoff = strike > spot ? strike - spot : 0.0
        spot *= u_sq
        payoff
      end
    end
    private_class_method :terminal_put_payoffs

    def backward_induct!(am_vals, eu_vals, params, n_steps, tree)
      saved = {}
      step = n_steps - 1
      while step >= 0
        spot = params.s0 * Math.exp(-step * tree[:ln_u])
        step_backward!(am_vals, eu_vals, step, spot, params.k, tree)
        capture_greek_slice(saved, am_vals, params.s0, tree[:ln_u], step) if step <= 2
        step -= 1
      end

      compute_greeks(saved, tree[:dt])
    end
    private_class_method :backward_induct!

    def step_backward!(am_vals, eu_vals, step, spot0, strike, tree)
      spot = spot0
      disc_p = tree[:disc_p]
      disc_q = tree[:disc_q]
      u_sq = tree[:u_sq]
      j = 0
      while j <= step
        eu_vals[j] = (disc_p * eu_vals[j + 1]) + (disc_q * eu_vals[j])
        cont = (disc_p * am_vals[j + 1]) + (disc_q * am_vals[j])
        intr = strike > spot ? strike - spot : 0.0
        am_vals[j] = [cont, intr].max
        spot *= u_sq
        j += 1
      end
    end
    private_class_method :step_backward!

    def capture_greek_slice(saved, am_vals, spot0, ln_u, step)
      u_sq = Math.exp(2.0 * ln_u)
      base_spot = spot0 * Math.exp(-step * ln_u)
      saved[step] = {
        vals: am_vals[0..step],
        spots: Array.new(step + 1) { |j| base_spot * (u_sq**j) }
      }
    end
    private_class_method :capture_greek_slice

    def compute_greeks(saved, step_dt)
      v1, s1 = saved[1].values_at(:vals, :spots)
      v2, s2 = saved[2].values_at(:vals, :spots)
      {
        delta: (v1[1] - v1[0]) / (s1[1] - s1[0]),
        gamma: compute_gamma(v2, s2),
        theta: (v2[1] - saved[0][:vals][0]) / (2.0 * step_dt)
      }
    end
    private_class_method :compute_greeks

    def compute_gamma(vals, spots)
      d_up = (vals[2] - vals[1]) / (spots[2] - spots[1])
      d_dn = (vals[1] - vals[0]) / (spots[1] - spots[0])
      (d_up - d_dn) / (0.5 * (spots[2] - spots[0]))
    end
    private_class_method :compute_gamma
  end
end
