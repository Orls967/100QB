# frozen_string_literal: true

require_relative 'samplers'
require_relative 'stats'

module PopGen
  # Hardy-Weinberg equilibrium testing and partial selfing / inbreeding dynamics.
  module Inbreeding
    module_function

    # Sample n_diploid genotypes from HWE with allele frequency p and run chi-square test (df = 1)
    def hwe_sample_and_test(n_diploid, p, rng)
      q = 1.0 - p
      probs = [p * p, 2.0 * p * q, q * q]
      obs = Samplers.multinomial(n_diploid, probs, rng)
      n_aa, n_ab, _n_bb = obs
      p_hat = (2.0 * n_aa + n_ab) / (2.0 * n_diploid)
      q_hat = 1.0 - p_hat
      exp = [n_diploid * p_hat * p_hat, n_diploid * 2.0 * p_hat * q_hat, n_diploid * q_hat * q_hat]
      chi2 = 0.0
      3.times do |k|
        chi2 += (obs[k] - exp[k])**2 / exp[k] if exp[k] > 0.0
      end
      # 3 genotype classes minus 1 estimated parameter p_hat minus 1 sum constraint => df = 1
      { obs: obs, exp: exp, p_hat: p_hat, chi2: chi2, p_value: Stats.chi2_p_value(chi2, 1) }
    end

    # Theoretical inbreeding recursion under partial selfing rate s:
    # F_{t+1} = (s / 2) * (1 + F_t), converging to F* = s / (2 - s)
    def selfing_recursion(s, f0, gens)
      fs = [f0.to_f]
      f = f0.to_f
      gens.times do
        f = 0.5 * s * (1.0 + f)
        fs << f
      end
      fs
    end

    # Individual-based diploid population with partial selfing rate s and random outcrossing (1 - s).
    # Tracks empirical fixation index F_hat = 1 - H_obs / (2*p*(1-p)) over generations.
    def simulate_partial_selfing(n_diploid, p0, s_self, gens, rng)
      # Each individual: [allele_1, allele_2]
      pop = Array.new(n_diploid) do
        [rng.rand < p0 ? 1 : 0, rng.rand < p0 ? 1 : 0]
      end
      f_trace = []

      (0..gens).each do |g|
        n_het = 0
        n_a = 0
        pop.each do |a1, a2|
          n_a += a1 + a2
          n_het += 1 if a1 != a2
        end
        p_cur = n_a.to_f / (2.0 * n_diploid)
        h_exp = 2.0 * p_cur * (1.0 - p_cur)
        h_obs = n_het.to_f / n_diploid
        f_hat = h_exp > 1e-12 ? (1.0 - h_obs / h_exp) : 0.0
        f_trace << { gen: g, p: p_cur, h_obs: h_obs, h_exp: h_exp, f_hat: f_hat }
        break if g == gens

        next_pop = Array.new(n_diploid)
        n_diploid.times do |i|
          if rng.rand < s_self
            # Self-fertilization of one parent
            par = pop[rng.rand(n_diploid)]
            next_pop[i] = [par[rng.rand(2)], par[rng.rand(2)]]
          else
            # Outcrossing between two randomly drawn parents
            p1 = pop[rng.rand(n_diploid)]
            p2 = pop[rng.rand(n_diploid)]
            next_pop[i] = [p1[rng.rand(2)], p2[rng.rand(2)]]
          end
        end
        pop = next_pop
      end
      f_trace
    end
  end
end
