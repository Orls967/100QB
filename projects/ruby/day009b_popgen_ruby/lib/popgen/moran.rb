# frozen_string_literal: true

module PopGen
  # Moran birth-death process simulator with fitness r >= 0.
  module Moran
    module_function

    def simulate_to_absorption(n, init_count, r, rng, max_steps: 500_000)
      i = init_count
      steps = 0
      while i.positive? && i < n && steps < max_steps
        w_tot = r * i + (n - i)
        p_birth_a = (r * i) / w_tot
        birth_a = rng.rand < p_birth_a
        death_a = rng.rand < (i.to_f / n)
        if birth_a && !death_a
          i += 1
        elsif !birth_a && death_a
          i -= 1
        end
        steps += 1
      end
      { fixed: i == n, lost: i.zero?, steps: steps }
    end

    # Track heterozygosity H_t = 2*p*(1-p) over total_steps Moran steps
    def simulate_heterozygosity(n, init_count, total_steps, sample_every, rng)
      i = init_count
      samples = [2.0 * (i.to_f / n) * (1.0 - i.to_f / n)]
      (1..total_steps).each do |st|
        if i.positive? && i < n
          birth_a = rng.rand < (i.to_f / n)
          death_a = rng.rand < (i.to_f / n)
          if birth_a && !death_a
            i += 1
          elsif !birth_a && death_a
            i -= 1
          end
        end
        if (st % sample_every).zero?
          p = i.to_f / n
          samples << (2.0 * p * (1.0 - p))
        end
      end
      samples
    end
  end
end
