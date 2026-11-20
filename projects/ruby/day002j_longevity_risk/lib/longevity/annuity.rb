# frozen_string_literal: true

module Longevity
  module Annuity
    START_AGE = 65
    MAX_AGE = 100
    RATE = 0.03
    DISCOUNT = 1.0 / (1.0 + RATE)

    module_function

    # Extracts cohort q_{65+j, 2020+j} for j = 0..34 (age 65..99), with q_100 = 1.0.
    # k_future[j] corresponds to year 2020 + j (i.e. T + 1 + j).
    def cohort_qx(a_hat, b_hat, k_future, ages: LeeCarter::AGES)
      age_offset = START_AGE - ages.first
      n_steps = MAX_AGE - START_AGE
      Array.new(n_steps) do |j|
        idx = age_offset + j
        m = Math.exp(a_hat[idx] + (b_hat[idx] * k_future[j]))
        1.0 - Math.exp(-m)
      end
    end

    # Exact backward recursion and forward sum for deterministic life annuity-due/immediate (end of year).
    def deterministic_moments(qx_list, disc: DISCOUNT)
      px_list = qx_list.map { |q| 1.0 - q }
      n = px_list.size

      # Backward recursion: a_x = v * p_x * (1 + a_{x+1})
      ax = 0.0
      (n - 1).downto(0) do |j|
        ax = disc * px_list[j] * (1.0 + ax)
      end

      # Forward sum: sum_{k=1..n} v^k * _k p_x
      fwd = 0.0
      surv = 1.0
      vk = 1.0
      n.times do |j|
        surv *= px_list[j]
        vk *= disc
        fwd += vk * surv
      end

      # Exact E[Z] and E[Z^2] via backward recursion
      ez = 0.0
      ez2 = 0.0
      v2 = disc * disc
      (n - 1).downto(0) do |j|
        p = px_list[j]
        ez2_new = p * v2 * (1.0 + (2.0 * ez) + ez2)
        ez_new = p * disc * (1.0 + ez)
        ez = ez_new
        ez2 = ez2_new
      end

      {
        mean_backward: ax,
        mean_forward: fwd,
        variance: ez2 - (ez * ez)
      }
    end

    # Precomputes cumulative survival probabilities [_1 p_65, ..., _35 p_65] and annuity payoffs.
    def build_cdf_table(qx_list, disc: DISCOUNT)
      n = qx_list.size
      surv_probs = Array.new(n)
      annuity_payoffs = Array.new(n + 1, 0.0)
      surv = 1.0
      pv = 0.0
      vk = 1.0
      n.times do |j|
        surv *= (1.0 - qx_list[j])
        surv_probs[j] = surv
        vk *= disc
        pv += vk
        annuity_payoffs[j + 1] = pv
      end
      { surv_probs: surv_probs, payoffs: annuity_payoffs }
    end

    # Binary search on decreasing survival curve: returns curtate lifetime K in 0..35.
    # Surv_probs[k-1] = P(K >= k). Draw u ~ U(0,1); K = max { k : surv_probs[k-1] > u }.
    def sample_curtate_lifetime(surv_probs, uniform_u)
      lo = 0
      hi = surv_probs.size - 1
      return 0 if uniform_u >= surv_probs[0]

      while lo < hi
        mid = (lo + hi + 1) / 2
        if surv_probs[mid] > uniform_u
          lo = mid
        else
          hi = mid - 1
        end
      end
      lo + 1
    end

    def simulate_individual_lives(qx_list, num_lives:, base_seed:, workers: 1, disc: DISCOUNT, chunk_size: 2000)
      table = build_cdf_table(qx_list, disc: disc)
      surv_probs = table[:surv_probs]
      payoffs = table[:payoffs]

      n_chunks = (num_lives + chunk_size - 1) / chunk_size

      chunks = Parallel.map_indices(n_chunks, workers: workers, base_seed: base_seed) do |cid, cseed|
        rng = Rng.new(cseed)
        start_i = cid * chunk_size
        end_i = [num_lives, start_i + chunk_size].min
        vals = Array.new(end_i - start_i)
        vals.each_index do |idx|
          k_life = sample_curtate_lifetime(surv_probs, rng.next_uniform)
          vals[idx] = payoffs[k_life]
        end
        vals
      end

      all_vals = chunks.flatten(1)
      se = Stats.batch_means_se(all_vals, num_batches: 25)
      {
        mean: Stats.mean(all_vals),
        variance: Stats.variance(all_vals, ddof: 1),
        mean_se: se[:mean_se],
        var_se: se[:var_se]
      }
    end
  end
end
