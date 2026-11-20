# frozen_string_literal: true

module Longevity
  module LeeCarter
    AGES = (30..100).to_a.freeze
    YEARS = (1980..2019).to_a.freeze
    TRUE_GM_A = 0.0006
    TRUE_GM_B = 0.00002
    TRUE_GM_C = 1.105
    TRUE_DRIFT = -1.0
    TRUE_SIGMA_K = 0.8
    TRUE_NOISE_SD = 0.03

    module_function

    def true_parameters(seed: 42)
      rng = Rng.new(seed)
      a_base = AGES.map { |x| Math.log(TRUE_GM_A + (TRUE_GM_B * (TRUE_GM_C**x))) }
      b_raw = AGES.map { |x| Math.exp(-(x - 30) / 50.0) }
      b_sum = b_raw.sum(0.0)
      b = b_raw.map { |v| v / b_sum }

      k = Array.new(YEARS.size, 0.0)
      (1...YEARS.size).each do |t|
        k[t] = k[t - 1] + TRUE_DRIFT + (TRUE_SIGMA_K * rng.next_normal)
      end
      k_mean = Stats.mean(k)
      k.map! { |v| v - k_mean }
      a = a_base.each_with_index.map { |ax, i| ax + (b[i] * k_mean) }

      { ages: AGES, years: YEARS, a: a, b: b, k: k, k_raw_mean: k_mean }
    end

    def simulate_log_mortality(truth, noise_sd: TRUE_NOISE_SD, seed: 42)
      rng = Rng.new(seed)
      a = truth[:a]
      b = truth[:b]
      k = truth[:k]
      n_ages = a.size
      n_years = k.size

      Array.new(n_ages) do |i|
        Array.new(n_years) do |t|
          eps = noise_sd.zero? ? 0.0 : noise_sd * rng.next_normal
          a[i] + (b[i] * k[t]) + eps
        end
      end
    end

    def fit(ln_m)
      n_ages = ln_m.size
      n_years = ln_m[0].size
      a_hat = ln_m.map { |row| Stats.mean(row) }
      centered = Array.new(n_ages) do |i|
        ai = a_hat[i]
        Array.new(n_years) { |t| ln_m[i][t] - ai }
      end

      svd_res = Linalg.svd(centered)
      u1 = svd_res[:u].map { |row| row[0] }
      v1 = svd_res[:v].map { |row| row[0] }
      s1 = svd_res[:s][0]

      sum_u = u1.sum(0.0)
      b_hat = u1.map { |val| val / sum_u }
      k_hat = v1.map { |val| s1 * sum_u * val }

      k_mean = Stats.mean(k_hat)
      k_hat.map! { |val| val - k_mean }
      a_hat.each_index { |i| a_hat[i] += b_hat[i] * k_mean }

      { a: a_hat, b: b_hat, k: k_hat }
    end
  end
end
