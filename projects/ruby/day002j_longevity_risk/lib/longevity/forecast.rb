# frozen_string_literal: true

module Longevity
  module Forecast
    module_function

    def estimate_rw_drift(k_series)
      n_years = k_series.size
      d_hat = (k_series[-1] - k_series[0]) / (n_years - 1).to_f
      diffs = (1...n_years).map { |t| k_series[t] - k_series[t - 1] - d_hat }
      sigma_hat = Math.sqrt(diffs.sum(0.0) { |x| x * x } / (n_years - 2).to_f)
      {
        k_last: k_series[-1],
        d_hat: d_hat,
        sigma_hat: sigma_hat,
        t_years: n_years
      }
    end

    # Simulates N paths of length H=40 with optional parameter uncertainty on drift d.
    # Uses fixed chunks of 500 paths so results are strictly invariant to worker count.
    def simulate_k_paths(rw_params, horizon: 40, num_paths: 20_000, param_uncertainty: false,
                         base_seed: 1001, workers: 1, chunk_size: 500)
      k_t = rw_params[:k_last]
      d_hat = rw_params[:d_hat]
      sigma_hat = rw_params[:sigma_hat]
      d_se = sigma_hat / Math.sqrt((rw_params[:t_years] - 1).to_f)

      n_chunks = (num_paths + chunk_size - 1) / chunk_size
      chunks = Parallel.map_indices(n_chunks, workers: workers, base_seed: base_seed) do |cid, cseed|
        rng = Rng.new(cseed)
        start_i = cid * chunk_size
        end_i = [num_paths, start_i + chunk_size].min
        Array.new(end_i - start_i) do
          d_path = param_uncertainty ? d_hat + (d_se * rng.next_normal) : d_hat
          curr = k_t
          Array.new(horizon) do
            curr += d_path + (sigma_hat * rng.next_normal)
            curr
          end
        end
      end
      chunks.flatten(1)
    end

    def check_moments(paths, rw_params, horizons: [1, 10, 40], param_uncertainty: false)
      k_t = rw_params[:k_last]
      d_hat = rw_params[:d_hat]
      sigma_hat = rw_params[:sigma_hat]
      t_len = rw_params[:t_years]
      n_paths = paths.size

      horizons.map do |h|
        vals = paths.map { |p| p[h - 1] }
        emp_mean = Stats.mean(vals)
        emp_var = Stats.variance(vals, ddof: 1)
        theo_mean = k_t + (h * d_hat)
        theo_var = h * (sigma_hat**2)
        theo_var += (h * h) * (sigma_hat**2) / (t_len - 1).to_f if param_uncertainty
        se = Stats.batch_means_se(vals, num_batches: 25)
        mean_se = Math.sqrt(emp_var / n_paths.to_f)
        {
          h: h,
          emp_mean: emp_mean,
          theo_mean: theo_mean,
          mean_z: (emp_mean - theo_mean) / mean_se,
          emp_var: emp_var,
          theo_var: theo_var,
          var_se: se[:var_se],
          var_z: (emp_var - theo_var) / se[:var_se]
        }
      end
    end

    # Parametric bootstrap R=200 refits of Lee-Carter to measure empirical SD(d_hat).
    def parametric_bootstrap_refit(a_hat, b_hat, k_hat, rw_params, reps: 200,
                                   noise_sd: LeeCarter::TRUE_NOISE_SD, base_seed: 3003, workers: 1)
      n_ages = a_hat.size
      n_years = k_hat.size
      d_hat = rw_params[:d_hat]
      sigma_hat = rw_params[:sigma_hat]

      d_boots = Parallel.map_indices(reps, workers: workers, base_seed: base_seed) do |_r, rseed|
        rng = Rng.new(rseed)
        k_sim = Array.new(n_years, 0.0)
        (1...n_years).each { |t| k_sim[t] = k_sim[t - 1] + d_hat + (sigma_hat * rng.next_normal) }
        k_mean = Stats.mean(k_sim)
        k_sim.map! { |v| v - k_mean }

        ln_m_sim = Array.new(n_ages) do |i|
          ai = a_hat[i]
          bi = b_hat[i]
          Array.new(n_years) { |t| ai + (bi * k_sim[t]) + (noise_sd * rng.next_normal) }
        end

        fit_r = LeeCarter.fit(ln_m_sim)
        (fit_r[:k][-1] - fit_r[:k][0]) / (n_years - 1).to_f
      end

      {
        reps: reps,
        mean_d: Stats.mean(d_boots),
        sd_d: Stats.std_dev(d_boots, ddof: 1),
        theo_sd_d: sigma_hat / Math.sqrt((n_years - 1).to_f)
      }
    end

    # Cohort 2020 life expectancy trajectory at age 65 (curtate expected remaining years from age 65+j)
    # or remaining curtate life expectancy at age 65 across forecast paths for fan chart.
    def cohort_life_expectancy_fan(a_hat, b_hat, paths)
      age_offset = Annuity::START_AGE - LeeCarter::AGES.first
      n_steps = Annuity::MAX_AGE - Annuity::START_AGE

      # For each path, compute remaining curtate life expectancy e_{65} = sum_{k=1..35} _k p_65
      # and also survival curve S(j) = _j p_65 for j = 1..35
      surv_by_j = Array.new(n_steps) { [] }
      e65_vals = paths.map do |k_path|
        surv = 1.0
        e65 = 0.0
        n_steps.times do |j|
          idx = age_offset + j
          m = Math.exp(a_hat[idx] + (b_hat[idx] * k_path[j]))
          p = Math.exp(-m)
          surv *= p
          e65 += surv
          surv_by_j[j] << surv
        end
        e65
      end

      fan_survival = surv_by_j.each_with_index.map do |vals, j|
        s = vals.sort
        {
          age: Annuity::START_AGE + j + 1,
          year: 2020 + j + 1,
          q05: Stats.quantile_sorted(s, 0.05),
          q50: Stats.quantile_sorted(s, 0.50),
          q95: Stats.quantile_sorted(s, 0.95)
        }
      end

      e_sorted = e65_vals.sort
      {
        e65_q05: Stats.quantile_sorted(e_sorted, 0.05),
        e65_q50: Stats.quantile_sorted(e_sorted, 0.50),
        e65_q95: Stats.quantile_sorted(e_sorted, 0.95),
        e65_mean: Stats.mean(e65_vals),
        survival_fan: fan_survival
      }
    end
  end
end
