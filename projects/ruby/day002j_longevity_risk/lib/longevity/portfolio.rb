# frozen_string_literal: true

module Longevity
  module Portfolio
    module_function

    # Computes conditional (mu_s, var_s) for each systematic k scenario s = 0...S-1.
    def scenario_moments(a_hat, b_hat, k_paths, workers: 1, chunk_size: 250)
      n_scen = k_paths.size
      n_chunks = (n_scen + chunk_size - 1) / chunk_size

      chunks = Parallel.map_indices(n_chunks, workers: workers, base_seed: 5005) do |cid, _cseed|
        start_i = cid * chunk_size
        end_i = [n_scen, start_i + chunk_size].min
        Array.new(end_i - start_i) do |offset|
          k_fut = k_paths[start_i + offset]
          qx = Annuity.cohort_qx(a_hat, b_hat, k_fut)
          mom = Annuity.deterministic_moments(qx)
          [mom[:mean_backward], mom[:variance]]
        end
      end

      flat = chunks.flatten(1)
      {
        mus: flat.map { |pair| pair[0] },
        vars: flat.map { |pair| pair[1] }
      }
    end

    def total_variance_summary(scen_moments, n_policies_list)
      mus = scen_moments[:mus]
      vars = scen_moments[:vars]
      mu_bar = Stats.mean(mus)
      e_var = Stats.mean(vars)
      var_mu = Stats.variance(mus, ddof: 1)

      rows = n_policies_list.map do |n_pol|
        if n_pol == Float::INFINITY
          cv = Math.sqrt(var_mu) / mu_bar
          { n: 'inf', mean: nil, var: nil, cv: cv }
        else
          mean_y = n_pol * mu_bar
          var_y = (n_pol * e_var) + ((n_pol * n_pol) * var_mu)
          cv = Math.sqrt((e_var / (n_pol * mu_bar * mu_bar)) + (var_mu / (mu_bar * mu_bar)))
          { n: n_pol, mean: mean_y, var: var_y, cv: cv }
        end
      end

      {
        mu_bar: mu_bar,
        e_sigma2: e_var,
        var_mu: var_mu,
        by_n: rows
      }
    end

    # Full individual-lives simulation across S scenarios, each with N policies.
    def simulate_portfolio_full(a_hat, b_hat, k_paths, num_policies:, base_seed:,
                                workers: 1, chunk_size: 100)
      n_scen = k_paths.size
      n_chunks = (n_scen + chunk_size - 1) / chunk_size

      chunks = Parallel.map_indices(n_chunks, workers: workers, base_seed: base_seed) do |cid, cseed|
        rng = Rng.new(cseed)
        start_i = cid * chunk_size
        end_i = [n_scen, start_i + chunk_size].min
        Array.new(end_i - start_i) do |offset|
          qx = Annuity.cohort_qx(a_hat, b_hat, k_paths[start_i + offset])
          table = Annuity.build_cdf_table(qx)
          surv_probs = table[:surv_probs]
          payoffs = table[:payoffs]
          total_y = 0.0
          num_policies.times do
            k_life = Annuity.sample_curtate_lifetime(surv_probs, rng.next_uniform)
            total_y += payoffs[k_life]
          end
          total_y
        end
      end

      ys = chunks.flatten(1)
      se = Stats.batch_means_se(ys, num_batches: 25)
      risk = Stats.var_and_es(ys)
      mean_y = Stats.mean(ys)
      var_y = Stats.variance(ys, ddof: 1)
      {
        n: num_policies,
        scenarios: n_scen,
        mean: mean_y,
        variance: var_y,
        cv: Math.sqrt(var_y) / mean_y,
        mean_se: se[:mean_se],
        var_se: se[:var_se],
        var995: risk[:var],
        es99: risk[:es],
        capital_pct: 100.0 * (risk[:var] - mean_y) / mean_y,
        samples: ys
      }
    end

    # Conditional normal approximation: Y | path ~ N(N * mu(path), N * sigma^2(path)).
    # Draws `draws_per_scen` normal samples per systematic scenario for smooth tail quantiles.
    def capital_conditional_normal(scen_moments, num_policies:, base_seed:, draws_per_scen: 10,
                                   workers: 1, chunk_size: 500)
      mus = scen_moments[:mus]
      vars = scen_moments[:vars]
      n_scen = mus.size
      n_chunks = (n_scen + chunk_size - 1) / chunk_size

      chunks = Parallel.map_indices(n_chunks, workers: workers, base_seed: base_seed) do |cid, cseed|
        rng = Rng.new(cseed)
        start_i = cid * chunk_size
        end_i = [n_scen, start_i + chunk_size].min
        out = []
        (start_i...end_i).each do |idx|
          m = num_policies * mus[idx]
          sd = Math.sqrt([0.0, num_policies * vars[idx]].max)
          draws_per_scen.times { out << (m + (sd * rng.next_normal)) }
        end
        out
      end

      ys = chunks.flatten(1)
      mean_y = Stats.mean(ys)
      risk = Stats.var_and_es(ys)
      {
        n: num_policies,
        mean: mean_y,
        var995: risk[:var],
        es99: risk[:es],
        capital_pct: 100.0 * (risk[:var] - mean_y) / mean_y
      }
    end

    # Deterministic longevity stress: shift k_{T+h} to its 0.5% quantile (z = -2.5758293035489004).
    def deterministic_stress_p005(a_hat, b_hat, rw_params, horizon: 35)
      z005 = -2.5758293035489004
      k_t = rw_params[:k_last]
      d_hat = rw_params[:d_hat]
      sigma_hat = rw_params[:sigma_hat]

      k_stress = Array.new(horizon) do |j|
        h = j + 1
        k_t + (h * d_hat) + (z005 * Math.sqrt(h.to_f) * sigma_hat)
      end
      qx = Annuity.cohort_qx(a_hat, b_hat, k_stress)
      Annuity.deterministic_moments(qx)[:mean_backward]
    end
  end
end
