# frozen_string_literal: true

module FactorPortfolioLab
  # Large-sample importance-sampling Monte Carlo verification of Black-Litterman Bayes parity.
  module BlStudy
    module_function

    def run(master_seed:, n_mc: 500_000, n_assets: 6, tau: 0.05)
      world = FactorWorld.build_world(n_assets: n_assets, n_factors: 2, seed: Rng.derive_seed(master_seed, 600, 0))
      sigma = world[:sigma]
      w_mkt = Array.new(n_assets, 1.0 / n_assets)
      pi = BlackLitterman.implied_equilibrium(sigma, w_mkt, risk_aversion: 2.5)

      p_mat = build_views_matrix(n_assets)
      q_vec = [pi[0] - pi[1] + 0.012, pi[2] + 0.008]
      omega = BlackLitterman.he_litterman_omega(sigma, p_mat, tau: tau)
      exact = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: tau)

      mc = run_importance_mc(master_seed, sigma, pi, p_mat, q_vec, omega, tau, n_mc)
      idzorek_sweep = run_idzorek_sweep(sigma, pi, p_mat, q_vec, tau)

      mu_err = exact[:mu_post].each_with_index.map { |m, i| (m - mc[:mu_mc][i]).abs }.max
      cov_err = max_mat_diff(exact[:sigma_mu_post], mc[:cov_mc])

      {
        n_mc: n_mc,
        n_assets: n_assets,
        tau: tau,
        pi: pi,
        mu_exact: exact[:mu_post],
        mu_mc: mc[:mu_mc],
        max_mu_diff: mu_err,
        max_cov_diff: cov_err,
        ess: mc[:ess],
        idzorek_sweep: idzorek_sweep
      }
    end

    def build_views_matrix(n_assets)
      row0 = Array.new(n_assets, 0.0)
      row0[0] = 1.0
      row0[1] = -1.0
      row1 = Array.new(n_assets, 0.0)
      row1[2] = 1.0
      [row0, row1]
    end

    def run_importance_mc(master_seed, sigma, pi, p_mat, q_vec, omega, tau, n_mc)
      n = sigma.length
      tau_sigma = BlackLitterman.scale_mat(sigma, tau)
      chol = LinAlg.cholesky(tau_sigma)
      omega_inv = LinAlg.invert(omega)
      n_chunks = 50
      per_chunk = n_mc / n_chunks

      chunks = Parallel.map_items(n_items: n_chunks, chunk_size: 5) do |c_idx|
        seed = Rng.derive_seed(master_seed, 601, c_idx)
        mc_chunk(per_chunk, n, pi, chol, p_mat, q_vec, omega_inv, seed)
      end

      aggregate_chunks(chunks, n)
    end

    def mc_chunk(count, n, pi, chol, p_mat, q_vec, omega_inv, seed)
      rng = Rng::Sampler.new(seed)
      log_ws = Array.new(count)
      mus = Array.new(count)
      max_lw = -1.0e30

      count.times do |r|
        z = Array.new(n) { rng.next_normal }
        delta = LinAlg.mat_vec(chol, z)
        mu = Array.new(n) { |i| pi[i] + delta[i] }
        p_mu = LinAlg.mat_vec(p_mat, mu)
        err = [p_mu[0] - q_vec[0], p_mu[1] - q_vec[1]]
        lw = -0.5 * LinAlg.quad_form(err, omega_inv)
        max_lw = lw if lw > max_lw
        log_ws[r] = lw
        mus[r] = mu
      end

      sum_w = 0.0
      sum_w2 = 0.0
      sum_mu = Array.new(n, 0.0)
      sum_m2 = Array.new(n) { Array.new(n, 0.0) }
      count.times do |r|
        w = Math.exp(log_ws[r] - max_lw)
        sum_w += w
        sum_w2 += w * w
        mu = mus[r]
        n.times do |i|
          wmi = w * mu[i]
          sum_mu[i] += wmi
          (i...n).each { |j| sum_m2[i][j] += wmi * mu[j] }
        end
      end

      { max_lw: max_lw, sum_w: sum_w, sum_w2: sum_w2, sum_mu: sum_mu, sum_m2: sum_m2 }
    end

    def aggregate_chunks(chunks, n)
      g_max = chunks.map { |c| c[:max_lw] }.max
      tot_w = 0.0
      tot_w2 = 0.0
      tot_mu = Array.new(n, 0.0)
      tot_m2 = Array.new(n) { Array.new(n, 0.0) }

      chunks.each do |c|
        scale = Math.exp(c[:max_lw] - g_max)
        tot_w += c[:sum_w] * scale
        tot_w2 += c[:sum_w2] * scale * scale
        n.times do |i|
          tot_mu[i] += c[:sum_mu][i] * scale
          (i...n).each { |j| tot_m2[i][j] += c[:sum_m2][i][j] * scale }
        end
      end

      mu_mc = tot_mu.map { |v| v / tot_w }
      cov_mc = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        (i...n).each do |j|
          val = (tot_m2[i][j] / tot_w) - (mu_mc[i] * mu_mc[j])
          cov_mc[i][j] = val
          cov_mc[j][i] = val
        end
      end

      { mu_mc: mu_mc, cov_mc: cov_mc, ess: (tot_w * tot_w) / tot_w2 }
    end

    def run_idzorek_sweep(sigma, pi, p_mat, q_vec, tau)
      [0.1, 0.25, 0.5, 0.75, 0.9].map do |conf|
        omega = BlackLitterman.idzorek_omega(sigma, p_mat, [conf, conf], tau: tau)
        post = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: tau)
        {
          confidence: conf,
          view0_implied: post[:mu_post][0] - post[:mu_post][1],
          view1_implied: post[:mu_post][2],
          shift_norm: Math.sqrt(post[:mu_post].each_with_index.sum { |m, i| (m - pi[i])**2 })
        }
      end
    end

    def max_mat_diff(a, b)
      m = 0.0
      a.each_with_index do |row, i|
        row.each_with_index do |v, j|
          d = (v - b[i][j]).abs
          m = d if d > m
        end
      end
      m
    end
  end
end
