# frozen_string_literal: true

module FactorPortfolioLab
  # Michaud (1998) resampled efficient frontier vs plug-in Markowitz and Shrinkage.
  module Resampling
    module_function

    def run_michaud(master_seed:, n_assets: 15, t_len: 120, n_boot: 500, kappa: 0.5)
      world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 800, 0))
      sim = FactorWorld.simulate(world, t_len: t_len, seed: Rng.derive_seed(master_seed, 801, 0))
      rows = sim[:returns]

      mu_hat = Pca.column_means(rows, n_assets)
      sigma_hat = LinAlg.cov_matrix(rows)
      chol_hat = LinAlg.cholesky(sigma_hat)

      plugin = Markowitz.solve_qp(sigma_hat, mu: mu_hat, kappa: kappa)[:weights]
      lw = Shrinkage.shrink_cov_lw(rows)
      mu_js = Shrinkage.shrink_mu_js(mu_hat, lw[:sigma], t_len)
      shrunk = Markowitz.solve_qp(lw[:sigma], mu: mu_js, kappa: kappa)[:weights]

      boot_weights = Parallel.map_items(n_items: n_boot, chunk_size: 25) do |b_idx|
        seed = Rng.derive_seed(master_seed, 802, b_idx)
        single_bootstrap_weights(mu_hat, chol_hat, n_assets, t_len, kappa, seed)
      end

      resampled = average_weights(boot_weights, n_assets)
      true_opt = Markowitz.solve_qp(world[:sigma], mu: world[:mu], kappa: kappa)[:weights]

      {
        n_boot: n_boot,
        n_assets: n_assets,
        t_len: t_len,
        true_opt: evaluate_against_true(true_opt, world),
        plugin: evaluate_against_true(plugin, world).merge(l1_to_true: l1_dist(plugin, true_opt)),
        shrunk: evaluate_against_true(shrunk, world).merge(l1_to_true: l1_dist(shrunk, true_opt)),
        resampled: evaluate_against_true(resampled, world).merge(l1_to_true: l1_dist(resampled, true_opt)),
        boot_weight_dispersion: weight_dispersion(boot_weights, resampled, n_assets)
      }
    end

    def single_bootstrap_weights(mu_hat, chol_hat, n_assets, t_len, kappa, seed)
      rng = Rng::Sampler.new(seed)
      boot_rows = Array.new(t_len) do
        z = Array.new(n_assets) { rng.next_normal }
        delta = LinAlg.mat_vec(chol_hat, z)
        Array.new(n_assets) { |i| mu_hat[i] + delta[i] }
      end
      mu_b = Pca.column_means(boot_rows, n_assets)
      sig_b = LinAlg.cov_matrix(boot_rows)
      Markowitz.solve_qp(sig_b, mu: mu_b, kappa: kappa)[:weights]
    end

    def average_weights(list, n_assets)
      inv = 1.0 / list.length
      Array.new(n_assets) { |i| list.sum { |w| w[i] } * inv }
    end

    def evaluate_against_true(weights, world)
      er = LinAlg.dot(weights, world[:mu])
      vol = Math.sqrt([LinAlg.quad_form(weights, world[:sigma]), 0.0].max)
      {
        expected_return: er,
        volatility: vol,
        sharpe: vol <= 1e-12 ? 0.0 : er / vol,
        hhi: weights.sum { |w| w * w },
        n_active: weights.count { |w| w > 1.0e-4 },
        weights: weights
      }
    end

    def l1_dist(w1, w2)
      w1.each_with_index.sum { |v, i| (v - w2[i]).abs }
    end

    def weight_dispersion(boot_weights, mean_w, n_assets)
      n_assets.times.sum do |i|
        boot_weights.sum { |w| (w[i] - mean_w[i])**2 } / boot_weights.length
      end
    end
  end
end
