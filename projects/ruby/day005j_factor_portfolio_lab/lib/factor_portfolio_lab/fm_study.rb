# frozen_string_literal: true

module FactorPortfolioLab
  # Monte Carlo verification of Fama-MacBeth vs Shanken (1992) EIV standard errors.
  module FmStudy
    module_function

    def run(master_seed:, n_assets: 40, n_factors: 3, t_len: 80, n_reps: 1000, mu_f_scale: 4.0)
      world = FactorWorld.build_world(
        n_assets: n_assets,
        n_factors: n_factors,
        seed: Rng.derive_seed(master_seed, 400, 0),
        alpha_scale: 0.0
      )
      scaled_mu_f = world[:mu_f].map { |m| m * mu_f_scale }
      scaled_world = world.merge(
        mu_f: scaled_mu_f,
        mu: FactorWorld.compute_expected_returns(world[:alpha], world[:betas], scaled_mu_f)
      )

      reps = Parallel.map_items(n_items: n_reps, chunk_size: 40) do |rep|
        seed = Rng.derive_seed(master_seed, 401, rep)
        sim = FactorWorld.simulate(scaled_world, t_len: t_len, seed: seed)
        fm = FamaMacbeth.fit(sim[:returns], sim[:factors], add_intercept: true)
        sh = Shanken.correct(fm, sim[:factors])
        {
          lambda_hat: fm[:lambda_hat],
          se_fm: fm[:se_fm],
          se_sh: sh[:se_shanken],
          c_sh: sh[:c_shanken],
          cs_r2: fm[:cs_r2]
        }
      end

      summarize(scaled_world, reps, n_reps)
    end

    def summarize(world, reps, n_reps)
      true_lambda = [0.0] + world[:mu_f]
      z_crit = 1.959963984540054
      stats = true_lambda.each_with_index.map { |tv, j| summarize_param(reps, tv, j, z_crit, n_reps) }

      {
        n_reps: n_reps,
        true_lambda: true_lambda,
        mean_lambda: stats.map { |s| s[:mean] },
        empirical_sd: stats.map { |s| s[:emp_sd] },
        mean_se_fm: stats.map { |s| s[:se_fm] },
        mean_se_shanken: stats.map { |s| s[:se_sh] },
        coverage_fm_95: stats.map { |s| s[:cov_fm] },
        coverage_shanken_95: stats.map { |s| s[:cov_sh] },
        mean_c_shanken: reps.sum { |r| r[:c_sh] } / n_reps,
        mean_cs_r2: reps.sum { |r| r[:cs_r2] } / n_reps
      }
    end

    def summarize_param(reps, true_val, idx, z_crit, n_reps)
      vals = reps.map { |r| r[:lambda_hat][idx] }
      mean = vals.sum / n_reps
      emp_sd = Math.sqrt(vals.sum { |x| (x - mean)**2 } / [n_reps - 1, 1].max)
      {
        mean: mean,
        emp_sd: emp_sd,
        se_fm: reps.sum { |r| r[:se_fm][idx] } / n_reps,
        se_sh: reps.sum { |r| r[:se_sh][idx] } / n_reps,
        cov_fm: reps.count { |r| (r[:lambda_hat][idx] - true_val).abs <= z_crit * r[:se_fm][idx] }.to_f / n_reps,
        cov_sh: reps.count { |r| (r[:lambda_hat][idx] - true_val).abs <= z_crit * r[:se_sh][idx] }.to_f / n_reps
      }
    end
  end
end
