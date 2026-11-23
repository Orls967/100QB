# frozen_string_literal: true

module FactorPortfolioLab
  # Large-scale Monte Carlo study comparing exact GRS F test vs asymptotic Wald chi-squared test.
  module GrsStudy
    module_function

    def run(master_seed:, n_reps: 5000, regimes: nil)
      regimes ||= [
        { label: 'N10_T120_H0', n_assets: 10, t_len: 120, n_factors: 3, alpha_scale: 0.0 },
        { label: 'N25_T120_H0', n_assets: 25, t_len: 120, n_factors: 3, alpha_scale: 0.0 },
        { label: 'N50_T120_H0', n_assets: 50, t_len: 120, n_factors: 3, alpha_scale: 0.0 },
        { label: 'N50_T600_H0', n_assets: 50, t_len: 600, n_factors: 3, alpha_scale: 0.0 },
        { label: 'N25_T120_H1', n_assets: 25, t_len: 120, n_factors: 3, alpha_scale: 0.004 }
      ]

      results = regimes.each_with_index.map do |reg, idx|
        world = FactorWorld.build_world(
          n_assets: reg[:n_assets],
          n_factors: reg[:n_factors],
          seed: Rng.derive_seed(master_seed, 200, idx),
          alpha_scale: reg[:alpha_scale]
        )
        run_regime(reg, world, master_seed, idx, n_reps)
      end

      { n_reps: n_reps, regimes: results }
    end

    def run_regime(reg, world, master_seed, idx, n_reps)
      reps = Parallel.map_items(n_items: n_reps, chunk_size: 100) do |rep|
        seed = Rng.derive_seed(master_seed, 300 + idx, rep)
        sim = FactorWorld.simulate(world, t_len: reg[:t_len], seed: seed)
        res = Grs.test(sim[:returns], sim[:factors])
        [res[:grs_f], res[:p_value_f], res[:wald_chi2], res[:p_value_wald]]
      end

      f_rej = reps.count { |r| r[1] < 0.05 }.to_f / n_reps
      w_rej = reps.count { |r| r[3] < 0.05 }.to_f / n_reps
      {
        label: reg[:label],
        n_assets: reg[:n_assets],
        t_len: reg[:t_len],
        n_factors: reg[:n_factors],
        alpha_scale: reg[:alpha_scale],
        grs_rejection_05: f_rej,
        grs_rejection_se: Math.sqrt(f_rej * (1.0 - f_rej) / n_reps),
        wald_rejection_05: w_rej,
        wald_rejection_se: Math.sqrt(w_rej * (1.0 - w_rej) / n_reps),
        mean_grs_f: reps.sum { |r| r[0] } / n_reps,
        mean_wald_chi2: reps.sum { |r| r[2] } / n_reps
      }
    end
  end
end
