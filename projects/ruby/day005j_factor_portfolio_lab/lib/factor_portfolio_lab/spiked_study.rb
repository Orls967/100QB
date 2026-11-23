# frozen_string_literal: true

module FactorPortfolioLab
  # Baik-Ben Arous-Peche (BBP) spiked covariance Monte Carlo verification study.
  module SpikedStudy
    module_function

    def run(master_seed:, n_assets: 40, t_len: 160, n_reps: 500, spikes: [0.0, 0.3, 0.5, 1.0, 2.0, 5.0])
      gamma = n_assets.to_f / t_len
      mp = Pca.marchenko_pastur_bounds(gamma, sigma2: 1.0)
      v_true = Array.new(n_assets, 0.0)
      v_true[0] = 1.0

      levels = spikes.each_with_index.map do |spike, idx|
        theory = Pca.bbp_prediction(spike, gamma, sigma2: 1.0)
        reps = Parallel.map_items(n_items: n_reps, chunk_size: 25) do |rep|
          seed = Rng.derive_seed(master_seed, 100 + idx, rep)
          run_single_rep(n_assets, t_len, spike, v_true, seed)
        end
        summarize_level(spike, theory, reps)
      end

      {
        n_assets: n_assets,
        t_len: t_len,
        gamma: gamma,
        n_reps: n_reps,
        marchenko_pastur: mp,
        levels: levels
      }
    end

    def run_single_rep(n_assets, t_len, spike, v_true, seed)
      rng = Rng::Sampler.new(seed)
      sd0 = Math.sqrt(1.0 + spike)
      rows = Array.new(t_len) do
        Array.new(n_assets) { |j| (j.zero? ? sd0 : 1.0) * rng.next_normal }
      end
      fit = Pca.fit(rows, center: false, ddof: 0)
      lam1 = fit[:eigenvalues][0]
      v1 = fit[:eigenvectors].map { |row| row[0] }
      align = Pca.eigenvector_alignment_sq(v1, v_true)
      {
        lambda1: lam1,
        lambda2: fit[:eigenvalues][1],
        lambda_min: fit[:eigenvalues][-1],
        alignment_sq: align,
        eff_rank: Pca.effective_rank(fit[:eigenvalues])
      }
    end

    def summarize_level(spike, theory, reps)
      n = reps.length.to_f
      l1 = reps.map { |r| r[:lambda1] }
      al = reps.map { |r| r[:alignment_sq] }
      l1_mean = l1.sum / n
      al_mean = al.sum / n
      l1_se = Math.sqrt(l1.sum { |x| (x - l1_mean)**2 } / [n - 1.0, 1.0].max / n)
      al_se = Math.sqrt(al.sum { |x| (x - al_mean)**2 } / [n - 1.0, 1.0].max / n)

      {
        spike: spike,
        above_threshold: theory[:above_threshold],
        theory_lambda1: theory[:sample_eigenvalue],
        emp_lambda1_mean: l1_mean,
        emp_lambda1_se: l1_se,
        theory_alignment_sq: theory[:alignment_sq],
        emp_alignment_sq_mean: al_mean,
        emp_alignment_sq_se: al_se,
        lambda2_mean: reps.sum { |r| r[:lambda2] } / n,
        lambda_min_mean: reps.sum { |r| r[:lambda_min] } / n,
        eff_rank_mean: reps.sum { |r| r[:eff_rank] } / n
      }
    end
  end
end
