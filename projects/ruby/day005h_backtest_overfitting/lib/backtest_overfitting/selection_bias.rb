# frozen_string_literal: true

module BacktestOverfitting
  module SelectionBias
    Z_95_ONE = 1.6448536269514722

    module_function

    def expected_max_exact(n_trials, n_gl: 160)
      nodes, weights = Special.gauss_legendre(n_gl, -8.5, 8.5)
      nodes.zip(weights).sum(0.0) do |x, w|
        cdf = Special.norm_cdf(x)
        pdf = Special.norm_pdf(x)
        w * x * n_trials * pdf * (cdf**(n_trials - 1))
      end
    end

    def expected_max_approx(n_trials)
      g = Special::EULER_GAMMA
      ((1.0 - g) * Special.norm_inv(1.0 - (1.0 / n_trials))) +
        (g * Special.norm_inv(1.0 - (1.0 / (n_trials * Math::E))))
    end

    def deflated_sharpe_ratio(sr_hat, sr0, n_obs, skew: 0.0, kurt: 3.0)
      denom = Math.sqrt([1.0 - (skew * sr_hat) + (0.25 * (kurt - 1.0) * sr_hat * sr_hat), 1.0e-12].max)
      z = (sr_hat - sr0) * Math.sqrt(n_obs - 1.0) / denom
      Special.norm_cdf(z)
    end

    def effective_trials_from_corr(n_trials, avg_corr)
      rho = avg_corr.clamp(0.0, 0.999)
      rho + ((1.0 - rho) * n_trials)
    end

    def average_pairwise_correlation(matrix_by_strat)
      n = matrix_by_strat.length
      t_len = matrix_by_strat.first.length
      zs = normalize_columns(matrix_by_strat, t_len)
      sum_c = Array.new(t_len, 0.0)
      zs.each do |col|
        t_len.times { |t| sum_c[t] += col[t] }
      end
      s2 = sum_c.sum(0.0) { |v| v * v } / t_len
      (s2 - n) / (n * (n - 1.0))
    end

    def normalize_columns(matrix_by_strat, t_len)
      matrix_by_strat.map do |col|
        m = col.sum(0.0) / t_len
        var = col.sum(0.0) { |x| (x - m)**2 } / t_len
        sd = var > 1.0e-24 ? Math.sqrt(var) : 1.0
        col.map { |x| (x - m) / sd }
      end
    end
    private_class_method :normalize_columns

    def run_study(seed: MAIN_SEED, workers: nil)
      table = [10, 50, 100, 300, 1_000].map do |n|
        ex = expected_max_exact(n)
        ap = expected_max_approx(n)
        { n_trials: n, exact: ex, approx: ap, abs_err: (ex - ap).abs }
      end

      lib_res = correlated_library_noise_study(n_reps: 120, n_obs: 252, seed: seed, workers: workers)
      dsr_cal = dsr_null_calibration_mc(n_trials: 300, n_obs: 252, n_reps: 2_000, seed: seed, workers: workers)
      { expected_max_table: table, library_noise_study: lib_res, dsr_calibration: dsr_cal }
    end

    def correlated_library_noise_study(n_reps:, n_obs:, seed:, workers:)
      lib = Strategy.build_library
      n_trials = lib.length
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 12) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 600, cidx))
        len.times.map do |k|
          rets = DataGenerator.iid_gaussian(n_obs, sigma: 0.01, rng: rng)
          evals = BacktestEngine.evaluate_library(lib, rets, cost_bps: 0.0)
          best_sr = evals.max_by { |e| e[:sharpe] }[:sharpe]
          avg_rho = cidx.zero? && k.zero? ? average_pairwise_correlation(evals.map { |e| e[:net_returns] }) : nil
          [best_sr, avg_rho]
        end
      end.flatten(1)

      bests = chunks.map { |x| x[0] }
      avg_rho = chunks.find { |x| !x[1].nil? }[1]
      summarize_library_noise(bests, avg_rho, n_trials, n_obs)
    end
    private_class_method :correlated_library_noise_study

    def summarize_library_noise(bests, avg_rho, n_trials, n_obs)
      m_best = bests.sum(0.0) / bests.length
      se_best = Math.sqrt((bests.sum(0.0) { |x| (x - m_best)**2 } / bests.length) / bests.length)
      n_eff = effective_trials_from_corr(n_trials, avg_rho)
      sigma_sr = 1.0 / Math.sqrt(n_obs)
      {
        n_trials: n_trials, n_obs: n_obs, n_reps: bests.length, avg_pairwise_corr: avg_rho,
        n_eff_corr_formula: n_eff, emp_mean_best_sharpe: m_best, emp_se_best_sharpe: se_best,
        theo_max_iid_300: sigma_sr * expected_max_exact(n_trials),
        theo_max_n_eff: sigma_sr * expected_max_exact(n_eff.round.clamp(2, n_trials)),
        best_sharpes_sample: bests.first(200)
      }
    end
    private_class_method :summarize_library_noise

    def dsr_null_calibration_mc(n_trials:, n_obs:, n_reps:, seed:, workers:)
      emax_n = expected_max_exact(n_trials)
      inv_sqrt_t = 1.0 / Math.sqrt(n_obs)
      rows = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 100) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 650, cidx))
        len.times.map { single_dsr_null_rep(n_trials, n_obs, emax_n, rng) }
      end.flatten(1)

      dsrs = rows.map { |r| r[0] }.sort
      naive_rate = rows.count { |r| r[1] }.to_f / n_reps
      dsr_rate = rows.count { |r| r[2] }.to_f / n_reps
      bests = rows.map { |r| r[3] }
      mean_best = bests.sum(0.0) / n_reps
      {
        n_trials: n_trials, n_obs: n_obs, n_reps: n_reps,
        emp_iid_mean_best_sharpe: mean_best, theo_iid_mean_best_sharpe: inv_sqrt_t * emax_n,
        naive_pass_rate_95: naive_rate, dsr_pass_rate_95: dsr_rate,
        dsr_median: dsrs[n_reps / 2], dsr_p90: dsrs[(0.90 * n_reps).to_i],
        dsr_p95: dsrs[(0.95 * n_reps).to_i], dsr_values_sample: dsrs.each_slice(10).map(&:first)
      }
    end
    private_class_method :dsr_null_calibration_mc

    def single_dsr_null_rep(n_trials, n_obs, emax_n, rng)
      best_sr = -1.0e9
      best_k = 3.0
      best_s = 0.0
      srs = Array.new(n_trials) do
        rets = DataGenerator.iid_gaussian(n_obs, sigma: 0.01, rng: rng)
        mom = Metrics.summary(rets)
        sr = mom[:sharpe_per_period]
        if sr > best_sr
          best_sr = sr
          best_s = mom[:skewness]
          best_k = mom[:kurtosis]
        end
        sr
      end
      m_sr = srs.sum(0.0) / n_trials
      sd_sr = Math.sqrt(srs.sum(0.0) { |x| (x - m_sr)**2 } / (n_trials - 1))
      dsr = deflated_sharpe_ratio(best_sr, sd_sr * emax_n, n_obs, skew: best_s, kurt: best_k)
      naive_z = best_sr / SharpeStats.se_sharpe(best_sr, n_obs)
      [dsr, naive_z > Z_95_ONE, dsr > 0.95, best_sr]
    end
    private_class_method :single_dsr_null_rep
  end
end
