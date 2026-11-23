# frozen_string_literal: true

module BacktestOverfitting
  module SharpeStats
    Z_95 = 1.959963984540054

    module_function

    def se_sharpe(sr, n_obs, skew: 0.0, kurt: 3.0)
      var_term = (1.0 - (skew * sr) + (0.25 * (kurt - 1.0) * sr * sr)) / n_obs
      Math.sqrt([var_term, 1.0e-18].max)
    end

    def lo_autocorr_factor(q, phi)
      return 1.0 if q <= 1

      denom = q.to_f
      (1...(q)).each { |k| denom += 2.0 * (q - k) * (phi**k) }
      q.to_f / denom
    end

    def stationary_bootstrap_sharpe_ci(returns, mean_block_len:, n_boot:, rng:, alpha: 0.05)
      n = returns.length
      p_new = 1.0 / mean_block_len
      boots = Array.new(n_boot) do
        idx = rng.next_int(n)
        s1 = 0.0
        s2 = 0.0
        n.times do |k|
          idx = (rng.next_float < p_new ? rng.next_int(n) : (idx + 1) % n) if k.positive?
          r = returns[idx]
          s1 += r
          s2 += r * r
        end
        m = s1 / n
        sd = Math.sqrt([((s2 - (n * m * m)) / (n - 1)), 1.0e-24].max)
        m / sd
      end.sort!
      lo_i = [(0.5 * alpha * n_boot).floor, 0].max
      hi_i = [((1.0 - (0.5 * alpha)) * n_boot).floor, n_boot - 1].min
      [boots[lo_i], boots[hi_i]]
    end

    def run_mc_study(seed: MAIN_SEED, workers: nil)
      gauss_rows = [250, 1_000, 5_000].map do |t|
        gaussian_coverage_mc(n_obs: t, sr_true: 0.10, n_reps: 20_000, seed: seed, workers: workers)
      end
      nongauss_rows = [
        student_t_coverage_mc(n_obs: 1_000, df: 5, sr_true: 0.25, n_reps: 20_000, seed: seed, workers: workers),
        lognormal_coverage_mc(n_obs: 1_000, sigma_log: 0.45, sr_true: 0.25, n_reps: 20_000,
                              seed: seed, workers: workers)
      ]
      ar1_check = verify_lo_autocorr_mc(seed: seed, workers: workers)
      { gaussian_table: gauss_rows, non_gaussian_table: nongauss_rows, autocorr_lo: ar1_check }
    end

    def gaussian_coverage_mc(n_obs:, sr_true:, n_reps:, seed:, workers:)
      res = run_sharpe_mc_loop(n_reps: n_reps, stream_id: 200 + n_obs, seed: seed, workers: workers) do |rng|
        sigma = 0.01
        mu = sr_true * sigma
        rets = DataGenerator.iid_gaussian(n_obs, sigma: sigma, mu: mu, rng: rng)
        mom = Metrics.summary(rets)
        se_hat = se_sharpe(mom[:sharpe_per_period], n_obs, skew: 0.0, kurt: 3.0)
        [mom[:sharpe_per_period], se_hat, (mom[:sharpe_per_period] - sr_true).abs <= (Z_95 * se_hat) ? 1.0 : 0.0]
      end
      theo_se = se_sharpe(sr_true, n_obs, skew: 0.0, kurt: 3.0)
      { n_obs: n_obs, sr_true: sr_true, theo_se: theo_se, emp_sd: res[:emp_sd],
        mean_est_se: res[:mean_se], coverage_95: res[:coverage], nominal: 0.95 }
    end

    def student_t_coverage_mc(n_obs:, df:, sr_true:, n_reps:, seed:, workers:)
      kurt_true = 3.0 + (6.0 / (df - 4.0))
      res = run_nongaussian_mc(n_obs: n_obs, sr_true: sr_true, n_reps: n_reps,
                               stream_id: 300 + df, seed: seed, workers: workers) do |rng|
        sigma = 0.01
        DataGenerator.student_t_returns(n_obs, df: df, sigma: sigma, mu: sr_true * sigma, rng: rng)
      end
      theo_se = se_sharpe(sr_true, n_obs, skew: 0.0, kurt: kurt_true)
      res.merge(dist: "Student-t(#{df})", n_obs: n_obs, sr_true: sr_true, theo_se: theo_se,
                naive_gauss_se: se_sharpe(sr_true, n_obs))
    end

    def lognormal_coverage_mc(n_obs:, sigma_log:, sr_true:, n_reps:, seed:, workers:)
      w = Math.exp(sigma_log * sigma_log)
      sd_raw = Math.sqrt(w * (w - 1.0))
      mu_log = Math.log((sr_true * sd_raw) + 1.0) - (0.5 * sigma_log * sigma_log)
      exact_sr = (Math.exp(mu_log + (0.5 * sigma_log * sigma_log)) - 1.0) /
                 (Math.exp(mu_log) * sd_raw)
      skew = (w + 2.0) * Math.sqrt(w - 1.0)
      kurt = (w**4) + (2.0 * (w**3)) + (3.0 * w * w) - 3.0
      res = run_nongaussian_mc(n_obs: n_obs, sr_true: exact_sr, n_reps: n_reps,
                               stream_id: 410, seed: seed, workers: workers) do |rng|
        DataGenerator.lognormal_returns(n_obs, mu_log: mu_log, sigma_log: sigma_log, rng: rng)
      end
      res.merge(dist: "Lognormal(#{sigma_log})", n_obs: n_obs, sr_true: exact_sr,
                theo_se: se_sharpe(exact_sr, n_obs, skew: skew, kurt: kurt),
                naive_gauss_se: se_sharpe(exact_sr, n_obs))
    end

    def verify_lo_autocorr_mc(seed:, workers:)
      phi = 0.3
      n_obs = 500
      sr_true = 0.05
      mean_blk = 10.0
      agg_q = verify_eta_aggregation(phi: phi, q_list: [5, 21], n_reps: 2_000, seed: seed, workers: workers)
      t_sweep = eta_finite_sample_sweep(phi: phi, q_list: [5, 21], t_list: [250, 1_000, 5_000, 20_000],
                                        seed: seed, workers: workers)
      cov = ar1_bootstrap_coverage_mc(phi: phi, n_obs: n_obs, sr_true: sr_true,
                                      mean_blk: mean_blk, n_reps: 800, seed: seed, workers: workers)
      { phi: phi, n_obs: n_obs, mean_block_len: mean_blk, q_scaling: agg_q,
        finite_sample_sweep: t_sweep }.merge(cov)
    end

    def eta_finite_sample_sweep(phi:, q_list:, t_list:, seed:, workers:, n_reps: 400)
      t_list.map do |t_len|
        chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 40) do |cidx, len|
          rng = Rng::Sampler.new(Rng.derive_seed(seed, 530 + t_len, cidx))
          run_finite_t_chunk(len, t_len, phi, q_list, rng)
        end
        summarize_finite_t_row(t_len, phi, q_list, chunks, n_reps)
      end
    end

    def run_finite_t_chunk(len, t_len, phi, q_list, rng)
      acc = { v1: 0.0 }
      q_list.each { |q| acc[q] = 0.0 }
      len.times do
        rets = DataGenerator.ar1_gaussian(t_len, phi: phi, sigma: 0.01, mu: 0.0, rng: rng)
        acc[:v1] += sample_var(rets)
        q_list.each do |q|
          blocks = (0...(t_len / q)).map { |b| rets[b * q, q].sum(0.0) }
          acc[q] += sample_var(blocks)
        end
      end
      acc
    end
    private_class_method :run_finite_t_chunk

    def sample_var(arr)
      n = arr.length
      m = arr.sum(0.0) / n
      arr.sum(0.0) { |x| (x - m)**2 } / (n - 1)
    end
    private_class_method :sample_var

    def summarize_finite_t_row(t_len, phi, q_list, chunks, n_reps)
      mv1 = chunks.sum(0.0) { |c| c[:v1] } / n_reps
      qs = {}
      q_list.each do |q|
        mvq = chunks.sum(0.0) { |c| c[q] } / n_reps
        n_blk = t_len / q
        eta_asymp = lo_autocorr_factor(q, phi)
        num = (t_len.to_f / (t_len - 1)) * (1.0 - (1.0 / (t_len * lo_autocorr_factor(t_len, phi))))
        den = (n_blk.to_f / (n_blk - 1)) * (1.0 - (eta_asymp / (n_blk * lo_autocorr_factor(n_blk * q, phi))))
        qs[q] = { eta_asymp: eta_asymp, eta_finite_theory: eta_asymp * (num / den), eta_mc_ratio: q * mv1 / mvq }
      end
      { n_obs: t_len, by_q: qs }
    end
    private_class_method :summarize_finite_t_row

    def verify_eta_aggregation(phi:, q_list:, n_reps:, seed:, workers:)
      sigma = 0.01
      sr1_true = 0.05
      mu = sr1_true * sigma
      all_q = [1] + q_list
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 200) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 510, cidx))
        sums = Hash.new(0.0)
        len.times do
          rets = DataGenerator.ar1_gaussian(420, phi: phi, sigma: sigma, mu: mu, rng: rng)
          all_q.each { |q| accumulate_q_blocks(rets, q, sums) }
        end
        sums
      end
      sr1_emp = sr_from_bucket(1, chunks)
      q_list.map { |q| summarize_q_row(q, phi, sr1_true, sr1_emp, chunks) }
    end
    private_class_method :verify_eta_aggregation

    def accumulate_q_blocks(rets, q, sums)
      (0...(rets.length / q)).each do |b|
        s = rets[(b * q)...((b + 1) * q)].sum(0.0)
        sums[q] += s
        sums[q + 1000] += s * s
        sums[q + 2000] += 1.0
      end
    end
    private_class_method :accumulate_q_blocks

    def sr_from_bucket(q, chunks)
      n = chunks.sum(0.0) { |h| h[q + 2000] }
      m = chunks.sum(0.0) { |h| h[q] } / n
      var = (chunks.sum(0.0) { |h| h[q + 1000] } / n) - (m * m)
      m / Math.sqrt(var)
    end
    private_class_method :sr_from_bucket

    def summarize_q_row(q, phi, sr1_true, sr1_emp, chunks)
      sr_q_mc = sr_from_bucket(q, chunks)
      eta_theor = lo_autocorr_factor(q, phi)
      eta_proper = (sr_q_mc / (Math.sqrt(q) * sr1_emp))**2
      eta_old = (sr_q_mc / (Math.sqrt(q) * sr1_true))**2
      { q: q, eta_theoretical: eta_theor, eta_mc: eta_proper, eta_mc_pooled_true_sr1: eta_old,
        sr1_emp: sr1_emp, sr1_true: sr1_true, sr_q_mc: sr_q_mc,
        sr_q_theor: sr1_true * Math.sqrt(q * eta_theor), sr_q_naive: sr1_true * Math.sqrt(q) }
    end
    private_class_method :summarize_q_row

    def ar1_bootstrap_coverage_mc(phi:, n_obs:, sr_true:, mean_blk:, n_reps:, seed:, workers:)
      sigma = 0.01
      mu = sr_true * sigma
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 40) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 550, cidx))
        cov_naive = 0
        cov_boot = 0
        len.times do
          rets = DataGenerator.ar1_gaussian(n_obs, phi: phi, sigma: sigma, mu: mu, rng: rng)
          sr = Metrics.sharpe_per_period(rets)
          se_naive = se_sharpe(sr, n_obs)
          cov_naive += 1 if (sr - sr_true).abs <= (Z_95 * se_naive)
          lo, hi = stationary_bootstrap_sharpe_ci(rets, mean_block_len: mean_blk, n_boot: 250, rng: rng)
          cov_boot += 1 if sr_true.between?(lo, hi)
        end
        [cov_naive, cov_boot]
      end
      p_naive = chunks.sum(0.0) { |c| c[0] } / n_reps
      p_boot = chunks.sum(0.0) { |c| c[1] } / n_reps
      se_boot = Math.sqrt(p_boot * (1.0 - p_boot) / n_reps)
      { n_reps: n_reps, naive_coverage_95: p_naive,
        naive_coverage_se: Math.sqrt(p_naive * (1.0 - p_naive) / n_reps),
        stationary_bootstrap_coverage_95: p_boot, boot_coverage_se: se_boot,
        margin_above_091_se: (p_boot - 0.91) / se_boot,
        effective_bootstrap_phi: (1.0 - (1.0 / mean_blk)) * phi, nominal: 0.95 }
    end
    private_class_method :ar1_bootstrap_coverage_mc

    def run_nongaussian_mc(n_obs:, sr_true:, n_reps:, stream_id:, seed:, workers:, &gen)
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 500) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, cidx))
        s1 = 0.0
        s2 = 0.0
        cov_ng = 0
        cov_na = 0
        len.times do
          mom = Metrics.summary(gen.call(rng))
          sr = mom[:sharpe_per_period]
          s1 += sr
          s2 += sr * sr
          cov_ng += 1 if (sr - sr_true).abs <= (Z_95 * se_sharpe(sr, n_obs, skew: mom[:skewness], kurt: mom[:kurtosis]))
          cov_na += 1 if (sr - sr_true).abs <= (Z_95 * se_sharpe(sr, n_obs))
        end
        [s1, s2, cov_ng, cov_na]
      end
      m = chunks.sum(0.0) { |c| c[0] } / n_reps
      sd = Math.sqrt([((chunks.sum(0.0) { |c| c[1] } / n_reps) - (m * m)), 0.0].max)
      { emp_sd: sd, coverage_95_moment_adj: chunks.sum(0.0) { |c| c[2] } / n_reps,
        coverage_95_naive_gauss: chunks.sum(0.0) { |c| c[3] } / n_reps }
    end
    private_class_method :run_nongaussian_mc

    def run_sharpe_mc_loop(n_reps:, stream_id:, seed:, workers:)
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 500) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, cidx))
        s1 = 0.0
        s2 = 0.0
        sse = 0.0
        cov = 0.0
        len.times do
          sr, se, hit = yield(rng)
          s1 += sr
          s2 += sr * sr
          sse += se
          cov += hit
        end
        [s1, s2, sse, cov]
      end
      m = chunks.sum(0.0) { |c| c[0] } / n_reps
      sd = Math.sqrt([((chunks.sum(0.0) { |c| c[1] } / n_reps) - (m * m)), 0.0].max)
      { emp_sd: sd, mean_se: chunks.sum(0.0) { |c| c[2] } / n_reps, coverage: chunks.sum(0.0) { |c| c[3] } / n_reps }
    end
    private_class_method :run_sharpe_mc_loop
  end
end
