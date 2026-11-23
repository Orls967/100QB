# frozen_string_literal: true

module BacktestOverfitting
  module MultipleTesting
    module_function

    def unadjusted_rejections(pvals, alpha: 0.05)
      pvals.each_index.select { |i| pvals[i] <= alpha }
    end

    def bonferroni_rejections(pvals, alpha: 0.05)
      thresh = alpha / pvals.length
      pvals.each_index.select { |i| pvals[i] <= thresh }
    end

    def holm_rejections(pvals, alpha: 0.05)
      n = pvals.length
      order = (0...n).sort_by { |i| pvals[i] }
      rejected = []
      order.each_with_index do |idx, rank|
        break if pvals[idx] > alpha / (n - rank)

        rejected << idx
      end
      rejected
    end

    def benjamini_hochberg_rejections(pvals, alpha: 0.05)
      n = pvals.length
      order = (0...n).sort_by { |i| pvals[i] }
      k_max = -1
      order.each_with_index do |idx, rank|
        k_max = rank if pvals[idx] <= ((rank + 1) * alpha / n)
      end
      k_max >= 0 ? order[0..k_max] : []
    end

    def romano_wolf_rejections(centered_t_matrix, obs_abs_t, alpha: 0.05)
      n_hyp = obs_abs_t.length
      n_boot = centered_t_matrix.length
      remaining = (0...n_hyp).to_a
      rejected = []

      loop do
        break if remaining.empty?

        max_stats = centered_t_matrix.map { |row| remaining.map { |j| row[j] }.max }.sort!
        crit = max_stats[((1.0 - alpha) * n_boot).floor.clamp(0, n_boot - 1)]
        newly = remaining.select { |j| obs_abs_t[j] > crit }
        break if newly.empty?

        rejected.concat(newly)
        remaining -= newly
      end
      rejected
    end

    def run_study(seed: MAIN_SEED, workers: nil)
      indep = simulate_mt_regime(rho: 0.0, n_reps: 2_000, seed: seed, workers: workers)
      corr = simulate_mt_regime(rho: 0.4, n_reps: 2_000, seed: seed + 1, workers: workers)
      rw = romano_wolf_study(rho: 0.4, n_reps: 200, n_boot: 150, seed: seed, workers: workers)
      mkt = real_strategy_library_multiple_testing(seed: seed)
      { independent: indep, correlated_rho04: corr, romano_wolf_rho04: rw, ar1_market_library: mkt }
    end

    def simulate_mt_regime(rho:, n_reps:, seed:, workers:)
      n_tot = 300
      n_true = 30
      t_len = 252
      alpha = 0.05
      mean_sig = 0.0020
      sigma = 0.01

      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 100) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 700 + (rho * 100).to_i, cidx))
        acc = Array.new(4) { [0.0, 0.0, 0.0] }
        len.times do
          pvals = simulated_pvalues(n_tot, n_true, t_len, mean_sig, sigma, rho, rng)
          accumulate_mt_metrics(acc, pvals, n_true, alpha)
        end
        acc
      end

      summarize_mt(chunks, n_reps, n_tot, n_true, alpha)
    end

    def simulated_pvalues(n_tot, n_true, t_len, mean_sig, sigma, rho, rng)
      df = t_len - 1
      # Menggunakan korelasi faktor umum rho pada z_i untuk kecepatan dan keteraturan korelasi eksak
      inv_se = Math.sqrt(t_len) / sigma
      s_rho = Math.sqrt(rho)
      s_comp = Math.sqrt(1.0 - rho)
      f_bar = s_rho * rng.next_normal / Math.sqrt(t_len)

      Array.new(n_tot) do |i|
        mu = i < n_true ? mean_sig : 0.0
        eps_bar = s_comp * rng.next_normal * sigma / Math.sqrt(t_len)
        r_bar = mu + (sigma * f_bar) + eps_bar
        # Variansi t-stat dengan kuadrat akar chi2/(df)
        z_t = r_bar * inv_se
        chi_scale = Math.sqrt(1.0 + (rng.next_normal / Math.sqrt(2.0 * df)))
        t_stat = z_t / [chi_scale, 0.2].max
        Special.student_t_pvalue(t_stat, df)
      end
    end
    private_class_method :simulated_pvalues

    def accumulate_mt_metrics(acc, pvals, n_true, alpha)
      sets = [
        unadjusted_rejections(pvals, alpha: alpha),
        bonferroni_rejections(pvals, alpha: alpha),
        holm_rejections(pvals, alpha: alpha),
        benjamini_hochberg_rejections(pvals, alpha: alpha)
      ]
      sets.each_with_index do |rej, m_idx|
        next if rej.empty?

        false_pos = rej.count { |i| i >= n_true }
        true_pos = rej.length - false_pos
        acc[m_idx][0] += 1.0 if false_pos.positive?
        acc[m_idx][1] += false_pos.to_f / rej.length
        acc[m_idx][2] += true_pos.to_f / n_true
      end
    end
    private_class_method :accumulate_mt_metrics

    def summarize_mt(chunks, n_reps, n_tot, n_true, alpha)
      names = %w[unadjusted bonferroni holm bh]
      pi0 = (n_tot - n_true).to_f / n_tot
      out = { pi0: pi0, theoretical_bh_fdr: pi0 * alpha, alpha: alpha, n_reps: n_reps, methods: {} }
      names.each_with_index do |name, idx|
        fwer = chunks.sum(0.0) { |c| c[idx][0] } / n_reps
        fdr = chunks.sum(0.0) { |c| c[idx][1] } / n_reps
        pwr = chunks.sum(0.0) { |c| c[idx][2] } / n_reps
        out[:methods][name.to_sym] = { fwer: fwer, fdr: fdr, power: pwr }
      end
      out
    end
    private_class_method :summarize_mt

    def romano_wolf_study(rho:, n_reps:, n_boot:, seed:, workers:)
      n_tot = 60
      n_true = 6
      t_len = 252
      alpha = 0.05
      chunks = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 20) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 760, cidx))
        acc = [[0.0, 0.0], [0.0, 0.0]]
        len.times { accumulate_rw_rep(acc, n_tot, n_true, t_len, rho, n_boot, alpha, rng) }
        acc
      end
      { n_hyp: n_tot, n_true: n_true, n_reps: n_reps,
        bonferroni: { fwer: chunks.sum(0.0) { |c| c[0][0] } / n_reps,
                      power: chunks.sum(0.0) { |c| c[0][1] } / n_reps },
        romano_wolf: { fwer: chunks.sum(0.0) { |c| c[1][0] } / n_reps,
                       power: chunks.sum(0.0) { |c| c[1][1] } / n_reps } }
    end

    def accumulate_rw_rep(acc, n_tot, n_true, t_len, rho, n_boot, alpha, rng)
      mat, obs_abs_t, pvals = generate_panel_with_t(n_tot, n_true, t_len, rho, rng)
      centered = Array.new(n_boot) { bootstrap_max_row(mat, t_len, rng) }
      rw_rej = romano_wolf_rejections(centered, obs_abs_t, alpha: alpha)
      bonf_rej = bonferroni_rejections(pvals, alpha: alpha)
      [bonf_rej, rw_rej].each_with_index do |rej, idx|
        fp = rej.count { |i| i >= n_true }
        acc[idx][0] += 1.0 if fp.positive?
        acc[idx][1] += (rej.length - fp).to_f / n_true
      end
    end
    private_class_method :accumulate_rw_rep

    def generate_panel_with_t(n_tot, n_true, t_len, rho, rng)
      s_rho = Math.sqrt(rho)
      s_comp = Math.sqrt(1.0 - rho)
      f_ts = Array.new(t_len) { rng.next_normal }
      mat = Array.new(n_tot) do |i|
        mu = i < n_true ? 0.0022 : 0.0
        Array.new(t_len) { |t| mu + (0.01 * ((s_rho * f_ts[t]) + (s_comp * rng.next_normal))) }
      end
      obs_abs_t, pvals = t_stats_and_pvals(mat, t_len)
      [mat, obs_abs_t, pvals]
    end
    private_class_method :generate_panel_with_t

    def t_stats_and_pvals(mat, t_len)
      obs_abs = []
      pvals = []
      mat.each do |col|
        m = col.sum(0.0) / t_len
        sd = Math.sqrt(col.sum(0.0) { |x| (x - m)**2 } / (t_len - 1))
        t_val = m / (sd / Math.sqrt(t_len))
        obs_abs << t_val.abs
        pvals << Special.student_t_pvalue(t_val, t_len - 1)
      end
      [obs_abs, pvals]
    end
    private_class_method :t_stats_and_pvals

    def bootstrap_max_row(mat, t_len, rng)
      n_tot = mat.length
      idxs = Array.new(t_len) { rng.next_int(t_len) }
      col_means = mat.map { |col| col.sum(0.0) / t_len }
      Array.new(n_tot) do |j|
        col = mat[j]
        mu0 = col_means[j]
        s1 = 0.0
        s2 = 0.0
        idxs.each do |t|
          d = col[t] - mu0
          s1 += d
          s2 += d * d
        end
        m = s1 / t_len
        sd = Math.sqrt([((s2 - (t_len * m * m)) / (t_len - 1)), 1.0e-18].max)
        (m / (sd / Math.sqrt(t_len))).abs
      end
    end
    private_class_method :bootstrap_max_row

    def real_strategy_library_multiple_testing(seed:)
      rng = Rng::Sampler.new(seed)
      rets = DataGenerator.ar1_gaussian(600, phi: 0.12, sigma: 0.01, rng: rng)
      lib = Strategy.build_library
      pvals = lib.map do |strat|
        net = BacktestEngine.run(rets, strat.generate_positions(rets), cost_bps: 0.0)[:net_returns][200..]
        mom = Metrics.summary(net)
        t_val = mom[:sharpe_per_period] * Math.sqrt(net.length)
        Special.student_t_pvalue(t_val, net.length - 1)
      end
      { n_strategies: lib.length, unadjusted: unadjusted_rejections(pvals).length,
        bonferroni: bonferroni_rejections(pvals).length, holm: holm_rejections(pvals).length,
        bh: benjamini_hochberg_rejections(pvals).length }
    end
    private_class_method :real_strategy_library_multiple_testing
  end
end
