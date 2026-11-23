# frozen_string_literal: true

module GarchFamily
  module MisspecPortfolio
    BREAK_BASE = { mu: 0.0, omega1: 3.0e-6, omega2: 1.8e-5, alpha: 0.08, beta: 0.77 }.freeze
    PORT_WEIGHTS = [0.25, 0.20, 0.20, 0.20, 0.15].freeze
    PORT_GJR = [
      { mu: 0.0002, omega: 3.0e-6, alpha: 0.03, gamma: 0.07, beta: 0.87, nu: 6.0 },
      { mu: 0.0003, omega: 5.0e-6, alpha: 0.04, gamma: 0.09, beta: 0.84, nu: 5.5 },
      { mu: 0.0004, omega: 7.0e-6, alpha: 0.05, gamma: 0.10, beta: 0.82, nu: 5.0 },
      { mu: 0.0001, omega: 1.5e-6, alpha: 0.04, gamma: 0.04, beta: 0.88, nu: 7.0 },
      { mu: 0.0001, omega: 4.0e-6, alpha: 0.05, gamma: 0.05, beta: 0.86, nu: 6.0 }
    ].freeze

    module_function

    def simulate_break_series(n_obs: 3000, seed: MAIN_SEED)
      rng = Rng::Sampler.new(seed)
      break_idx = n_obs / 2
      h = BREAK_BASE[:omega1] / (1.0 - BREAK_BASE[:alpha] - BREAK_BASE[:beta])
      rets = Array.new(n_obs)
      eps = 0.0
      (n_obs + 400).times do |t|
        w = t >= (break_idx + 400) ? BREAK_BASE[:omega2] : BREAK_BASE[:omega1]
        h = w + (BREAK_BASE[:alpha] * eps * eps) + (BREAK_BASE[:beta] * h) if t.positive?
        eps = Math.sqrt(h) * rng.next_normal
        rets[t - 400] = eps if t >= 400
      end
      { returns: rets, break_idx: break_idx }
    end

    def inclan_tiao_cusum(rets)
      n = rets.length
      c_all = rets.sum(0.0) { |r| r * r }
      cum = 0.0
      best_k = 0
      best_val = 0.0
      (n - 1).times do |k|
        cum += rets[k] * rets[k]
        dk = ((cum / c_all) - ((k + 1).to_f / n)).abs
        stat = Math.sqrt(0.5 * n) * dk
        if stat > best_val
          best_val = stat
          best_k = k + 1
        end
      end
      { cusum_stat: best_val, detected_idx: best_k, critical_5pct: 1.358, break_detected: best_val > 1.358 }
    end

    def run_structural_break_study(seed: MAIN_SEED, workers: nil, n_reps: 300)
      g = Garch11.new
      norm = Innovations::Normal.new
      rows = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 5) do |idx|
        sim = simulate_break_series(n_obs: 3000, seed: Rng.derive_seed(seed, 800, idx))
        rets = sim[:returns]
        b = sim[:break_idx]
        f_full = Optim.fit(g, rets, dist: norm, max_iter: 120)
        f_pre = Optim.fit(g, rets[0...b], dist: norm, max_iter: 120)
        f_post = Optim.fit(g, rets[b..], dist: norm, max_iter: 120)
        it_test = inclan_tiao_cusum(rets)
        { full_p: g.persistence(f_full[:params]), pre_p: g.persistence(f_pre[:params]),
          post_p: g.persistence(f_post[:params]), cusum: it_test[:cusum_stat],
          tau: it_test[:detected_idx], detected: it_test[:break_detected] }
      end
      summarize_break_mc(rows)
    end

    def summarize_break_mc(rows)
      n = rows.length.to_f
      full_vals = rows.map { |r| r[:full_p] }
      sub_vals = rows.map { |r| 0.5 * (r[:pre_p] + r[:post_p]) }
      full_dist = dist_stats(full_vals)
      sub_dist = dist_stats(sub_vals)
      { true_persistence: 0.85, n_reps: rows.length,
        full_sample_persistence_mean: full_dist[:mean], full_sample_persistence_se: full_dist[:se],
        subsample_persistence_mean: sub_dist[:mean], subsample_persistence_se: sub_dist[:se],
        full_persistence_dist: full_dist, subsample_persistence_dist: sub_dist,
        spurious_inflation: full_dist[:mean] - sub_dist[:mean],
        cusum_power: rows.count { |r| r[:detected] } / n,
        detected_idx_mean: rows.sum(0.0) { |r| r[:tau] } / n }
    end
    private_class_method :summarize_break_mc

    def dist_stats(vals)
      s = vals.sort
      n = s.length.to_f
      mean = s.sum(0.0) / n
      sd = Math.sqrt(s.sum(0.0) { |v| (v - mean)**2 } / [n - 1.0, 1.0].max)
      { mean: mean, sd: sd, se: sd / Math.sqrt(n), min: s.first,
        p05: s[(0.05 * n).floor], p25: s[(0.25 * n).floor], median: s[(0.50 * n).floor],
        p75: s[(0.75 * n).floor], p95: s[(0.95 * n).floor], max: s.last }
    end
    private_class_method :dist_stats

    def run_five_asset_portfolio(seed: MAIN_SEED)
      sim = simulate_five_asset_dcc(n_obs: 2500, seed: Rng.derive_seed(seed, 850, 0))
      train_end = 1500
      rets = sim[:returns]
      gjr = GjrGarch11.new
      t_fit = Innovations::StudentT.new(nu: 6.0)
      fits, z_mat, h_mat = fit_five_marginals(gjr, t_fit, rets, train_end)
      d_fit = Dcc.fit_stage2(z_mat[0...train_end])
      r_seq = Dcc.filter_corr(d_fit[:a], d_fit[:b], z_mat, q_bar: d_fit[:q_bar])[:corr_seq]
      oos_res = evaluate_portfolio_oos(rets, h_mat, r_seq, fits, train_end)
      avg_nu = fits.sum(0.0) { |f| f[:nu] } / 5.0
      stress = evaluate_crisis_stress(h_mat[train_end - 1], r_seq[train_end - 1], avg_nu)
      { dcc_fit: { a: d_fit[:a], b: d_fit[:b] }, oos_backtest: oos_res, stress_test: stress }
    end

    def simulate_five_asset_dcc(n_obs:, seed:)
      rng = Rng::Sampler.new(seed)
      k = 5
      qb = Array.new(k) { |i| Array.new(k) { |j| i == j ? 1.0 : 0.35 } }
      q = qb.map(&:dup)
      gjr = GjrGarch11.new
      h = PORT_GJR.map { |p| gjr.unconditional_variance(p) }
      eps = Array.new(k, 0.0)
      z = Array.new(k, 0.0)
      out = Array.new(n_obs)
      (n_obs + 300).times do |t|
        r_mat = Dcc.normalize_corr(q)
        l_mat = Dcc.cholesky_lower(r_mat)
        u = Array.new(k) { rng.next_std_t(6.0) }
        z = Array.new(k) { |i| (0..i).sum(0.0) { |j| l_mat[i][j] * u[j] } }
        eps = Array.new(k) { |i| Math.sqrt(h[i]) * z[i] }
        out[t - 300] = Array.new(k) { |i| PORT_GJR[i][:mu] + eps[i] } if t >= 300
        h, q = step_five_asset(h, q, qb, eps, z, t)
      end
      { returns: out }
    end

    def step_five_asset(h, q, qb, eps, z, t_idx)
      k = 5
      h_new = Array.new(k) do |i|
        p = PORT_GJR[i]
        ind = eps[i].negative? ? 1.0 : 0.0
        p[:omega] + ((p[:alpha] + (p[:gamma] * ind)) * eps[i] * eps[i]) + (p[:beta] * h[i])
      end
      # Introduce crisis correlation persistence bump in final OOS segment
      a_d = t_idx > 1800 ? 0.08 : 0.04
      b_d = 0.90
      c_d = 1.0 - a_d - b_d
      q_new = Array.new(k) { |i| Array.new(k) { |j| (c_d * qb[i][j]) + (a_d * z[i] * z[j]) + (b_d * q[i][j]) } }
      [h_new, q_new]
    end
    private_class_method :step_five_asset

    def fit_five_marginals(gjr, t_fit, rets, train_end)
      k = 5
      t_tot = rets.length
      fits = []
      h_mat = Array.new(t_tot) { Array.new(k, 0.0) }
      z_mat = Array.new(t_tot) { Array.new(k, 0.0) }
      k.times do |i|
        col = rets.map { |row| row[i] }
        fit = Optim.fit(gjr, col[0...train_end], dist: t_fit, init_params: PORT_GJR[i], max_iter: 100)[:params]
        fits << fit
        vars = gjr.filter(fit, col)
        t_tot.times do |t|
          h_mat[t][i] = vars[t]
          z_mat[t][i] = (col[t] - fit[:mu]) / Math.sqrt(vars[t])
        end
      end
      [fits, z_mat, h_mat]
    end
    private_class_method :fit_five_marginals

    def evaluate_portfolio_oos(rets, h_mat, r_seq, fits, train_end)
      w = PORT_WEIGHTS
      static_cov = sample_cov(rets[0...train_end])
      static_var = Linalg.quad_form(w, static_cov)
      mu_p = w.each_with_index.sum(0.0) { |wi, i| wi * fits[i][:mu] }
      losses = rets[train_end..].map { |row| -w.each_with_index.sum(0.0) { |wi, i| wi * row[i] } }
      avg_nu = fits.sum(0.0) { |f| f[:nu] } / 5.0
      td = Innovations::StudentT.new(nu: avg_nu)
      nd = Innovations::Normal.new
      res_dyn = eval_dynamic_port(losses, h_mat[train_end..], r_seq[train_end..], w, mu_p, td)
      res_stat = eval_static_port(losses, Math.sqrt(static_var), mu_p, nd)
      { dynamic_gjr_dcc: res_dyn, static_normal: res_stat }
    end
    private_class_method :evaluate_portfolio_oos

    def eval_dynamic_port(losses, h_oos, r_oos, w, mu_p, td)
      sig_t = Array.new(losses.length) do |i|
        s = Array.new(5) { |r| Array.new(5) { |c| r_oos[i][r][c] * Math.sqrt(h_oos[i][r] * h_oos[i][c]) } }
        Math.sqrt(Linalg.quad_form(w, s))
      end
      { var1_p: run_var_es(losses, sig_t, mu_p, 0.01, td), var5_p: run_var_es(losses, sig_t, mu_p, 0.05, td) }
    end
    private_class_method :eval_dynamic_port

    def eval_static_port(losses, sig_c, mu_p, nd)
      sig_arr = Array.new(losses.length, sig_c)
      { var1_p: run_var_es(losses, sig_arr, mu_p, 0.01, nd), var5_p: run_var_es(losses, sig_arr, mu_p, 0.05, nd) }
    end
    private_class_method :eval_static_port

    def run_var_es(losses, sigmas, mu_p, alpha, dist)
      q = dist.var_quantile(alpha)
      es_m = dist.expected_shortfall(alpha)
      vars = sigmas.map { |s| (-mu_p) + (s * q) }
      ess = sigmas.map { |s| (-mu_p) + (s * es_m) }
      az = RiskBacktest.acerbi_szekely_z2(losses, vars, ess, alpha)
      RiskBacktest.christoffersen_test(az[:hits], alpha).merge(acerbi_z2: az[:z2])
    end
    private_class_method :run_var_es

    def evaluate_crisis_stress(h_calm, r_calm, avg_nu)
      w = PORT_WEIGHTS
      cov_calm = Array.new(5) { |i| Array.new(5) { |j| r_calm[i][j] * Math.sqrt(h_calm[i] * h_calm[j]) } }
      cov_stress = Array.new(5) do |i|
        Array.new(5) { |j| (i == j ? 1.0 : 0.85) * 6.25 * Math.sqrt(h_calm[i] * h_calm[j]) }
      end
      calm_rep = euler_risk_decomp(w, cov_calm, avg_nu)
      stress_rep = euler_risk_decomp(w, cov_stress, avg_nu)
      { calm: calm_rep, stress: stress_rep,
        var1_amplification: stress_rep[:port_var1] / calm_rep[:port_var1],
        div_ratio_drop: calm_rep[:div_ratio] - stress_rep[:div_ratio] }
    end
    private_class_method :evaluate_crisis_stress

    def euler_risk_decomp(w, cov, avg_nu)
      sig_w = Array.new(5) { |i| Linalg.dot(cov[i], w) }
      port_sd = Math.sqrt(Linalg.dot(w, sig_w))
      rc = Array.new(5) { |i| w[i] * sig_w[i] / port_sd }
      weighted_sd = (0..4).sum(0.0) { |i| w[i] * Math.sqrt(cov[i][i]) }
      td = Innovations::StudentT.new(nu: avg_nu)
      { port_sd: port_sd, euler_sum_err: (rc.sum(0.0) - port_sd).abs,
        rc_pct: rc.map { |v| v / port_sd }, div_ratio: weighted_sd / port_sd,
        port_var1: port_sd * td.var_quantile(0.01), port_es1: port_sd * td.expected_shortfall(0.01) }
    end
    private_class_method :euler_risk_decomp

    def sample_cov(mat)
      n = mat.length.to_f
      means = Array.new(5) { |j| mat.sum(0.0) { |row| row[j] } / n }
      Array.new(5) do |i|
        Array.new(5) { |j| mat.sum(0.0) { |r| (r[i] - means[i]) * (r[j] - means[j]) } / (n - 1.0) }
      end
    end
    private_class_method :sample_cov
  end
end
