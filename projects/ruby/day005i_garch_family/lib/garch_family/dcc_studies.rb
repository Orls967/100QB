# frozen_string_literal: true

module GarchFamily
  module Dcc
    module_function

    def simulate_bivariate_dcc(n_obs, seed:, rho_bar: RHO_BAR, dcc_p: TRUE_DCC, burn: 400)
      rng = Rng::Sampler.new(seed)
      st = init_biv_state(rho_bar)
      r1 = []
      r2 = []
      rhos = []
      (n_obs + burn).times do |t|
        step_biv_dcc!(st, rho_bar, dcc_p, rng, t.positive?)
        next if t < burn

        r1 << st[:ep1]
        r2 << st[:ep2]
        rhos << st[:rho]
      end
      { r1: r1, r2: r2, true_rhos: rhos }
    end

    def init_biv_state(rho_bar)
      g = Garch11.new
      { h1: g.unconditional_variance(M1_PARAMS), h2: g.unconditional_variance(M2_PARAMS),
        q: [[1.0, rho_bar], [rho_bar, 1.0]], ep1: 0.0, ep2: 0.0, z1: 0.0, z2: 0.0, rho: rho_bar }
    end
    private_class_method :init_biv_state

    def step_biv_dcc!(st, rho_bar, dcc_p, rng, advance)
      update_biv_hq!(st, rho_bar, dcc_p) if advance
      q = st[:q]
      rho = (q[0][1] / Math.sqrt([q[0][0] * q[1][1], 1.0e-15].max)).clamp(-0.995, 0.995)
      u1 = rng.next_normal
      u2 = rng.next_normal
      st[:rho] = rho
      st[:z1] = u1
      st[:z2] = (rho * u1) + (Math.sqrt(1.0 - (rho * rho)) * u2)
      st[:ep1] = Math.sqrt(st[:h1]) * st[:z1]
      st[:ep2] = Math.sqrt(st[:h2]) * st[:z2]
    end
    private_class_method :step_biv_dcc!

    def update_biv_hq!(st, rho_bar, dcc_p)
      st[:h1] = M1_PARAMS[:omega] + (M1_PARAMS[:alpha] * st[:ep1] * st[:ep1]) + (M1_PARAMS[:beta] * st[:h1])
      st[:h2] = M2_PARAMS[:omega] + (M2_PARAMS[:alpha] * st[:ep2] * st[:ep2]) + (M2_PARAMS[:beta] * st[:h2])
      a, b = dcc_p.values_at(:a, :b)
      c = 1.0 - a - b
      q = st[:q]
      z1, z2 = st.values_at(:z1, :z2)
      q01 = (c * rho_bar) + (a * z1 * z2) + (b * q[0][1])
      st[:q] = [[(c * 1.0) + (a * z1 * z1) + (b * q[0][0]), q01],
                [q01, (c * 1.0) + (a * z2 * z2) + (b * q[1][1])]]
    end
    private_class_method :update_biv_hq!

    def run_study(seed: MAIN_SEED, workers: nil, n_reps: 200)
      { ccc_pd_anchor: ccc_reduction_and_pd_anchor(seed),
        mc_recovery: mc_dcc_recovery(seed: seed, workers: workers, n_reps: n_reps),
        portfolio_var_study: bivariate_portfolio_study(seed: seed) }
    end

    def ccc_reduction_and_pd_anchor(seed)
      sim = simulate_bivariate_dcc(1500, seed: Rng.derive_seed(seed, 600, 0))
      z_rows = std_biv_resids(sim, M1_PARAMS, M2_PARAMS)
      ccc = filter_corr(0.0, 0.0, z_rows)
      max_ccc_err = ccc[:corr_seq].map { |r| (r[0][1] - ccc[:q_bar][0][1]).abs }.max
      fit = fit_stage2(z_rows)
      f_seq = filter_corr(fit[:a], fit[:b], z_rows, q_bar: fit[:q_bar])[:corr_seq]
      rhos_hat = f_seq.map { |r| r[0][1] }
      { max_ccc_err: max_ccc_err, min_det: f_seq.map { |r| (r[0][0] * r[1][1]) - (r[0][1] * r[1][0]) }.min,
        max_abs_rho: rhos_hat.map(&:abs).max, sample_fit: fit,
        sample_corr_true_vs_filt: Special.pearson_corr(sim[:true_rhos], rhos_hat),
        sample_true_rhos: sim[:true_rhos][0, 120], sample_filt_rhos: rhos_hat[0, 120] }
    end

    def std_biv_resids(sim, p1, p2)
      g = Garch11.new
      v1 = g.filter(p1, sim[:r1])
      v2 = g.filter(p2, sim[:r2])
      Array.new(sim[:r1].length) { |t| [sim[:r1][t] / Math.sqrt(v1[t]), sim[:r2][t] / Math.sqrt(v2[t])] }
    end
    private_class_method :std_biv_resids

    def mc_dcc_recovery(seed:, workers:, n_reps:)
      g = Garch11.new
      norm = Innovations::Normal.new
      rows = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 15) do |idx|
        sim = simulate_bivariate_dcc(2000, seed: Rng.derive_seed(seed, 610, idx))
        f1 = Optim.fit(g, sim[:r1], dist: norm, init_params: M1_PARAMS, max_iter: 80)
        f2 = Optim.fit(g, sim[:r2], dist: norm, init_params: M2_PARAMS, max_iter: 80)
        d_fit = fit_stage2(std_biv_resids(sim, f1[:params], f2[:params]))
        { ok: f1[:converged] && f2[:converged] && d_fit[:converged], a: d_fit[:a], b: d_fit[:b] }
      end
      good = rows.select { |r| r[:ok] }
      { n_reps: n_reps, n_converged: good.length, n_failed: n_reps - good.length,
        a: summarize_scalar(good.map { |r| r[:a] }, TRUE_DCC[:a]),
        b: summarize_scalar(good.map { |r| r[:b] }, TRUE_DCC[:b]) }
    end

    def summarize_scalar(vals, target)
      n = vals.length.to_f
      m = vals.sum(0.0) / n
      rmse = Math.sqrt(vals.sum(0.0) { |v| (v - target)**2 } / n)
      { truth: target, mean: m, bias: m - target, se_mean: rmse / Math.sqrt(n), rmse: rmse }
    end
    private_class_method :summarize_scalar

    def bivariate_portfolio_study(seed:)
      sim = simulate_bivariate_dcc(
        3000, seed: Rng.derive_seed(seed, 620, 0), rho_bar: 0.72, dcc_p: { a: 0.08, b: 0.89 }
      )
      oos_loss, v1_oos, v2_oos, r_dcc, d_fit = fit_biv_oos(sim, 1200)
      eval_biv_var(oos_loss, v1_oos, v2_oos, r_dcc, d_fit)
    end
    private_class_method :bivariate_portfolio_study

    def fit_biv_oos(sim, train)
      g = Garch11.new
      norm = Innovations::Normal.new
      p1 = Optim.fit(g, sim[:r1][0, train], dist: norm)[:params]
      p2 = Optim.fit(g, sim[:r2][0, train], dist: norm)[:params]
      v1 = g.filter(p1, sim[:r1])
      v2 = g.filter(p2, sim[:r2])
      z_all = std_biv_resids(sim, p1, p2)
      d_fit = fit_stage2(z_all[0...train])
      r_dcc = filter_corr(d_fit[:a], d_fit[:b], z_all, q_bar: d_fit[:q_bar])[:corr_seq][train..]
      oos_loss = Array.new(3000 - train) { |i| -0.5 * (sim[:r1][train + i] + sim[:r2][train + i]) }
      [oos_loss, v1[train..], v2[train..], r_dcc, d_fit]
    end
    private_class_method :fit_biv_oos

    def eval_biv_var(oos_loss, v1_oos, v2_oos, r_dcc, d_fit)
      q5 = Innovations::Normal.new.var_quantile(0.05)
      v_dcc = Array.new(oos_loss.length) do |i|
        s12 = r_dcc[i][0][1] * Math.sqrt(v1_oos[i] * v2_oos[i])
        q5 * Math.sqrt(0.25 * (v1_oos[i] + v2_oos[i] + (2.0 * s12)))
      end
      v_unc = Array.new(oos_loss.length) { |i| q5 * Math.sqrt(0.25 * (v1_oos[i] + v2_oos[i])) }
      { dcc_5pct: RiskBacktest.christoffersen_test(oos_loss.each_with_index.map { |l, i| l > v_dcc[i] }, 0.05),
        uncorr_5pct: RiskBacktest.christoffersen_test(oos_loss.each_with_index.map { |l, i| l > v_unc[i] }, 0.05),
        rho_ccc: d_fit[:q_bar][0][1], dcc_a: d_fit[:a], dcc_b: d_fit[:b] }
    end
    private_class_method :eval_biv_var
  end
end
