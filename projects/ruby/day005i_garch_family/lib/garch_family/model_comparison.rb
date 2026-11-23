# frozen_string_literal: true

module GarchFamily
  module ModelComparison
    SYM_TRUE = { mu: 0.0, omega: 5.0e-6, alpha: 0.08, beta: 0.85 }.freeze
    ASYM_TRUE = { mu: 0.0, omega: 5.0e-6, alpha: 0.03, gamma: 0.10, beta: 0.85 }.freeze

    module_function

    def qlike_loss(h_fc, proxy_sq)
      h = [h_fc, 1.0e-15].max
      r = [proxy_sq, 1.0e-18].max / h
      r - Math.log(r) - 1.0
    end

    def mse_loss(h_fc, proxy_sq) = (h_fc - proxy_sq)**2

    def diebold_mariano(loss_a, loss_b, lags: nil)
      n = loss_a.length
      d = Array.new(n) { |i| loss_a[i] - loss_b[i] }
      mean_d = d.sum(0.0) / n
      nw = Linalg.ols_with_nw(Array.new(n) { [1.0] }, d, lags: lags)
      se = nw[:se_nw][0]
      t_stat = se.positive? ? (mean_d / se) : 0.0
      { mean_diff: mean_d, se_nw: se, dm_stat: t_stat, p_value: 2.0 * (1.0 - Special.norm_cdf(t_stat.abs)) }
    end

    def run_study(seed: MAIN_SEED, workers: nil, n_reps: 200)
      anchor = qlike_anchor_test(seed)
      dgp_sym = compare_on_dgp(Garch11.new, SYM_TRUE, stream_id: 500, seed: seed)
      dgp_asym = compare_on_dgp(GjrGarch11.new, ASYM_TRUE, stream_id: 501, seed: seed)
      lev_mc = leverage_test_mc(seed: seed, workers: workers, n_reps: n_reps)
      { qlike_anchor: anchor, symmetric_dgp: dgp_sym, asymmetric_dgp: dgp_asym, leverage_mc: lev_mc }
    end

    def qlike_anchor_test(seed)
      g = Garch11.new
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 510, 0))
      sim = g.simulate(SYM_TRUE, 40_000, Innovations::Normal.new, rng)
      sq = sim[:returns].map { |r| r * r }
      ht = sim[:true_vars]
      v_bar = g.unconditional_variance(SYM_TRUE)
      losses = {
        true_h: ht.each_with_index.sum(0.0) { |h, i| qlike_loss(h, sq[i]) } / sq.length,
        scaled_up_h: ht.each_with_index.sum(0.0) { |h, i| qlike_loss(1.25 * h, sq[i]) } / sq.length,
        scaled_dn_h: ht.each_with_index.sum(0.0) { |h, i| qlike_loss(0.80 * h, sq[i]) } / sq.length,
        const_vbar: sq.sum(0.0) { |r2| qlike_loss(v_bar, r2) } / sq.length
      }
      losses[:true_is_lowest] = losses[:true_h] < [losses[:scaled_up_h], losses[:scaled_dn_h], losses[:const_vbar]].min
      losses
    end

    def compare_on_dgp(dgp_model, dgp_params, stream_id:, seed:)
      dist = Innovations::Normal.new
      rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, 0))
      rets = dgp_model.simulate(dgp_params, 3500, dist, rng)[:returns]
      train = rets[0, 2500]
      oos = rets[2500, 1000]
      fits = fit_family(train, dist)
      oos_eval = evaluate_oos_models(fits, train, oos)
      nic_x = (-12..12).map { |k| k * 0.0025 }
      { in_sample: summarize_in_sample(fits, train.length), oos_qlike: oos_eval[:qlike],
        oos_mse: oos_eval[:mse], dm_garch_vs_gjr: oos_eval[:dm],
        nic_grid: nic_x, nic_curves: compute_nic(fits, nic_x) }
    end

    def fit_family(train, dist)
      models = { garch11: Garch11.new, gjr: GjrGarch11.new, egarch: Egarch11.new }
      models.transform_values { |m| { model: m, fit: Optim.fit(m, train, dist: dist, max_iter: 120) } }
    end

    def summarize_in_sample(fits, n_obs)
      norm = Innovations::Normal.new
      fits.transform_values do |entry|
        m = entry[:model]
        ll = entry[:fit][:loglik]
        k = m.param_names(norm).length
        { loglik: ll, k: k, aic: (-2.0 * ll) + (2.0 * k),
          bic: (-2.0 * ll) + (k * Math.log(n_obs)), params: entry[:fit][:params] }
      end
    end
    private_class_method :summarize_in_sample

    def evaluate_oos_models(fits, train, oos)
      full = train + oos
      t0 = train.length
      oos_sq = oos.map { |r| r * r }
      q_series = {}
      m_series = {}
      fits.each do |key, entry|
        hv = entry[:model].filter(entry[:fit][:params], full)[t0..]
        q_series[key] = hv.each_with_index.map { |h, i| qlike_loss(h, oos_sq[i]) }
        m_series[key] = hv.each_with_index.map { |h, i| mse_loss(h, oos_sq[i]) }
      end
      ew_h = Ewma.new.filter({}, full)[t0..]
      q_series[:ewma] = ew_h.each_with_index.map { |h, i| qlike_loss(h, oos_sq[i]) }
      m_series[:ewma] = ew_h.each_with_index.map { |h, i| mse_loss(h, oos_sq[i]) }
      { qlike: q_series.transform_values { |v| v.sum(0.0) / v.length },
        mse: m_series.transform_values { |v| v.sum(0.0) / v.length },
        dm: diebold_mariano(q_series[:garch11], q_series[:gjr]) }
    end

    def compute_nic(fits, nic_x)
      base_v = fits[:garch11][:model].unconditional_variance(fits[:garch11][:fit][:params])
      fits.transform_values do |entry|
        vals = entry[:model].news_impact(entry[:fit][:params], nic_x, base_var: base_v)
        { values: vals, ratio_neg_pos_2sd: vals.first / [vals.last, 1.0e-15].max }
      end
    end
    private_class_method :compute_nic

    def leverage_test_mc(seed:, workers:, n_reps:)
      { size_h0_garch: run_lev_mc_cell(Garch11.new, SYM_TRUE, 1500, n_reps, 530, seed, workers),
        power_h1_gjr: run_lev_mc_cell(GjrGarch11.new, ASYM_TRUE, 1500, n_reps, 531, seed, workers) }
    end

    def run_lev_mc_cell(dgp_model, dgp_p, n_obs, n_reps, stream_id, seed, workers)
      dist = Innovations::Normal.new
      g11 = Garch11.new
      gjr = GjrGarch11.new
      rows = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 15) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream_id, idx))
        rets = dgp_model.simulate(dgp_p, n_obs, dist, rng)[:returns]
        single_lev_test(g11, gjr, rets, dist)
      end
      n = rows.length.to_f
      lr_rate = rows.count { |r| r[:lr_p] < 0.05 } / n
      w_rate = rows.count { |r| r[:wald_p] < 0.05 } / n
      { n_reps: rows.length, lr_rej_5pct: lr_rate, lr_se: Math.sqrt(lr_rate * (1.0 - lr_rate) / n),
        wald_rej_5pct: w_rate, wald_se: Math.sqrt(w_rate * (1.0 - w_rate) / n) }
    end
    private_class_method :run_lev_mc_cell

    def single_lev_test(g11, gjr, rets, dist)
      f0 = Optim.fit(g11, rets, dist: dist, init_params: SYM_TRUE, max_iter: 90)
      init_j = f0[:params].merge(gamma: 0.002)
      extra_j = [f0[:params].merge(gamma: 0.04), f0[:params].merge(gamma: 0.08)]
      f1 = Optim.fit(gjr, rets, dist: dist, init_params: init_j, extra_starts: extra_j, max_iter: 110)
      lr = Inference.lr_test(f1[:loglik], f0[:loglik], mixture_zero_one: true)
      rep = Inference.covariance_report(gjr, f1[:params], rets, dist)
      w_p = one_sided_gamma_wald_p(f1[:params][:gamma], rep)
      { lr_p: lr[:p_value], wald_p: w_p }
    end
    private_class_method :single_lev_test

    def one_sided_gamma_wald_p(gamma_hat, rep)
      return 1.0 unless rep

      se = Math.sqrt([rep[:cov_sandwich][3][3], 0.0].max)
      return 1.0 unless se.positive? && se.finite?

      z = gamma_hat / se
      gamma_hat > 1.0e-5 ? 1.0 - Special.norm_cdf(z) : 1.0
    end
    private_class_method :one_sided_gamma_wald_p
  end
end
