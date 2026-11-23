# frozen_string_literal: true

module GarchFamily
  module HarRv
    module_function

    def simulate_intraday_sv(n_days:, m_bars: 78, seed: MAIN_SEED, burn: 200)
      rng = Rng::Sampler.new(seed)
      tot = n_days + burn
      iv_days = Array.new(n_days)
      rv_days = Array.new(n_days)
      daily_rets = Array.new(n_days)
      # Two-factor persistent log-variance + mean-reverting intraday volatility
      theta_lv = -9.4
      lv_slow = theta_lv
      lv_fast = theta_lv
      tot.times do |d|
        lv_slow = theta_lv + (0.98 * (lv_slow - theta_lv)) + (0.08 * rng.next_normal)
        lv_fast = theta_lv + (0.70 * (lv_fast - theta_lv)) + (0.18 * rng.next_normal)
        v_bar = Math.exp((0.55 * lv_slow) + (0.45 * lv_fast))
        iv_sum, rv_sum, ret_sum = simulate_day_bars(v_bar, m_bars, rng)
        next if d < burn

        idx = d - burn
        iv_days[idx] = iv_sum
        rv_days[idx] = rv_sum
        daily_rets[idx] = ret_sum
      end
      { iv: iv_days, rv: rv_days, returns: daily_rets }
    end

    def simulate_day_bars(v_bar, m_bars, rng)
      dt = 1.0 / m_bars
      v = v_bar
      iv_sum = rv_sum = ret_sum = 0.0
      m_bars.times do
        v = [v + (6.0 * (v_bar - v) * dt) + (0.5 * v * Math.sqrt(dt) * rng.next_normal), 1.0e-9].max
        r_bar = Math.sqrt(v * dt) * rng.next_normal
        iv_sum += v * dt
        rv_sum += r_bar * r_bar
        ret_sum += r_bar
      end
      [iv_sum, rv_sum, ret_sum]
    end
    private_class_method :simulate_day_bars

    def run_study(seed: MAIN_SEED)
      infill = infill_anchor_study(seed)
      har_cmp = har_vs_garch_study(seed)
      { infill_anchor: infill, har_vs_garch: har_cmp }
    end

    def infill_anchor_study(seed)
      m_grid = [13, 39, 78, 390]
      rows = m_grid.each_with_index.map do |m, idx|
        sim = simulate_intraday_sv(n_days: 2000, m_bars: m, seed: Rng.derive_seed(seed, 700 + idx, 0))
        diffs = sim[:rv].each_with_index.map { |rv, i| rv - sim[:iv][i] }
        mean_d = diffs.sum(0.0) / diffs.length
        sd_d = Math.sqrt(diffs.sum(0.0) { |x| (x - mean_d)**2 } / (diffs.length - 1))
        se_d = sd_d / Math.sqrt(diffs.length)
        mse = diffs.sum(0.0) { |x| x * x } / diffs.length
        { m_bars: m, mean_rv_minus_iv: mean_d, se_diff: se_d, z_unbiased: mean_d.abs / se_d, mse_rv_iv: mse }
      end
      { grid: rows, mse_ratio_13_to_390: rows.first[:mse_rv_iv] / rows.last[:mse_rv_iv] }
    end

    def har_design(rv_series, start_idx, end_idx)
      x_mat = []
      y_vec = []
      (start_idx...end_idx).each do |t|
        rv_d = rv_series[t - 1]
        rv_w = rv_series[(t - 5)...t].sum(0.0) / 5.0
        rv_m = rv_series[(t - 22)...t].sum(0.0) / 22.0
        x_mat << [1.0, rv_d, rv_w, rv_m]
        y_vec << rv_series[t]
      end
      [x_mat, y_vec]
    end

    def har_vs_garch_study(seed)
      sim = simulate_intraday_sv(n_days: 2500, m_bars: 78, seed: Rng.derive_seed(seed, 720, 0))
      train_end = 1500
      x_tr, y_tr = har_design(sim[:rv], 22, train_end)
      har_ols = Linalg.ols_with_nw(x_tr, y_tr)
      x_oos, = har_design(sim[:rv], train_end, 2500)
      har_fc = x_oos.map { |row| [Linalg.dot(row, har_ols[:beta]), 1.0e-12].max }
      g_fc = daily_garch_fc(sim[:returns], train_end)
      rv_naive_fc = sim[:rv][(train_end - 1)...2499]
      ewma_rv_fc = ewma_rv_filter(sim[:rv], 0.94)[train_end..]
      garch_x_fc = fit_and_filter_garch_x(sim[:returns], sim[:rv], train_end)
      iv_oos = sim[:iv][train_end..]
      summarize_har_benchmarks(har_ols, iv_oos, har_fc, g_fc, rv_naive_fc, ewma_rv_fc, garch_x_fc)
    end

    def daily_garch_fc(rets, train_end)
      g = Garch11.new
      g_fit = Optim.fit(g, rets[0, train_end], dist: Innovations::Normal.new, max_iter: 120)
      g.filter(g_fit[:params], rets)[train_end..]
    end
    private_class_method :daily_garch_fc

    def ewma_rv_filter(rv, lam)
      out = Array.new(rv.length, rv[0])
      (1...rv.length).each { |t| out[t] = (lam * out[t - 1]) + ((1.0 - lam) * rv[t - 1]) }
      out
    end
    private_class_method :ewma_rv_filter

    def fit_and_filter_garch_x(rets, rv, train_end)
      r_tr = rets[0, train_end]
      rv_tr = rv[0, train_end]
      v0 = r_tr.sum(0.0) { |x| x * x } / train_end
      th0 = [Math.log(0.10 * v0), Math.log(0.05), Math.log(0.75), Math.log(0.15)]
      opt = Optim.bfgs(th0, ->(th) { garch_x_nll_grad(th, r_tr, rv_tr) }, max_iter: 120, gtol: 1.0e-5)
      filter_garch_x(unpack_gx(opt[:x]), rets, rv)[train_end..]
    end
    private_class_method :fit_and_filter_garch_x

    def unpack_gx(th)
      [Math.exp(th[0].clamp(-25.0, 2.0)), Math.exp(th[1].clamp(-15.0, -0.1)),
       Math.exp(th[2].clamp(-15.0, -0.05)), Math.exp(th[3].clamp(-15.0, 2.0))]
    end
    private_class_method :unpack_gx

    def filter_garch_x(params, rets, rv)
      w, a, b, gx = params
      h = [w / [1.0 - a - b, 0.05].max, 1.0e-9].max
      out = Array.new(rets.length, h)
      (1...rets.length).each do |t|
        h = [w + (a * rets[t - 1] * rets[t - 1]) + (b * h) + (gx * rv[t - 1]), 1.0e-12].max
        out[t] = h
      end
      out
    end
    private_class_method :filter_garch_x

    def garch_x_nll_grad(th, rets, rv)
      f0 = garch_x_avg_nll(th, rets, rv)
      step = 1.0e-5
      grad = Array.new(4) do |j|
        tp = th.dup
        tm = th.dup
        tp[j] += step
        tm[j] -= step
        (garch_x_avg_nll(tp, rets, rv) - garch_x_avg_nll(tm, rets, rv)) / (2.0 * step)
      end
      [f0, grad]
    end
    private_class_method :garch_x_nll_grad

    def garch_x_avg_nll(th, rets, rv)
      hs = filter_garch_x(unpack_gx(th), rets, rv)
      s = 0.0
      rets.each_with_index { |r, t| s += 0.5 * (Math.log(2.0 * Math::PI) + Math.log(hs[t]) + ((r * r) / hs[t])) }
      s / rets.length
    end
    private_class_method :garch_x_avg_nll

    def summarize_har_benchmarks(har_ols, iv_oos, har_fc, g_fc, rv_n, ew_rv, gx_fc)
      ql = { har: mean_ql(har_fc, iv_oos), garch: mean_ql(g_fc, iv_oos),
             rv_naive: mean_ql(rv_n, iv_oos), ewma_rv: mean_ql(ew_rv, iv_oos), garch_x: mean_ql(gx_fc, iv_oos) }
      dm = { garch_minus_har: dm_vs_har(g_fc, har_fc, iv_oos),
             rv_naive_minus_har: dm_vs_har(rv_n, har_fc, iv_oos),
             ewma_rv_minus_har: dm_vs_har(ew_rv, har_fc, iv_oos),
             garch_x_minus_har: dm_vs_har(gx_fc, har_fc, iv_oos) }
      { beta: har_ols[:beta], se_nw: har_ols[:se_nw],
        qlike_har_vs_iv: ql[:har], qlike_garch_vs_iv: ql[:garch],
        qlike_benchmarks: ql, dm_garch_minus_har: dm[:garch_minus_har], dm_benchmarks: dm,
        sample_iv_oos: iv_oos[0, 120], sample_har_fc: har_fc[0, 120], sample_garch_fc: g_fc[0, 120] }
    end
    private_class_method :summarize_har_benchmarks

    def mean_ql(fc, iv_oos)
      fc.each_with_index.sum(0.0) { |h, i| ModelComparison.qlike_loss(h, iv_oos[i]) } / iv_oos.length
    end
    private_class_method :mean_ql

    def dm_vs_har(alt_fc, har_fc, iv_oos)
      la = alt_fc.each_with_index.map { |h, i| ModelComparison.qlike_loss(h, iv_oos[i]) }
      lh = har_fc.each_with_index.map { |h, i| ModelComparison.qlike_loss(h, iv_oos[i]) }
      ModelComparison.diebold_mariano(la, lh)
    end
    private_class_method :dm_vs_har
  end
end
