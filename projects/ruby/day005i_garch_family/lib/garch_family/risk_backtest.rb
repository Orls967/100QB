# frozen_string_literal: true

module GarchFamily
  module RiskBacktest
    ALPHA_LEVELS = [0.01, 0.05].freeze

    module_function

    def kupiec_test(hits, alpha)
      t = hits.length
      x = hits.count(true)
      p_hat = x.to_f / t
      ll_alt = safe_binom_ll(x, t, p_hat)
      ll_null = safe_binom_ll(x, t, alpha)
      lr = [2.0 * (ll_alt - ll_null), 0.0].max
      { violations: x, n_obs: t, hit_rate: p_hat, lr_uc: lr, p_uc: Special.chi2_sf(lr, 1) }
    end

    def christoffersen_test(hits, alpha)
      kup = kupiec_test(hits, alpha)
      c = transition_counts(hits)
      n0 = c[:n00] + c[:n01]
      n1 = c[:n10] + c[:n11]
      pi01 = n0.positive? ? c[:n01].to_f / n0 : 0.0
      pi11 = n1.positive? ? c[:n11].to_f / n1 : 0.0
      pi_all = (c[:n01] + c[:n11]).to_f / [n0 + n1, 1].max
      ll_ind = safe_binom_ll(c[:n01] + c[:n11], n0 + n1, pi_all)
      ll_mark = safe_binom_ll(c[:n01], n0, pi01) + safe_binom_ll(c[:n11], n1, pi11)
      lr_ind = [2.0 * (ll_mark - ll_ind), 0.0].max
      lr_cc = kup[:lr_uc] + lr_ind
      kup.merge(lr_ind: lr_ind, p_ind: Special.chi2_sf(lr_ind, 1),
                lr_cc: lr_cc, p_cc: Special.chi2_sf(lr_cc, 2), transitions: c)
    end

    def acerbi_szekely_z2(losses, vars, ess, alpha)
      t = losses.length
      sum_ratio = 0.0
      hits = Array.new(t, false)
      t.times do |i|
        next unless losses[i] > vars[i]

        hits[i] = true
        sum_ratio += losses[i] / [ess[i], 1.0e-12].max
      end
      z2 = 1.0 - (sum_ratio / (alpha * t))
      { z2: z2, hits: hits }
    end

    def run_study(seed: MAIN_SEED, workers: nil)
      rolling = rolling_oos_study(seed: seed, n_total: 4000, win: 1000, refit_every: 20, workers: workers)
      mc = mc_size_power_study(seed: seed, workers: workers, n_reps: 500, n_oos: 1000)
      { rolling_oos: rolling, mc_size_power: mc }
    end

    def rolling_oos_study(seed:, n_total:, win:, refit_every:, workers: nil)
      model = Garch11.new
      t5 = Innovations::StudentT.new(nu: 5.0)
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 400, 0))
      rets = model.simulate(InnovationStudy::TRUE_T_PARAMS, n_total, t5, rng)[:returns]
      fc = build_rolling_forecasts(rets, win, refit_every, workers: workers)
      losses = rets[win...n_total].map(&:-@)
      results = {}
      fc.each do |mname, series|
        ALPHA_LEVELS.each do |alpha|
          results["#{mname}_#{(alpha * 100).to_i}pct"] = eval_oos_series(losses, series, alpha)
        end
      end
      results
    end

    def build_rolling_forecasts(rets, win, refit_every, workers: nil)
      names = %i[const_norm ewma_norm garch_norm garch_t garch_fhs]
      out = names.to_h { |k| [k, []] }
      snap_times = (win...rets.length).step(refit_every).to_a
      snaps = Parallel.map_items(n_items: snap_times.length, workers: workers, chunk_size: 8) do |idx|
        t0 = snap_times[idx]
        refit_all_models(rets[(t0 - win)...t0])
      end
      replay_rolling_steps(rets, win, refit_every, snaps, out)
      out
    end
    private_class_method :build_rolling_forecasts

    def replay_rolling_steps(rets, win, refit_every, snaps, out)
      cur_fits = nil
      (win...rets.length).each do |t|
        offset = t - win
        cur_fits = deep_dup_fit(snaps[offset / refit_every]) if (offset % refit_every).zero?
        append_one_step!(out, rets[t - 1], cur_fits)
      end
    end
    private_class_method :replay_rolling_steps

    def deep_dup_fit(fit)
      fit.merge(ewma_v: fit[:ewma_v], hn: fit[:hn], ht: fit[:ht])
    end
    private_class_method :deep_dup_fit

    def refit_all_models(win_rets)
      g = Garch11.new
      mu_c = win_rets.sum(0.0) / win_rets.length
      var_c = win_rets.sum(0.0) { |r| (r - mu_c)**2 } / win_rets.length
      ew_v = Ewma.new(lambda_decay: 0.94).filter({}, win_rets).last
      fn = Optim.fit(g, win_rets, dist: Innovations::Normal.new, max_iter: 80)[:params]
      td = Innovations::StudentT.new(nu: 6.0)
      ft = Optim.fit(g, win_rets, dist: td, init_params: fn.merge(nu: 6.0), max_iter: 80)[:params]
      f_vars = g.filter(fn, win_rets)
      z_sorted = Array.new(win_rets.length) { |i| (win_rets[i] - fn[:mu]) / Math.sqrt(f_vars[i]) }.sort
      fhs = ALPHA_LEVELS.to_h { |al| [al, fhs_tail_moments(z_sorted, al)] }
      { const: [mu_c, var_c], ewma_v: ew_v, fn: fn, ft: ft, hn: f_vars.last,
        ht: g.filter(ft, win_rets).last, fhs: fhs }
    end
    private_class_method :refit_all_models

    def append_one_step!(out, prev_r, fits)
      out[:const_norm] << { mu: fits[:const][0], sigma: Math.sqrt(fits[:const][1]), kind: :norm }
      ew_next = (0.94 * fits[:ewma_v]) + (0.06 * prev_r * prev_r)
      fits[:ewma_v] = ew_next
      out[:ewma_norm] << { mu: 0.0, sigma: Math.sqrt(ew_next), kind: :norm }
      append_garch_steps!(out, prev_r, fits)
    end
    private_class_method :append_one_step!

    def append_garch_steps!(out, prev_r, fits)
      fits[:hn] = step_var(fits[:fn], prev_r, fits[:hn])
      fits[:ht] = step_var(fits[:ft], prev_r, fits[:ht])
      sig_n = Math.sqrt(fits[:hn])
      out[:garch_norm] << { mu: fits[:fn][:mu], sigma: sig_n, kind: :norm }
      out[:garch_t] << { mu: fits[:ft][:mu], sigma: Math.sqrt(fits[:ht]), kind: :t, nu: fits[:ft][:nu] }
      out[:garch_fhs] << { mu: fits[:fn][:mu], sigma: sig_n, kind: :fhs, fhs: fits[:fhs] }
    end
    private_class_method :append_garch_steps!

    def step_var(p, prev_r, prev_var)
      ep = prev_r - p[:mu]
      p[:omega] + (p[:alpha] * ep * ep) + (p[:beta] * prev_var)
    end
    private_class_method :step_var

    def fhs_tail_moments(sorted_z, alpha)
      k = [(alpha * sorted_z.length).floor, 1].max
      [-sorted_z[k - 1], -sorted_z[0...k].sum(0.0) / k]
    end
    private_class_method :fhs_tail_moments

    def eval_oos_series(losses, series, alpha)
      norm = Innovations::Normal.new
      tdist = Innovations::StudentT.new(nu: 5.0)
      vars = []
      ess = []
      series.each do |s|
        v, e = var_es_from_step(s, alpha, norm, tdist)
        vars << v
        ess << e
      end
      az = acerbi_szekely_z2(losses, vars, ess, alpha)
      christoffersen_test(az[:hits], alpha).merge(acerbi_z2: az[:z2])
    end
    private_class_method :eval_oos_series

    def var_es_from_step(s, alpha, norm, tdist)
      case s[:kind]
      when :norm
        [(-s[:mu]) + (s[:sigma] * norm.var_quantile(alpha)),
         (-s[:mu]) + (s[:sigma] * norm.expected_shortfall(alpha))]
      when :t
        [(-s[:mu]) + (s[:sigma] * tdist.var_quantile(alpha, s[:nu])),
         (-s[:mu]) + (s[:sigma] * tdist.expected_shortfall(alpha, s[:nu]))]
      when :fhs
        q_loss, es_loss = s[:fhs][alpha]
        [(-s[:mu]) + (s[:sigma] * q_loss), (-s[:mu]) + (s[:sigma] * es_loss)]
      end
    end
    private_class_method :var_es_from_step

    def mc_size_power_study(seed:, workers:, n_reps:, n_oos:)
      size_hits = Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 25) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 420, idx))
        hits = Array.new(n_oos) { rng.next_float < 0.05 }
        res = christoffersen_test(hits, 0.05)
        { uc: res[:p_uc] < 0.05, ind: res[:p_ind] < 0.05, cc: res[:p_cc] < 0.05 }
      end
      power_hits = run_power_mc(seed: seed, workers: workers, n_reps: n_reps, n_oos: n_oos)
      { size_nominal_5pct: summarize_rates(size_hits), power_misspecified_t5_const: summarize_rates(power_hits) }
    end

    def run_power_mc(seed:, workers:, n_reps:, n_oos:)
      g = Garch11.new
      t5 = Innovations::StudentT.new(nu: 5.0)
      q5 = Special.norm_inv(0.95)
      Parallel.map_items(n_items: n_reps, workers: workers, chunk_size: 25) do |idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 425, idx))
        rets = g.simulate(InnovationStudy::TRUE_T_PARAMS, n_oos, t5, rng)[:returns]
        sd = Math.sqrt(rets.sum(0.0) { |r| r * r } / n_oos)
        hits = rets.map { |r| -r > q5 * sd }
        res = christoffersen_test(hits, 0.05)
        { uc: res[:p_uc] < 0.05, ind: res[:p_ind] < 0.05, cc: res[:p_cc] < 0.05 }
      end
    end
    private_class_method :run_power_mc

    def summarize_rates(rows)
      n = rows.length.to_f
      %i[uc ind cc].each_with_object({}) do |k, acc|
        rate = rows.count { |r| r[k] } / n
        acc[k] = { rate: rate, se: Math.sqrt(rate * (1.0 - rate) / n) }
      end
    end
    private_class_method :summarize_rates

    def transition_counts(hits)
      c = { n00: 0, n01: 0, n10: 0, n11: 0 }
      (1...hits.length).each do |i|
        key = :"n#{hits[i - 1] ? 1 : 0}#{hits[i] ? 1 : 0}"
        c[key] += 1
      end
      c
    end
    private_class_method :transition_counts

    def safe_binom_ll(k, n, p)
      return 0.0 if n.zero? || p <= 0.0 || p >= 1.0

      a = k.positive? ? k * Math.log(p) : 0.0
      b = (n - k).positive? ? (n - k) * Math.log(1.0 - p) : 0.0
      a + b
    end
    private_class_method :safe_binom_ll
  end
end
