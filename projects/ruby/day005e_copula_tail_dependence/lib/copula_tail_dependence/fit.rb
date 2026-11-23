# frozen_string_literal: true

module CopulaTailDependence
  module Fit
    PHI = (Math.sqrt(5.0) - 1.0) * 0.5

    module_function

    def pseudo_observations(xs)
      n = xs.length
      ranks = Dependence.average_ranks(xs)
      denom = n + 1.0
      ranks.map { |r| r / denom }
    end

    def golden_max(lo, hi, tol: 1.0e-5, &obj)
      a = lo
      b = hi
      c = b - (PHI * (b - a))
      d = a + (PHI * (b - a))
      fc = obj.call(c)
      fd = obj.call(d)
      while (b - a) > tol
        if fc > fd
          b = d
          d = c
          fd = fc
          c = b - (PHI * (b - a))
          fc = obj.call(c)
        else
          a = c
          c = d
          fc = fd
          d = a + (PHI * (b - a))
          fd = obj.call(d)
        end
      end
      best_x = 0.5 * (a + b)
      [best_x, obj.call(best_x)]
    end

    def fit_gaussian_mle(u1, u2)
      z1 = u1.map { |u| Special.norm_inv(u.clamp(1.0e-12, 1.0 - 1.0e-12)) }
      z2 = u2.map { |u| Special.norm_inv(u.clamp(1.0e-12, 1.0 - 1.0e-12)) }
      n = z1.length
      s11 = 0.0
      s22 = 0.0
      s12 = 0.0
      i = 0
      while i < n
        s11 += z1[i] * z1[i]
        s22 += z2[i] * z2[i]
        s12 += z1[i] * z2[i]
        i += 1
      end
      ll_fn = lambda do |rho|
        r2 = rho * rho
        (-0.5 * n * Math.log(1.0 - r2)) - (((r2 * (s11 + s22)) - (2.0 * rho * s12)) / (2.0 * (1.0 - r2)))
      end
      rho_hat, ll = golden_max(-0.99, 0.99, &ll_fn)
      { rho: rho_hat, loglik: ll, aic: 2.0 - (2.0 * ll) }
    end

    def fit_t_for_df(u1, u2, df)
      return fit_gaussian_mle(u1, u2).merge(df: df) if df.infinite?

      n = u1.length
      t1 = u1.map { |u| Special.t_inv(u.clamp(1.0e-12, 1.0 - 1.0e-12), df) }
      t2 = u2.map { |u| Special.t_inv(u.clamp(1.0e-12, 1.0 - 1.0e-12), df) }
      c_const = Special.log_gamma(0.5 * (df + 2.0)) + Special.log_gamma(0.5 * df) -
                (2.0 * Special.log_gamma(0.5 * (df + 1.0)))
      marg = t1.each_with_index.sum(0.0) do |a, i|
        b = t2[i]
        0.5 * (df + 1.0) * (Math.log(1.0 + ((a * a) / df)) + Math.log(1.0 + ((b * b) / df)))
      end
      rho_hat, ll = golden_max(-0.99, 0.99, tol: 1.0e-4) do |rho|
        t_bivariate_loglik(n, t1, t2, rho, df, c_const, marg)
      end
      { df: df, rho: rho_hat, loglik: ll, aic: 4.0 - (2.0 * ll) }
    end
    private_class_method :fit_t_for_df

    def t_bivariate_loglik(n, t1, t2, rho, df, c_const, marg)
      r2 = rho * rho
      sum_q = 0.0
      i = 0
      while i < n
        a = t1[i]
        b = t2[i]
        q = ((a * a) - (2.0 * rho * a * b) + (b * b)) / (df * (1.0 - r2))
        sum_q += Math.log(1.0 + q)
        i += 1
      end
      (n * (c_const - (0.5 * Math.log(1.0 - r2)))) - (0.5 * (df + 2.0) * sum_q) + marg
    end
    private_class_method :t_bivariate_loglik

    def fit_t_profile_mle(u1, u2, df_grid: [2.5, 3.0, 4.0, 6.0, 10.0, 30.0, Float::INFINITY])
      profiles = df_grid.map { |df| fit_t_for_df(u1, u2, df) }
      best = profiles.max_by { |r| r[:loglik] }
      { best: best, profiles: profiles }
    end

    def fit_clayton_mle(u1, u2)
      n = u1.length
      lu = u1.map { |u| Math.log(u.clamp(1.0e-12, 1.0 - 1.0e-12)) }
      lv = u2.map { |u| Math.log(u.clamp(1.0e-12, 1.0 - 1.0e-12)) }
      sum_log = lu.sum(0.0) + lv.sum(0.0)
      ll_fn = lambda do |th|
        s = 0.0
        i = 0
        while i < n
          val = Math.exp(-th * lu[i]) + Math.exp(-th * lv[i]) - 1.0
          return -1.0e12 if val <= 0.0

          s += Math.log(val)
          i += 1
        end
        (n * Math.log(1.0 + th)) - ((1.0 + th) * sum_log) - ((2.0 + (1.0 / th)) * s)
      end
      th_hat, ll = golden_max(0.05, 10.0, tol: 1.0e-4, &ll_fn)
      { theta: th_hat, loglik: ll, aic: 2.0 - (2.0 * ll) }
    end

    def compare_tau_vs_mle(rho_true: 0.6, ns: [250, 1000, 5000], r_reps: 500, seed: MAIN_SEED, workers: nil)
      ns.map do |n|
        rows = Parallel.map_chunks(n_items: r_reps, workers: workers, chunk_size: 25) do |cidx, len|
          rng = Rng::Sampler.new(Rng.derive_seed(seed, 300 + n, cidx))
          out = []
          len.times do
            x = Array.new(n, 0.0)
            y = Array.new(n, 0.0)
            s_rho = Math.sqrt(1.0 - (rho_true * rho_true))
            n.times do |i|
              z1 = rng.next_normal
              z2 = (rho_true * z1) + (s_rho * rng.next_normal)
              x[i] = z1
              y[i] = z2
            end
            tau = Dependence.kendall_tau(x, y)
            rho_tau = Math.sin(0.5 * Math::PI * tau)
            u1 = pseudo_observations(x)
            u2 = pseudo_observations(y)
            rho_mle = fit_gaussian_mle(u1, u2)[:rho]
            out << [rho_tau, rho_mle]
          end
          out
        end.flatten(1)

        taus = rows.map(&:first)
        mles = rows.map(&:last)
        {
          n: n,
          r_reps: r_reps,
          rho_tau_bias: Stats.mean(taus) - rho_true,
          rho_tau_sd: Stats.sd(taus),
          rho_mle_bias: Stats.mean(mles) - rho_true,
          rho_mle_sd: Stats.sd(mles)
        }
      end
    end

    def transform_variate(zk, scale, copula, marginal, df)
      u = copula == :t ? Special.t_cdf(zk * scale, df) : Special.norm_cdf(zk)
      u_c = u.clamp(1.0e-12, 1.0 - 1.0e-12)
      return Special.t_inv(u_c, df) if marginal == :t

      copula == :t ? Special.norm_inv(u_c) : zk
    end
    private_class_method :transform_variate

    # ES99 portofolio d=5 ekui-bobot, ekui-korelasi rho pada kopula :gaussian atau :t.
    def equicorr_portfolio_es99(copula:, rho:, df: 4.0, marginal: :normal, n_paths: 40_000, seed: 777)
      d = 5
      w = 1.0 / d
      rng = Rng::Sampler.new(seed)
      losses = Array.new(n_paths, 0.0)
      s_rho = Math.sqrt(1.0 - rho)
      sqrt_rho = Math.sqrt(rho)
      i = 0
      while i < n_paths
        f = rng.next_normal
        scale = copula == :t ? Math.sqrt(df / rng.next_chi2(df)) : 1.0
        port = 0.0
        k = 0
        while k < d
          zk = (sqrt_rho * f) + (s_rho * rng.next_normal)
          port += -w * transform_variate(zk, scale, copula, marginal, df)
          k += 1
        end
        losses[i] = port
        i += 1
      end
      losses.sort!
      tail_start = (0.99 * n_paths).floor
      Stats.mean(losses[tail_start..])
    end

    def sample_bivariate_t(n_obs, rho, df, rng)
      x = Array.new(n_obs, 0.0)
      y = Array.new(n_obs, 0.0)
      s_r = Math.sqrt(1.0 - (rho * rho))
      n_obs.times do |i|
        z1 = rng.next_normal
        z2 = (rho * z1) + (s_r * rng.next_normal)
        sc = Math.sqrt(df / rng.next_chi2(df))
        x[i] = z1 * sc
        y[i] = z2 * sc
      end
      [x, y]
    end
    private_class_method :sample_bivariate_t

    def misspecification_study(rho_true: 0.5, df_true: 4.0, n_obs: 2000, r_reps: 200, seed: MAIN_SEED, workers: nil)
      x0, y0 = sample_bivariate_t(n_obs, rho_true, df_true, Rng::Sampler.new(Rng.derive_seed(seed, 350, 0)))
      u1 = pseudo_observations(x0)
      u2 = pseudo_observations(y0)
      fit_g = fit_gaussian_mle(u1, u2)
      fit_t = fit_t_profile_mle(u1, u2)
      fit_c = fit_clayton_mle(u1, u2)
      lrt_stat = 2.0 * (fit_t[:best][:loglik] - fit_g[:loglik])
      es_true = equicorr_portfolio_es99(copula: :t, rho: rho_true, df: df_true, n_paths: 80_000, seed: seed + 1)

      rep_rows = Parallel.map_chunks(n_items: r_reps, workers: workers, chunk_size: 10) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 360, cidx))
        Array.new(len) do |j|
          xr, yr = sample_bivariate_t(n_obs, rho_true, df_true, rng)
          rho_g = fit_gaussian_mle(pseudo_observations(xr), pseudo_observations(yr))[:rho].clamp(0.05, 0.95)
          equicorr_portfolio_es99(copula: :gaussian, rho: rho_g, n_paths: 10_000, seed: seed + (cidx * 100) + j)
        end
      end.flatten(1)

      es_gauss_mean = Stats.mean(rep_rows)
      {
        single_sample: {
          gaussian: fit_g, t_best: fit_t[:best], t_profiles: fit_t[:profiles], clayton: fit_c,
          lrt_t_vs_gauss: lrt_stat, gauss_rejected_at_5pct: lrt_stat > 3.841
        },
        es99_true_t4: es_true,
        es99_fitted_gauss_mean: es_gauss_mean,
        es99_fitted_gauss_sd: Stats.sd(rep_rows),
        es99_rel_bias: (es_gauss_mean - es_true) / es_true
      }
    end
  end
end
