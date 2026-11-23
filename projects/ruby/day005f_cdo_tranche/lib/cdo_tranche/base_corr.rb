# frozen_string_literal: true

module CdoTranche
  module BaseCorr
    module_function

    def bisection_root(lo, hi, tol: 1.0e-5, &fn)
      flo = fn.call(lo)
      fhi = fn.call(hi)
      return nil if flo * fhi > 0.0

      45.times do
        mid = 0.5 * (lo + hi)
        fmid = fn.call(mid)
        return mid if fmid.abs < tol || (hi - lo) < tol

        if flo * fmid <= 0.0
          hi = mid
          fhi = fmid
        else
          lo = mid
          flo = fmid
        end
      end
      0.5 * (lo + hi)
    end

    # Menghitung EL [0, D] dari harga pasar Student-t(df=4, rho=0.25) lewat kuadratur 2D (M, W).
    def t_copula_equity_el(n_names, hazard, recovery, t, rho, df, detach)
      nodes, weights = Special.gauss_legendre(32)
      half = 7.0
      p_uncond = Models.default_prob(hazard, t)
      lgd = (1.0 - recovery) / n_names
      # Integrasi atas chi2 via transformasi eksponensial atau Gauss-Legendre pada u in (0,1)?
      # Untuk kecepatan & akurasi, gunakan 16 titik GL pada log-chi2 atau 2000 titik QMC?
      # Lebih cepat & deterministik: 20-titik GL pada u in (0.001, 0.999) untuk Z_W, tidak perlu chi2:
      # W ~ chi2_df / df -> pakai 24-titik GL pada M dan 12-titik GL pada u in (0, 1) via Gamma quantile?
      # Kita pakai MC deterministik R=60_000 yang sangat cepat untuk kurva pasar t(4) sintetis!
      rng = Rng::Sampler.new(909_001)
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)
      thresh = Special.t_inv(p_uncond, df)
      sum = 0.0
      r_paths = 25_000
      r_paths.times do
        m = rng.next_normal
        sc = Math.sqrt(df / rng.next_chi2(df))
        p_c = Special.norm_cdf(((thresh / sc) - (sr * m)) / s1r)
        # Binomial expectation approx atau 1-draw binomial? Pakai eksak conditional expectation untuk smooth!
        sum += Asb.homogeneous_binomial_el(n_names, lgd, p_c, 0.0, detach)
      end
      _unused = [nodes, weights, half]
      sum / r_paths
    end

    def calibrate_base_correlation(market_equity_els, detaches, hazard: 0.015, recovery: 0.40, t: 5.0)
      p = Models.default_prob(hazard, t)
      prev_el = 0.0
      detaches.each_with_index.map do |d, idx|
        target = market_equity_els[idx]
        rho_base = bisection_root(0.01, 0.95) do |rho|
          Lhp.equity_el_bivariate(d, p, recovery, rho) - target
        end
        el_implied = rho_base ? Lhp.equity_el_bivariate(d, p, recovery, rho_base) : target
        tranche_el = el_implied - prev_el
        prev_el = el_implied
        { detach: d, target_equity_el: target, base_corr: rho_base, implied_tranche_el: tranche_el }
      end
    end

    def compound_correlation_scan(target_el, attach, detach, hazard: 0.015, recovery: 0.40, t: 5.0)
      p = Models.default_prob(hazard, t)
      grid = (1..19).map { |i| i * 0.05 }
      vals = grid.map do |rho|
        Lhp.equity_el_bivariate(detach, p, recovery,
                                rho) - Lhp.equity_el_bivariate(attach, p, recovery, rho) - target_el
      end
      roots = []
      (0...(grid.length - 1)).each do |i|
        next unless vals[i] * vals[i + 1] <= 0.0

        r = bisection_root(grid[i], grid[i + 1]) do |rho|
          Lhp.equity_el_bivariate(detach, p, recovery, rho) -
            Lhp.equity_el_bivariate(attach, p, recovery, rho) - target_el
        end
        roots << r if r
      end
      roots
    end
  end
end
