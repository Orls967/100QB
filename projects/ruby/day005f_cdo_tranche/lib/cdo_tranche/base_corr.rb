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

    # Menghitung EL [0, D] dari pasar Student-t(df, rho) lewat kuadratur 2D eksak (M, U = sqrt(W)).
    # Substitusi u = sqrt(w) menghilangkan singularitas turunan sqrt(w) di w = 0.
    def t_copula_equity_el(n_names, hazard, recovery, t, rho, df, detach, m_gl: 64, w_gl: 32)
      m_nodes, m_weights = Special.gauss_legendre(m_gl)
      u_nodes, u_weights = Special.gauss_legendre(w_gl)
      m_half = 6.5
      u_half = 2.5
      p_uncond = Models.default_prob(hazard, t)
      lgd = (1.0 - recovery) / n_names
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)
      thresh = Special.t_inv(p_uncond, df)
      half_df = 0.5 * df
      log_norm = (half_df * Math.log(half_df)) - Special.log_gamma(half_df)

      total = 0.0
      m_nodes.each_with_index do |xi_m, im|
        m = m_half * xi_m
        pdf_m = Special.norm_pdf(m)
        inner = 0.0
        u_nodes.each_with_index do |xi_u, iu|
          u = u_half * (xi_u + 1.0)
          next if u <= 1.0e-10

          w = u * u
          pdf_u = 2.0 * u * Math.exp(log_norm + ((half_df - 1.0) * Math.log(w)) - (half_df * w))
          p_c = Special.norm_cdf(((thresh * u) - (sr * m)) / s1r)
          inner += u_weights[iu] * pdf_u * Asb.homogeneous_binomial_el(n_names, lgd, p_c, 0.0, detach)
        end
        total += m_weights[im] * pdf_m * (u_half * inner)
      end
      m_half * total
    end

    def t_copula_sensitivity_grid(detaches: [0.03, 0.07, 0.10, 0.15, 0.30])
      [
        { df: 4.0, rho: 0.15 }, { df: 4.0, rho: 0.25 }, { df: 4.0, rho: 0.40 },
        { df: 10.0, rho: 0.25 }
      ].map do |cfg|
        els = detaches.map { |d| t_copula_equity_el(125, 0.015, 0.40, 5.0, cfg[:rho], cfg[:df], d) }
        bcs = calibrate_base_correlation(els, detaches).map { |r| r[:base_corr] }
        { df: cfg[:df], rho_market: cfg[:rho], base_corrs: bcs, skew_30_minus_3: bcs.last - bcs.first }
      end
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
      grid = (1..94).map { |i| 0.01 * i }
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
