# frozen_string_literal: true

module CdoTranche
  module Lhp
    module_function

    # CDF kerugian portofolio Large Homogeneous Pool (Vasicek):
    # P(L <= x) = Phi( (sqrt(1-rho)*Phi^-1(x/(1-R)) - Phi^-1(p)) / sqrt(rho) ).
    def loss_cdf(x, p_uncond, recovery, rho)
      lgd = 1.0 - recovery
      return 0.0 if x <= 0.0
      return 1.0 if x >= lgd

      c_p = Special.norm_inv(p_uncond.clamp(1.0e-12, 1.0 - 1.0e-12))
      c_x = Special.norm_inv((x / lgd).clamp(1.0e-12, 1.0 - 1.0e-12))
      Special.norm_cdf(((Math.sqrt(1.0 - rho) * c_x) - c_p) / Math.sqrt(rho))
    end

    # Ekspektasi kerugian equity tranche E[min(L, D)] dalam LHP = (1-R) * Phi_2(c_p, -Phi^-1(D/(1-R))...)?
    # Lebih tepat: int_0^D (1 - F_LHP(x)) dx via Gauss-Legendre 64-titik (akurat < 1e-13).
    def equity_el(detach, p_uncond, recovery, rho, n_gl: 96)
      lgd = 1.0 - recovery
      d_eff = [detach, lgd].min
      return 0.0 if d_eff <= 0.0

      nodes, weights = Special.gauss_legendre(n_gl)
      half = 0.5 * d_eff
      sum = 0.0
      nodes.each_with_index do |xi, idx|
        x = half * (xi + 1.0)
        sum += weights[idx] * (1.0 - loss_cdf(x, p_uncond, recovery, rho))
      end
      half * sum
    end

    def tranche_el(attach, detach, p_uncond, recovery, rho)
      equity_el(detach, p_uncond, recovery, rho) - equity_el(attach, p_uncond, recovery, rho)
    end

    # Anchor tertutup LHP via ortan bivariat Gaussian: E[max(L - A, 0)]?
    # E[min(L, D)] = D * (1 - F_LHP(D)) + (1 - R) * Phi_2(Phi^-1(p), -z_D, -sqrt(1 - rho)).
    def equity_el_bivariate(detach, p_uncond, recovery, rho)
      return 0.0 if detach <= 0.0

      lgd = 1.0 - recovery
      return p_uncond * lgd if detach >= lgd

      c_p = Special.norm_inv(p_uncond)
      z_d = Special.norm_inv(detach / lgd)
      # P(L > D) = P(M < (c_p - sqrt(1-rho)*z_D)/sqrt(rho))
      m_thresh = (c_p - (Math.sqrt(1.0 - rho) * z_d)) / Math.sqrt(rho)
      p_gt_d = Special.norm_cdf(m_thresh)
      # E[L * 1_{L <= D}] = (1-R) * P(Z_i <= c_p, M >= m_thresh) = (1-R) * Phi_2(c_p, -m_thresh, -sqrt(rho))
      e_l_le_d = lgd * Special.bivariate_norm_cdf(c_p, -m_thresh, -Math.sqrt(rho))
      (detach * p_gt_d) + e_l_le_d
    end
  end
end
