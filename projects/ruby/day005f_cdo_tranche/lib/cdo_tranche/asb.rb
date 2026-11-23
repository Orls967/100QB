# frozen_string_literal: true

module CdoTranche
  module Asb
    module_function

    # Distribusi kerugian bersyarat via rekursi Andersen-Sidenius-Basu (DP).
    # Unit diskretisasi loss_unit = 0.0004 (1/2500) cukup halus untuk portofolio heterogen.
    def conditional_loss_distribution(obligors, probs, loss_unit: 0.0004)
      max_units = obligors.sum { |ob| (ob[:lgd] / loss_unit).round }
      dist = Array.new(max_units + 1, 0.0)
      dist[0] = 1.0
      cur_max = 0

      obligors.each_with_index do |ob, idx|
        u = (ob[:lgd] / loss_unit).round
        p = probs[idx]
        q = 1.0 - p
        k = cur_max
        while k >= 0
          val = dist[k]
          if val.positive?
            dist[k + u] += val * p
            dist[k] = val * q
          end
          k -= 1
        end
        cur_max += u
      end
      [dist, loss_unit]
    end

    # Distribusi binomial eksak untuk portofolio homogen.
    def homogeneous_binomial_el(n_names, lgd_per_name, p_cond, attach, detach)
      el = 0.0
      comb = 1.0
      q = 1.0 - p_cond
      log_p = p_cond > 1.0e-15 ? Math.log(p_cond) : -700.0
      log_q = q > 1.0e-15 ? Math.log(q) : -700.0
      (0..n_names).each do |k|
        comb *= (n_names - k + 1).to_f / k if k.positive?
        prob = Math.exp(Math.lgamma(n_names + 1).first - Math.lgamma(k + 1).first -
                        Math.lgamma(n_names - k + 1).first + (k * log_p) + ((n_names - k) * log_q))
        loss = k * lgd_per_name
        pay = (loss - attach).clamp(0.0, detach - attach)
        el += prob * pay
      end
      el
    end

    def conditional_prob_1f(p_uncond, rho, factor_m, copula: :gaussian, df: 4.0, scale_w: 1.0)
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)
      if copula == :gaussian
        thresh = Special.norm_inv(p_uncond.clamp(1.0e-12, 1.0 - 1.0e-12))
        Special.norm_cdf((thresh - (sr * factor_m)) / s1r)
      else
        thresh = Special.t_inv(p_uncond.clamp(1.0e-10, 1.0 - 1.0e-10), df)
        Special.norm_cdf(((thresh / scale_w) - (sr * factor_m)) / s1r)
      end
    end

    def tranche_el_asb(obligors, t, rho, attach, detach, loss_unit: 0.0004, n_gl: 64)
      nodes, weights = Special.gauss_legendre(n_gl)
      half = 7.0
      uncond_probs = obligors.map { |ob| Models.default_prob(ob[:hazard], t) }
      total_el = 0.0

      nodes.each_with_index do |xi, idx|
        m = half * xi
        cond_probs = uncond_probs.map { |p| conditional_prob_1f(p, rho, m) }
        dist, unit = conditional_loss_distribution(obligors, cond_probs, loss_unit: loss_unit)
        cond_el = 0.0
        dist.each_with_index do |pr, k|
          next if pr <= 0.0

          loss = k * unit
          pay = (loss - attach).clamp(0.0, detach - attach)
          cond_el += pr * pay
        end
        total_el += weights[idx] * Special.norm_pdf(m) * cond_el
      end
      half * total_el
    end

    def homogeneous_tranche_el_exact(n_names, hazard, recovery, t, rho, attach, detach, n_gl: 64)
      nodes, weights = Special.gauss_legendre(n_gl)
      half = 7.0
      p_uncond = Models.default_prob(hazard, t)
      lgd_per_name = (1.0 - recovery) / n_names
      total_el = 0.0

      nodes.each_with_index do |xi, idx|
        m = half * xi
        p_c = conditional_prob_1f(p_uncond, rho, m)
        cond_el = homogeneous_binomial_el(n_names, lgd_per_name, p_c, attach, detach)
        total_el += weights[idx] * Special.norm_pdf(m) * cond_el
      end
      half * total_el
    end
  end
end
