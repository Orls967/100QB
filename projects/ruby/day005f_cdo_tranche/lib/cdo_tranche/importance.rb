# frozen_string_literal: true

module CdoTranche
  module Importance
    module_function

    # Pergeseran faktor optimal mu* untuk tranche [attach, detach] pada portofolio homogen:
    # P(L > attach | M = mu*) ~ 0.5 -> mu* = (c_p - sqrt(1-rho)*Phi^-1(attach/lgd)) / sqrt(rho).
    def optimal_factor_shift(p_uncond, recovery, rho, attach)
      lgd = 1.0 - recovery
      c_p = Special.norm_inv(p_uncond)
      q_target = (attach / lgd).clamp(1.0e-6, 1.0 - 1.0e-6)
      z_q = Special.norm_inv(q_target)
      (c_p - (Math.sqrt(1.0 - rho) * z_q)) / Math.sqrt(rho)
    end

    # Conditional MC + Exponential Twisting (Glasserman-Li 2005) untuk tranche [attach, detach].
    def compare_naive_vs_is(n_names: 125, hazard: 0.015, recovery: 0.40, t: 5.0, rho: 0.25,
                            attach: 0.15, detach: 0.30, n_paths: 100_000, seed: MAIN_SEED, workers: nil)
      exact = Asb.homogeneous_tranche_el_exact(n_names, hazard, recovery, t, rho, attach, detach, n_gl: 64)
      p_uncond = Models.default_prob(hazard, t)
      mu_star = optimal_factor_shift(p_uncond, recovery, rho, attach)
      c_p = Special.norm_inv(p_uncond)
      lgd_unit = (1.0 - recovery) / n_names
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)

      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 5_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 400 + (attach * 100).round, cidx))
        simulate_is_chunk(len, n_names, lgd_unit, c_p, sr, s1r, attach, detach, mu_star, rng)
      end

      summarize_is(chunks, n_paths, exact, mu_star, attach, detach)
    end

    def simulate_is_chunk(len, n_names, lgd_unit, c_p, sr, s1r, attach, detach, mu_star, rng)
      ns1 = 0.0
      ns2 = 0.0
      is1 = 0.0
      is2 = 0.0
      width = detach - attach
      i = 0
      while i < len
        # Naive MC (faktor M ~ N(0,1))
        m_n = rng.next_normal
        p_n = Special.norm_cdf((c_p - (sr * m_n)) / s1r)
        k_n = draw_binomial(n_names, p_n, rng)
        pay_n = ((k_n * lgd_unit) - attach).clamp(0.0, width)
        ns1 += pay_n
        ns2 += pay_n * pay_n

        # Two-stage IS: M ~ N(mu_star, 1) + exponential twisting per obligor jika E[L|M] < attach
        m_is = mu_star + rng.next_normal
        lr_m = Math.exp((-mu_star * m_is) + (0.5 * mu_star * mu_star))
        p_is = Special.norm_cdf((c_p - (sr * m_is)) / s1r)
        pay_is = twisted_tranche_payoff(n_names, lgd_unit, p_is, attach, width, rng) * lr_m
        is1 += pay_is
        is2 += pay_is * pay_is
        i += 1
      end
      [ns1, ns2, is1, is2]
    end
    private_class_method :simulate_is_chunk

    def draw_binomial(n_names, prob, rng)
      k = 0
      i = 0
      while i < n_names
        k += 1 if rng.next_uniform < prob
        i += 1
      end
      k
    end
    private_class_method :draw_binomial

    def twisted_tranche_payoff(n_names, lgd_unit, p_cond, attach, width, rng)
      q_target = attach / (n_names * lgd_unit)
      if p_cond >= q_target || p_cond <= 1.0e-12
        k = draw_binomial(n_names, p_cond, rng)
        return ((k * lgd_unit) - attach).clamp(0.0, width)
      end

      q_twist = q_target.clamp(1.0e-9, 1.0 - 1.0e-9)
      theta_l = Math.log((q_twist * (1.0 - p_cond)) / (p_cond * (1.0 - q_twist)))
      k = draw_binomial(n_names, q_twist, rng)
      psi = n_names * Math.log(1.0 - p_cond + (p_cond * Math.exp(theta_l)))
      lr_twist = Math.exp((-theta_l * k) + psi)
      ((k * lgd_unit) - attach).clamp(0.0, width) * lr_twist
    end
    private_class_method :twisted_tranche_payoff

    def summarize_is(chunks, n_paths, exact, mu_star, attach, detach)
      ns1 = chunks.sum(0.0) { |c| c[0] }
      ns2 = chunks.sum(0.0) { |c| c[1] }
      is1 = chunks.sum(0.0) { |c| c[2] }
      is2 = chunks.sum(0.0) { |c| c[3] }

      n_mean = ns1 / n_paths
      n_var = [(ns2 / n_paths) - (n_mean * n_mean), 0.0].max
      i_mean = is1 / n_paths
      i_var = [(is2 / n_paths) - (i_mean * i_mean), 0.0].max
      i_se = Math.sqrt(i_var / n_paths)
      {
        attach: attach, detach: detach, exact_el: exact, mu_star: mu_star,
        naive_el: n_mean, naive_se: Math.sqrt(n_var / n_paths), naive_var: n_var,
        is_el: i_mean, is_se: i_se, is_var: i_var,
        is_z: (i_mean - exact) / [i_se, 1.0e-18].max,
        variance_reduction: i_var.positive? ? n_var / i_var : Float::INFINITY
      }
    end
    private_class_method :summarize_is
  end
end
