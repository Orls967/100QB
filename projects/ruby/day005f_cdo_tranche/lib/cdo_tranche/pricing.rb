# frozen_string_literal: true

module CdoTranche
  module Pricing
    module_function

    def tranche_pv_from_el(el_curve, attach, detach, times: Models::PAYMENT_TIMES, rate: Models::RISK_FREE_RATE)
      width = detach - attach
      pv_prot = 0.0
      pv_prem_unit = 0.0
      prev_t = 0.0
      prev_el = 0.0

      times.each_with_index do |t, idx|
        dt = t - prev_t
        df = Math.exp(-rate * t)
        el_t = el_curve[idx]
        pv_prot += df * (el_t - prev_el)
        avg_rem = width - (0.5 * (el_t + prev_el))
        pv_prem_unit += dt * df * [avg_rem, 0.0].max
        prev_t = t
        prev_el = el_t
      end

      fair_bps = pv_prem_unit > 1.0e-15 ? 10_000.0 * pv_prot / pv_prem_unit : 0.0
      { attach: attach, detach: detach, el_terminal: el_curve.last,
        pv_prot: pv_prot, pv_prem_unit: pv_prem_unit, fair_spread_bps: fair_bps }
    end

    def lhp_tranche_pv(hazard, recovery, rho, attach, detach)
      el_curve = Models::PAYMENT_TIMES.map do |t|
        p = Models.default_prob(hazard, t)
        Lhp.equity_el_bivariate(detach, p, recovery, rho) - Lhp.equity_el_bivariate(attach, p, recovery, rho)
      end
      tranche_pv_from_el(el_curve, attach, detach)
    end

    def asb_tranche_pv(obligors, rho, attach, detach, loss_unit: 0.0004)
      el_curve = Models::PAYMENT_TIMES.map do |t|
        Asb.tranche_el_asb(obligors, t, rho, attach, detach, loss_unit: loss_unit, n_gl: 48)
      end
      tranche_pv_from_el(el_curve, attach, detach)
    end

    def homogeneous_exact_pv(n_names, hazard, recovery, rho, attach, detach)
      el_curve = Models::PAYMENT_TIMES.map do |t|
        Asb.homogeneous_tranche_el_exact(n_names, hazard, recovery, t, rho, attach, detach, n_gl: 48)
      end
      tranche_pv_from_el(el_curve, attach, detach)
    end

    def granularity_study(ns: [25, 50, 125, 500], hazard: 0.015, recovery: 0.40, rho: 0.25, t: 5.0)
      p = Models.default_prob(hazard, t)
      lhp_senior = Lhp.equity_el_bivariate(0.10, p, recovery, rho) - Lhp.equity_el_bivariate(0.07, p, recovery, rho)
      ns.map do |n|
        ex = Asb.homogeneous_tranche_el_exact(n, hazard, recovery, t, rho, 0.07, 0.10, n_gl: 48)
        { n: n, exact_el: ex, lhp_el: lhp_senior, abs_err: (ex - lhp_senior).abs }
      end
    end
  end
end
