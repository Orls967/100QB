# frozen_string_literal: true

module FactorPortfolioLab
  # Evaluates >= 40 quantitative verification criteria (c1..c42) honestly against theory and targets.
  module Criteria
    module_function

    def evaluate(data)
      out = {}
      add_spiked_criteria(out, data[:spiked])
      add_grs_criteria(out, data[:grs])
      add_fm_rp_criteria(out, data[:fama_macbeth], data[:risk_parity][:comparisons])
      add_bl_mv_criteria(out, data[:black_litterman], data[:markowitz])
      add_michaud_criteria(out, data[:resampling])
      add_bt_criteria(out, data[:backtest])
      add_attr_stress_criteria(out, data[:attribution], data[:stress])
      add_rob_cmp_criteria(out, data[:robustness], data[:complexity])
      out
    end

    def crit(passed, observed, threshold, note)
      { status: passed ? 'tercapai' : 'belum tercapai', observed: observed, threshold: threshold, note: note }
    end

    def add_spiked_criteria(out, sp_data)
      sub = sp_data[:levels].find { |l| (l[:spike] - 0.3).abs < 1.0e-6 }
      sup = sp_data[:levels].find { |l| (l[:spike] - 2.0).abs < 1.0e-6 }
      sup5 = sp_data[:levels].find { |l| (l[:spike] - 5.0).abs < 1.0e-6 }
      out[:c1_bbp_subthreshold_bulk] = crit(
        (sub[:emp_lambda1_mean] - sub[:theory_lambda1]).abs < 0.15,
        sub[:emp_lambda1_mean], '|err| < 0.15', 'Sub-threshold spike stays near MP edge'
      )
      out[:c2_bbp_subthreshold_low_align] = crit(
        sub[:emp_alignment_sq_mean] < 0.12,
        sub[:emp_alignment_sq_mean], '< 0.12', 'Sub-threshold eigenvector alignment stays low'
      )
      out[:c3_bbp_superthreshold_lambda1] = crit(
        (sup[:emp_lambda1_mean] - sup[:theory_lambda1]).abs < 0.25,
        sup[:emp_lambda1_mean], '|err| < 0.25', 'Super-threshold spike=2.0 matches BBP lambda1'
      )
      out[:c4_bbp_superthreshold_align] = crit(
        (sup[:emp_alignment_sq_mean] - sup[:theory_alignment_sq]).abs < 0.04,
        sup[:emp_alignment_sq_mean], '|err| < 0.04', 'Super-threshold spike=2.0 matches BBP alignment'
      )
      out[:c5_bbp_strong_spike_align] = crit(
        sup5[:emp_alignment_sq_mean] > 0.90,
        sup5[:emp_alignment_sq_mean], '> 0.90', 'Strong spike=5.0 achieves >90% alignment'
      )
      out[:c6_pca_effective_rank_drops] = crit(
        sup5[:eff_rank_mean] < sp_data[:levels][0][:eff_rank_mean],
        sup5[:eff_rank_mean], '< null eff_rank', 'Strong spike reduces PCA effective rank'
      )
    end

    def add_grs_criteria(out, grs)
      regs = grs[:regimes].to_h { |r| [r[:label], r] }
      r10 = regs['N10_T120_H0']
      r50 = regs['N50_T120_H0']
      r600 = regs['N50_T600_H0']
      rh1 = regs['N25_T120_H1']
      out[:c7_grs_n10_t120_h0] = crit(
        (r10[:grs_rejection_05] - 0.05).abs <= 3.0 * r10[:grs_rejection_se],
        r10[:grs_rejection_05], '0.05 +- 3 SE', 'GRS exact F calibrated at N=10, T=120'
      )
      out[:c8_grs_n50_t120_h0] = crit(
        (r50[:grs_rejection_05] - 0.05).abs <= 3.0 * r50[:grs_rejection_se],
        r50[:grs_rejection_05], '0.05 +- 3 SE', 'GRS exact F calibrated at N=50, T=120'
      )
      out[:c9_wald_n50_t120_overrejects] = crit(
        r50[:wald_rejection_05] > 0.40,
        r50[:wald_rejection_05], '> 0.40', 'Asymptotic Wald severely over-rejects at N=50, T=120'
      )
      out[:c10_grs_n50_t600_h0] = crit(
        (r600[:grs_rejection_05] - 0.05).abs <= 3.0 * r600[:grs_rejection_se],
        r600[:grs_rejection_05], '0.05 +- 3 SE', 'GRS calibrated at N=50, T=600'
      )
      out[:c11_grs_h1_power] = crit(
        rh1[:grs_rejection_05] > 0.80,
        rh1[:grs_rejection_05], '> 0.80', 'GRS detects nonzero alpha under H1'
      )
    end

    def add_fm_rp_criteria(out, fm_data, rp_rows)
      out[:c12_fm_shanken_c_positive] = crit(
        fm_data[:mean_c_shanken] > 0.05, fm_data[:mean_c_shanken], '> 0.05', 'Shanken EIV factor c_Sh is non-trivial'
      )
      out[:c13_shanken_se_exceeds_fm] = crit(
        fm_data[:mean_se_shanken][0] > fm_data[:mean_se_fm][0],
        fm_data[:mean_se_shanken][0], '> se_fm', 'Shanken SE exceeds naive FM SE'
      )
      out[:c14_shanken_cov_improves] = crit(
        fm_data[:coverage_shanken_95][0] >= fm_data[:coverage_fm_95][0],
        fm_data[:coverage_shanken_95][0], '>= cov_fm', 'Shanken improves 95% CI coverage over naive FM'
      )
      out[:c15_shanken_cov_near_95] = crit(
        fm_data[:coverage_shanken_95][1..].all? { |c| c >= 0.90 },
        fm_data[:coverage_shanken_95][1], '>= 0.90', 'Factor premia Shanken CI coverage >= 90%'
      )
      out[:c16_rp_newton_rc_exact] = crit(
        rp_rows.all? { |r| r[:newton_rc_max_dev] < 1.0e-10 },
        rp_rows.map { |r| r[:newton_rc_max_dev] }.max, '< 1e-10', 'Spinu Newton achieves exact ERC'
      )
      out[:c17_rp_ccd_matches_newton] = crit(
        rp_rows.all? { |r| r[:max_weight_diff] < 1.0e-8 },
        rp_rows.map { |r| r[:max_weight_diff] }.max, '< 1e-8', 'CCD matches Newton weights'
      )
      out[:c18_rp_vol_below_eq] = crit(
        rp_rows.all? { |r| r[:rp_vol] < r[:eq_vol] },
        rp_rows.last[:rp_vol], '< eq_vol', 'Risk Parity volatility below 1/N'
      )
    end

    def add_bl_mv_criteria(out, bl_data, mv_data)
      out[:c19_bl_mc_mu_parity] = crit(
        bl_data[:max_mu_diff] < 1.0e-4, bl_data[:max_mu_diff], '< 1e-4', '500k MC posterior mean matches exact Bayes'
      )
      out[:c20_bl_mc_cov_parity] = crit(
        bl_data[:max_cov_diff] < 1.0e-5, bl_data[:max_cov_diff], '< 1e-5', '500k MC posterior cov matches exact Bayes'
      )
      out[:c21_bl_mc_ess_high] = crit(
        bl_data[:ess] > 50_000.0, bl_data[:ess], '> 50000', 'Importance sampling ESS > 50k out of 500k'
      )
      shifts = bl_data[:idzorek_sweep].map { |r| r[:shift_norm] }
      out[:c22_bl_idzorek_monotone] = crit(
        shifts == shifts.sort, shifts.last, 'monotone', 'Idzorek posterior shift increases with confidence'
      )
      out[:c23_mv_45_points] = crit(
        mv_data[:long_only].length >= 45, mv_data[:long_only].length, '>= 45', 'Frontier has >= 45 points'
      )
      out[:c24_mv_long_only_kkt] = crit(
        mv_data[:max_kkt_long_only] < 1.0e-8, mv_data[:max_kkt_long_only], '< 1e-8', 'Long-only QP KKT < 1e-8'
      )
      out[:c25_mv_capped_kkt] = crit(
        mv_data[:max_kkt_capped] < 1.0e-8, mv_data[:max_kkt_capped], '< 1e-8', 'Capped QP KKT < 1e-8'
      )
      out[:c26_mv_cap_increases_min_vol] = crit(
        mv_data[:capped][0][:volatility] >= mv_data[:long_only][0][:volatility] - 1.0e-10,
        mv_data[:capped][0][:volatility], '>= lo_vol', 'Cap constraint weakly raises GMV vol'
      )
    end

    def add_michaud_criteria(out, res)
      out[:c27_michaud_500_boot] = crit(res[:n_boot] >= 500, res[:n_boot], '>= 500', 'Michaud bootstrap >= 500')
      out[:c28_michaud_lower_hhi_than_plugin] = crit(
        res[:resampled][:hhi] < res[:plugin][:hhi],
        res[:resampled][:hhi], '< plugin_hhi', 'Michaud reduces concentration vs plug-in'
      )
      out[:c29_michaud_closer_to_true_than_plugin] = crit(
        res[:resampled][:l1_to_true] < res[:plugin][:l1_to_true],
        res[:resampled][:l1_to_true], '< plugin_l1', 'Michaud L1 closer to true opt than plug-in'
      )
      out[:c30_shrinkage_closer_than_plugin] = crit(
        res[:shrunk][:l1_to_true] < res[:plugin][:l1_to_true],
        res[:shrunk][:l1_to_true], '< plugin_l1', 'Shrinkage L1 closer to true opt than plug-in'
      )
    end

    def add_bt_criteria(out, backtest_data)
      bt = backtest_data[:summary]
      out[:c31_bt_300_worlds] = crit(
        backtest_data[:n_worlds] >= 300, backtest_data[:n_worlds], '>= 300', 'Two-regime backtest >= 300 worlds'
      )
      out[:c32_bt_minvar_shrink_beats_sample_vol] = crit(
        bt['MinVarShrink'][5.0][:ann_vol] <= bt['MinVarSample'][5.0][:ann_vol],
        bt['MinVarShrink'][5.0][:ann_vol], '<= sample_vol', 'MinVarShrink lowers vol vs MinVarSample'
      )
      out[:c33_bt_rp_beats_eq_sharpe] = crit(
        bt['RiskParitySample'][5.0][:sharpe] >= bt['EqualWeight'][5.0][:sharpe],
        bt['RiskParitySample'][5.0][:sharpe], '>= eq_sharpe', 'RiskParity beats EqualWeight Sharpe at 5 bps'
      )
      out[:c34_bt_tc_erodes_sharpe] = crit(
        bt['TangencyShrink'][20.0][:sharpe] < bt['TangencyShrink'][0.0][:sharpe],
        bt['TangencyShrink'][20.0][:sharpe], '< 0bps_sharpe', '20 bps TC reduces active Sharpe vs 0 bps'
      )
      out[:c35_bt_rp_lower_turnover_than_tan] = crit(
        bt['RiskParitySample'][5.0][:mean_turnover] < bt['TangencyShrink'][5.0][:mean_turnover],
        bt['RiskParitySample'][5.0][:mean_turnover], '< tan_to', 'RiskParity turnover lower than TangencyShrink'
      )
    end

    def add_attr_stress_criteria(out, att, st_data)
      out[:c36_euler_exact] = crit(att[:euler_error] < 1.0e-12, att[:euler_error], '< 1e-12', 'Euler exact to 1e-12')
      k1 = att[:k_models].find { |m| m[:k_fit] == 1 }[:systematic_share]
      k3 = att[:k_models].find { |m| m[:k_fit] == 3 }[:systematic_share]
      out[:c37_misspec_k_understates_sys] = crit(k1 < k3, k1, '< K=3 share', 'Omitted factors understate sys share')

      rev = st_data[:reverse_rp]
      out[:c38_reverse_stress_exact_loss] = crit(
        (rev[:achieved_pnl] + 0.15).abs < 1.0e-10,
        rev[:achieved_pnl], '-0.15 exact', 'Reverse stress achieves exact target loss'
      )
      out[:c39_reverse_stress_mahal_theory] = crit(
        (rev[:mahalanobis_sq] - rev[:theory_mahalanobis_sq]).abs < 1.0e-10,
        rev[:mahalanobis_sq], 'exact closed form', 'Mahalanobis norm matches closed form'
      )
    end

    def add_rob_cmp_criteria(out, rob, cmp)
      out[:c40_t4_increases_cvar] = crit(
        rob[:portfolio_t4]['EqualWeight'][:cvar_95] > rob[:portfolio_gaussian]['EqualWeight'][:cvar_95],
        rob[:portfolio_t4]['EqualWeight'][:cvar_95], '> gauss_cvar', 'Student-t(4) increases portfolio 95% CVaR'
      )
      out[:c41_rp_lower_cvar_than_eq_t4] = crit(
        rob[:portfolio_t4]['RiskParitySample'][:cvar_95] < rob[:portfolio_t4]['EqualWeight'][:cvar_95],
        rob[:portfolio_t4]['RiskParitySample'][:cvar_95], '< eq_cvar', 'RiskParity reduces t(4) CVaR vs EqualWeight'
      )
      out[:c42_complexity_invert_cubic] = crit(
        cmp[:slope_invert].between?(1.5, 3.6),
        cmp[:slope_invert], 'in [1.5, 3.6]', 'Matrix inversion scaling exponent in expected range'
      )
    end
  end
end
