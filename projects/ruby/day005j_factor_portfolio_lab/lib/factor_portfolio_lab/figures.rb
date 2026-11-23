# frozen_string_literal: true

module FactorPortfolioLab
  # Generates 26 SVG figures under outputs/figures/day005j_*.svg.
  module Figures
    module_function

    def generate_all(data, fig_dir)
      spiked_figures(data[:spiked], fig_dir) +
        grs_fm_rp_figures(data, fig_dir) +
        bl_mv_figures(data, fig_dir) +
        michaud_bt_figures(data, fig_dir) +
        part9_12_figures(data, fig_dir)
    end

    def spiked_figures(sp_data, dir)
      spikes = sp_data[:levels].map { |l| l[:spike] }
      p1 = File.join(dir, 'day005j_01_bbp_eigenvalue_phase.svg')
      p2 = File.join(dir, 'day005j_02_bbp_eigenvector_alignment.svg')
      p3 = File.join(dir, 'day005j_03_pca_effective_rank.svg')

      Svg.line_chart(
        path: p1, title: 'BBP Spiked Eigenvalue: Empirical vs Theory', x_label: 'Spike (ell)', y_label: 'Lambda 1',
        series: [
          { name: 'Theory BBP', x: spikes, y: sp_data[:levels].map { |l| l[:theory_lambda1] }, dashed: true },
          { name: 'Empirical MC', x: spikes, y: sp_data[:levels].map { |l| l[:emp_lambda1_mean] } }
        ]
      )
      Svg.line_chart(
        path: p2, title: 'BBP Eigenvector Alignment |<v1, u1>|^2', x_label: 'Spike (ell)', y_label: 'Alignment^2',
        series: [
          { name: 'Theory BBP', x: spikes, y: sp_data[:levels].map { |l| l[:theory_alignment_sq] }, dashed: true },
          { name: 'Empirical MC', x: spikes, y: sp_data[:levels].map { |l| l[:emp_alignment_sq_mean] } }
        ]
      )
      Svg.line_chart(
        path: p3, title: 'PCA Effective Rank vs Spike Strength', x_label: 'Spike (ell)', y_label: 'Effective Rank',
        series: [{ name: 'Eff Rank', x: spikes, y: sp_data[:levels].map { |l| l[:eff_rank_mean] } }]
      )
      [p1, p2, p3]
    end

    def grs_fm_rp_figures(data, dir)
      grs = data[:grs]
      fm = data[:fama_macbeth]
      rp = data[:risk_parity]
      paths = [
        'day005j_04_grs_vs_wald_rejection.svg', 'day005j_05_wald_overrejection.svg',
        'day005j_06_fm_vs_shanken_coverage.svg', 'day005j_07_fm_naive_coverage.svg',
        'day005j_08_rp_solver_iterations.svg', 'day005j_09_rp_vs_eq_vol.svg'
      ].map { |f| File.join(dir, f) }

      Svg.bar_chart(
        path: paths[0], title: 'GRS Exact F Rejection Rate (alpha=5%)', y_label: 'Rejection',
        labels: grs[:regimes].map { |r| r[:label] }, values: grs[:regimes].map { |r| r[:grs_rejection_05] }
      )
      Svg.bar_chart(
        path: paths[1], title: 'Asymptotic Wald Rejection Rate (alpha=5%)', y_label: 'Rejection',
        labels: grs[:regimes].map { |r| r[:label] },
        values: grs[:regimes].map { |r| r[:wald_rejection_05] }, color: '#d62728'
      )
      Svg.bar_chart(
        path: paths[2], title: 'Shanken 95% CI Coverage by Parameter', y_label: 'Coverage',
        labels: %w[Intercept Factor1 Factor2 Factor3], values: fm[:coverage_shanken_95], color: '#2ca02c'
      )
      Svg.bar_chart(
        path: paths[3], title: 'Naive Fama-MacBeth 95% CI Coverage', y_label: 'Coverage',
        labels: %w[Intercept Factor1 Factor2 Factor3], values: fm[:coverage_fm_95], color: '#ff7f0e'
      )
      Svg.line_chart(
        path: paths[4], title: 'Risk Parity Iterations: Newton vs CCD', x_label: 'N Assets', y_label: 'Iterations',
        series: [
          { name: 'Spinu Newton', x: rp[:comparisons].map { |c| c[:n_assets] },
            y: rp[:comparisons].map { |c| c[:newton_iters] } },
          { name: 'Cyclic CD', x: rp[:comparisons].map { |c| c[:n_assets] },
            y: rp[:comparisons].map { |c| c[:ccd_iters] } }
        ]
      )
      Svg.bar_chart(
        path: paths[5], title: 'Portfolio Volatility: Risk Parity vs Equal Weight', y_label: 'Volatility',
        labels: rp[:comparisons].flat_map { |c| ["RP_N#{c[:n_assets]}", "EW_N#{c[:n_assets]}"] },
        values: rp[:comparisons].flat_map { |c| [c[:rp_vol], c[:eq_vol]] }, color: '#9467bd'
      )
      paths
    end

    def bl_mv_figures(data, dir)
      bl_figures(data[:black_litterman], dir) + mv_figures(data[:markowitz], dir)
    end

    def bl_figures(bl, dir)
      p10 = File.join(dir, 'day005j_10_bl_idzorek_sweep.svg')
      p11 = File.join(dir, 'day005j_11_bl_exact_vs_mc.svg')
      Svg.line_chart(
        path: p10, title: 'Black-Litterman Idzorek View Shift vs Confidence',
        x_label: 'Confidence', y_label: '||mu_post - pi||',
        series: [{ name: 'Posterior Shift', x: bl[:idzorek_sweep].map { |r| r[:confidence] },
                   y: bl[:idzorek_sweep].map { |r| r[:shift_norm] } }]
      )
      Svg.bar_chart(
        path: p11, title: 'Black-Litterman Posterior Mean: Exact vs 500k MC', y_label: 'Posterior Mean',
        labels: (0...bl[:mu_exact].length).flat_map { |i| ["Ex#{i}", "MC#{i}"] },
        values: bl[:mu_exact].each_with_index.flat_map { |m, i| [m, bl[:mu_mc][i]] }
      )
      [p10, p11]
    end

    def mv_figures(mv, dir)
      p12 = File.join(dir, 'day005j_12_markowitz_frontiers.svg')
      p13 = File.join(dir, 'day005j_13_markowitz_kkt_residuals.svg')
      Svg.line_chart(
        path: p12, title: 'Efficient Frontiers: Unconstrained, Long-Only, Capped',
        x_label: 'Volatility', y_label: 'Expected Return',
        series: [
          { name: 'Unconstrained', x: mv[:unconstrained].map { |p| p[:volatility] },
            y: mv[:unconstrained].map { |p| p[:expected_return] }, dashed: true },
          { name: 'Long-Only', x: mv[:long_only].map { |p| p[:volatility] },
            y: mv[:long_only].map { |p| p[:expected_return] } },
          { name: 'Capped 12%', x: mv[:capped].map { |p| p[:volatility] },
            y: mv[:capped].map { |p| p[:expected_return] } }
        ]
      )
      Svg.line_chart(
        path: p13, title: 'QP Solver Max KKT Residual Across Frontier', x_label: 'Kappa', y_label: 'Max KKT',
        series: [
          { name: 'Long-Only', x: mv[:long_only].map { |p| p[:kappa] }, y: mv[:long_only].map { |p| p[:max_kkt] } },
          { name: 'Capped 12%', x: mv[:capped].map { |p| p[:kappa] }, y: mv[:capped].map { |p| p[:max_kkt] } }
        ]
      )
      [p12, p13]
    end

    def michaud_bt_figures(data, dir)
      res = data[:resampling]
      bt = data[:backtest]
      strats = bt[:summary].keys
      paths = [
        'day005j_14_michaud_l1_error.svg', 'day005j_15_michaud_hhi.svg',
        'day005j_16_backtest_sharpe_5bps.svg', 'day005j_17_backtest_sharpe_20bps.svg',
        'day005j_18_backtest_turnover.svg', 'day005j_19_backtest_max_drawdown.svg'
      ].map { |f| File.join(dir, f) }

      Svg.bar_chart(
        path: paths[0], title: 'Distance to True Optimal Weights (L1)', y_label: 'L1 Distance',
        labels: %w[PlugIn Shrinkage Michaud],
        values: [res[:plugin][:l1_to_true], res[:shrunk][:l1_to_true], res[:resampled][:l1_to_true]], color: '#2ca02c'
      )
      Svg.bar_chart(
        path: paths[1], title: 'Portfolio Concentration (HHI)', y_label: 'HHI',
        labels: %w[TrueOpt PlugIn Shrinkage Michaud],
        values: [res[:true_opt][:hhi], res[:plugin][:hhi], res[:shrunk][:hhi], res[:resampled][:hhi]]
      )
      Svg.bar_chart(
        path: paths[2], title: 'Rolling Two-Regime Sharpe (5 bps TC, R=300)', y_label: 'Sharpe',
        labels: strats, values: strats.map { |s| bt[:summary][s][5.0][:sharpe] }
      )
      Svg.bar_chart(
        path: paths[3], title: 'Rolling Two-Regime Sharpe (20 bps TC, R=300)', y_label: 'Sharpe',
        labels: strats, values: strats.map { |s| bt[:summary][s][20.0][:sharpe] }, color: '#ff7f0e'
      )
      Svg.bar_chart(
        path: paths[4], title: 'Mean Rebalance Turnover Across Strategies', y_label: 'Turnover',
        labels: strats, values: strats.map { |s| bt[:summary][s][5.0][:mean_turnover] }, color: '#d62728'
      )
      Svg.bar_chart(
        path: paths[5], title: 'Mean Max Drawdown Across Strategies (5 bps)', y_label: 'Max DD',
        labels: strats, values: strats.map { |s| bt[:summary][s][5.0][:max_drawdown] }, color: '#8c564b'
      )
      paths
    end

    def part9_12_figures(data, dir)
      paths = [
        'day005j_20_attribution_misspec_k.svg', 'day005j_21_stress_scenarios_rp_vs_mv.svg',
        'day005j_22_reverse_stress_mahalanobis.svg', 'day005j_23_robustness_grs_gaussian_vs_t4.svg',
        'day005j_24_robustness_cvar_t4.svg', 'day005j_25_complexity_scaling_ms.svg',
        'day005j_26_complexity_loglog_slopes.svg'
      ].map { |f| File.join(dir, f) }

      draw_attr_stress_figs(data[:attribution], data[:stress], paths)
      draw_rob_cmp_figs(data[:robustness], data[:complexity], paths)
      paths
    end

    def draw_attr_stress_figs(att, st, paths)
      Svg.bar_chart(
        path: paths[0], title: 'Systematic Share Under Misspecified K (True K=3)', y_label: 'Sys Share',
        labels: att[:k_models].map { |m| "K=#{m[:k_fit]}" },
        values: att[:k_models].map { |m| m[:systematic_share] }, color: '#17becf'
      )
      Svg.bar_chart(
        path: paths[1], title: 'Stress Scenario Total PnL: Risk Parity vs Markowitz', y_label: 'PnL',
        labels: st[:scenarios].keys.flat_map { |k| ["#{k[0, 6]}_RP", "#{k[0, 6]}_MV"] },
        values: st[:scenarios].values.flat_map { |v| [v[:rp][:total_pnl], v[:mv][:total_pnl]] }, color: '#d62728'
      )
      Svg.bar_chart(
        path: paths[2], title: 'Reverse Stress Min Mahalanobis Distance (-15% Loss)', y_label: 'Mahal Dist',
        labels: %w[RiskParity Markowitz],
        values: [st[:reverse_rp][:mahalanobis_dist], st[:reverse_mv][:mahalanobis_dist]], color: '#9467bd'
      )
    end

    def draw_rob_cmp_figs(rob, cmp, paths)
      Svg.bar_chart(
        path: paths[3], title: 'GRS & Wald Rejection Rate: Gaussian vs Student-t(4)', y_label: 'Rejection',
        labels: %w[GRS_Gauss Wald_Gauss GRS_t4 Wald_t4],
        values: [rob[:grs_gaussian_rej][:grs_rej_05], rob[:grs_gaussian_rej][:wald_rej_05],
                 rob[:grs_t4_rej][:grs_rej_05], rob[:grs_t4_rej][:wald_rej_05]],
        color: '#ff7f0e'
      )
      Svg.bar_chart(
        path: paths[4], title: 'Portfolio 95% CVaR Under Student-t(4) Shocks', y_label: 'CVaR 95%',
        labels: rob[:portfolio_t4].keys, values: rob[:portfolio_t4].values.map { |v| v[:cvar_95] }, color: '#d62728'
      )
      ns = cmp[:rows].map { |r| r[:n_assets] }
      Svg.line_chart(
        path: paths[5], title: 'Runtime Scaling Across N Assets (ms)', x_label: 'N Assets', y_label: 'Time (ms)',
        series: [
          { name: 'Invert', x: ns, y: cmp[:rows].map { |r| r[:invert_ms] } },
          { name: 'RP Newton', x: ns, y: cmp[:rows].map { |r| r[:rp_newton_ms] } },
          { name: 'QP Markowitz', x: ns, y: cmp[:rows].map { |r| r[:qp_ms] } },
          { name: 'LW Shrink', x: ns, y: cmp[:rows].map { |r| r[:lw_ms] } }
        ]
      )
      Svg.bar_chart(
        path: paths[6], title: 'Empirical Log-Log Scaling Exponents', y_label: 'Exponent',
        labels: %w[Invert RP_Newton QP_Solve LW_Shrink],
        values: [cmp[:slope_invert], cmp[:slope_rp_newton], cmp[:slope_qp], cmp[:slope_lw]], color: '#2ca02c'
      )
    end
  end
end
