# frozen_string_literal: true

module CdoTranche
  module Figures
    module_function

    def write_all(payload, fig_dir)
      {
        granularity: write_granularity(payload[:part1_anchors], fig_dir),
        spreads: write_spreads(payload[:part2_valuation], fig_dir),
        base_corr: write_base_corr(payload[:part3_base_corr], fig_dir),
        vr_is: write_vr_is(payload[:part4_importance], fig_dir),
        model_risk: write_model_risk(payload[:part5_model_risk], fig_dir)
      }
    end

    def write_granularity(p1, fig_dir)
      xs = p1[:granularity].map { |r| Math.log10(r[:n]) }
      ys = p1[:granularity].map { |r| Math.log10(r[:abs_err]) }
      out = File.join(fig_dir, 'day005f_lhp_granularity.svg')
      Svg.write_line_chart(out, title: 'Konvergensi Galat LHP vs ASB Berhingga (log10 |Err| vs log10 N)',
                                x_label: 'log10(N)', y_label: 'log10(|EL_ASB - EL_LHP|)', series: [
                                  { label: 'Galat LHP 7-10%', x: xs, y: ys, color: '#2563eb', dashed: false }
                                ])
      out
    end

    def write_spreads(p2, fig_dir)
      xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
      out = File.join(fig_dir, 'day005f_tranche_spreads.svg')
      Svg.write_line_chart(out, title: 'Fair Spread Tranche CDO 5Y (bps): Homogen vs Heterogen',
                                x_label: 'Indeks Tranche (1..6)', y_label: 'Fair Spread (bps)', series: [
                                  { label: 'Homogen', x: xs, y: p2[:homogeneous].map { |r| r[:fair_spread_bps] },
                                    color: '#2563eb', dashed: false },
                                  { label: 'Heterogen', x: xs, y: p2[:heterogeneous].map { |r| r[:fair_spread_bps] },
                                    color: '#dc2626', dashed: true }
                                ])
      out
    end

    def write_base_corr(p3, fig_dir)
      rows = p3[:base_correlations]
      out = File.join(fig_dir, 'day005f_base_correlation_skew.svg')
      Svg.write_line_chart(out, title: 'Base Correlation Skew (Kalibrasi ke Pasar Student-t df=4)',
                                x_label: 'Detachment D', y_label: 'Base Correlation', series: [
                                  { label: 'Base Corr rho(D)', x: rows.map { |r| r[:detach] },
                                    y: rows.map { |r| r[:base_corr] }, color: '#059669', dashed: false }
                                ])
      out
    end

    def write_vr_is(p4, fig_dir)
      xs = p4.map { |r| r[:attach] * 100.0 }
      ys = p4.map { |r| Math.log10(r[:variance_reduction]) }
      out = File.join(fig_dir, 'day005f_is_variance_reduction.svg')
      Svg.write_line_chart(out, title: 'Log10 Penghematan Variansi Exponential Twisting IS',
                                x_label: 'Attachment (%)', y_label: 'log10(VR Factor)', series: [
                                  { label: 'log10(VR)', x: xs, y: ys, color: '#d97706', dashed: false }
                                ])
      out
    end

    def write_model_risk(p5, fig_dir)
      xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
      out = File.join(fig_dir, 'day005f_model_risk_spreads.svg')
      Svg.write_line_chart(out, title: 'Model Risk Fair Spread (bps): Gaussian vs t(4) vs Random Recovery',
                                x_label: 'Indeks Tranche (1..6)', y_label: 'Fair Spread (bps)', series: [
                                  { label: 'Gaussian', x: xs, y: p5[:gaussian].map { |r| r[:fair_spread_bps] },
                                    color: '#2563eb', dashed: false },
                                  { label: 'Student-t(4)', x: xs, y: p5[:t4].map { |r| r[:fair_spread_bps] },
                                    color: '#dc2626', dashed: false },
                                  { label: 'Random Recovery', x: xs,
                                    y: p5[:random_recovery].map { |r| r[:fair_spread_bps] },
                                    color: '#059669', dashed: true }
                                ])
      out
    end
  end
end
