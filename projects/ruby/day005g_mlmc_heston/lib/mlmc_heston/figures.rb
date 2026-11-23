# frozen_string_literal: true

module MlmcHeston
  module Figures
    module_function

    def write_all(payload, fig_dir)
      {
        var_decay: write_var_decay(payload[:part2_rates], fig_dir),
        mean_decay: write_mean_decay(payload[:part2_rates], fig_dir),
        complexity: write_complexity(payload[:part3_complexity], fig_dir),
        samples: write_samples(payload[:part3_complexity][:heston_european][:rows].last, fig_dir),
        antithetic: write_antithetic(payload[:part4_antithetic], fig_dir)
      }
    end

    def write_var_decay(p2, fig_dir)
      ls = p2[:gbm_euler][:rows].map { |r| r[:level].to_f }
      out = File.join(fig_dir, 'day005g_level_variance_decay.svg')
      Svg.write_line_chart(out, title: 'Peluruhan Variansi Level MLMC (log2 V_l vs Level l)',
                                x_label: 'Level l', y_label: 'log2(V_l)', series: [
                                  { label: 'GBM Euler', x: ls, y: p2[:gbm_euler][:rows].map { |r| r[:log2_var] },
                                    color: '#2563eb', dashed: true },
                                  { label: 'GBM Milstein', x: ls, y: p2[:gbm_milstein][:rows].map { |r| r[:log2_var] },
                                    color: '#059669', dashed: false },
                                  { label: 'Heston FT', x: ls,
                                    y: p2[:heston_feller_violated][:rows].map { |r| r[:log2_var] },
                                    color: '#dc2626', dashed: false }
                                ])
      out
    end

    def write_mean_decay(p2, fig_dir)
      ls = p2[:gbm_euler][:rows].map { |r| r[:level].to_f }
      out = File.join(fig_dir, 'day005g_level_mean_decay.svg')
      Svg.write_line_chart(out, title: 'Peluruhan Bias Lemah MLMC (log2 |E[Y_l]| vs Level l)',
                                x_label: 'Level l', y_label: 'log2(|E[Y_l]|)', series: [
                                  { label: 'GBM Milstein', x: ls, y: p2[:gbm_milstein][:rows].map { |r| r[:log2_mean] },
                                    color: '#059669', dashed: false },
                                  { label: 'Heston FT', x: ls,
                                    y: p2[:heston_feller_violated][:rows].map { |r| r[:log2_mean] },
                                    color: '#dc2626', dashed: false }
                                ])
      out
    end

    def write_complexity(p3, fig_dir)
      xs = p3[:heston_european][:rows].map { |r| Math.log10(1.0 / r[:eps]) }
      out = File.join(fig_dir, 'day005g_complexity_scaling.svg')
      Svg.write_line_chart(out, title: 'Kompleksitas Biaya MLMC vs MC Standar (log10 Cost vs log10(1/eps))',
                                x_label: 'log10(1/eps)', y_label: 'log10(Total Steps Cost)', series: [
                                  { label: 'Std MC Heston', x: xs,
                                    y: p3[:heston_european][:rows].map { |r| Math.log10(r[:std_mc_cost]) },
                                    color: '#111827', dashed: true },
                                  { label: 'MLMC Heston', x: xs,
                                    y: p3[:heston_european][:rows].map { |r| Math.log10(r[:total_cost]) },
                                    color: '#dc2626', dashed: false },
                                  { label: 'MLMC GBM Milstein', x: xs,
                                    y: p3[:gbm_milstein][:rows].map { |r| Math.log10(r[:total_cost]) },
                                    color: '#059669', dashed: false }
                                ])
      out
    end

    def write_samples(last_row, fig_dir)
      ls = (0..last_row[:max_level]).map(&:to_f)
      ys = last_row[:n_samples].map { |n| Math.log10(n) }
      out = File.join(fig_dir, 'day005g_mlmc_samples_per_level.svg')
      Svg.write_line_chart(out, title: 'Alokasi Sampel Optimal N_l per Level (Heston eps=0.01)',
                                x_label: 'Level l', y_label: 'log10(N_l)', series: [
                                  { label: 'log10(N_l)', x: ls, y: ys, color: '#2563eb', dashed: false }
                                ])
      out
    end

    def write_antithetic(p4, fig_dir)
      ls = p4.map { |r| r[:level].to_f }
      out = File.join(fig_dir, 'day005g_antithetic_variance.svg')
      Svg.write_line_chart(out, title: 'Variansi Level Heston: Standar vs Antithetic Giles-Szpruch',
                                x_label: 'Level l', y_label: 'log10(V_l)', series: [
                                  { label: 'Standard', x: ls, y: p4.map { |r| Math.log10(r[:var_standard]) },
                                    color: '#dc2626', dashed: true },
                                  { label: 'Antithetic', x: ls, y: p4.map { |r| Math.log10(r[:var_antithetic]) },
                                    color: '#059669', dashed: false }
                                ])
      out
    end
  end
end
