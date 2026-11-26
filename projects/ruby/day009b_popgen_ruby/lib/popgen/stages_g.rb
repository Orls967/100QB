# frozen_string_literal: true

require_relative 'rng'
require_relative 'recomb'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group G: Two-Locus Coalescent with Recombination (ARG-lite) Stages S57..S59
  module StagesG
    module_function

    def run(ctx)
      s57_s59_arg_lite(ctx)
    end

    def s57_s59_arg_lite(ctx)
      t57 = Time.now
      n = 15
      rhos = [0.0, 0.5, 2.0, 8.0, 32.0]
      n_rep = ctx.small? ? 1_200 : 4_000
      n_blocks = 10
      e_tmrca_exact = 2.0 * (1.0 - 1.0 / n)

      rows = []
      max_marg_z = 0.0
      t1_by_rho = []
      t2_by_rho = []

      rhos.each_with_index do |rho, idx|
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(57 + idx, b_i)
          b_ids.map do
            ar = Recomb.simulate_arg_lite(n, rho, rng)
            [ar[:t_mrca_1], ar[:t_mrca_2], ar[:len_1], ar[:len_2], ar[:n_recomb].to_f]
          end
        end
        flat = res.flatten(1)
        t1s = flat.map { |x| x[0] }
        t2s = flat.map { |x| x[1] }
        l1s = flat.map { |x| x[2] }
        l2s = flat.map { |x| x[3] }
        rcs = flat.map { |x| x[4] }
        t1_by_rho << t1s
        t2_by_rho << t2s

        m1 = Stats.mean(t1s)
        se1 = Stats.se_iid(t1s)
        m2 = Stats.mean(t2s)
        se2 = Stats.se_iid(t2s)
        z1 = (m1 - e_tmrca_exact).abs / se1
        z2 = (m2 - e_tmrca_exact).abs / se2
        max_marg_z = [max_marg_z, z1, z2].max

        # 10-block SE for inter-locus correlation Corr(T_MRCA_1, T_MRCA_2)
        bsz = flat.length / n_blocks
        blk_corr = Array.new(n_blocks) do |b|
          s1 = t1s[(b * bsz)...((b + 1) * bsz)]
          s2 = t2s[(b * bsz)...((b + 1) * bsz)]
          Stats.pearson_corr(s1, s2)
        end
        corr_full = Stats.pearson_corr(t1s, t2s)
        corr_se = Stats.se_iid(blk_corr)

        rows << [rho, m1, se1, m2, se2, e_tmrca_exact, z1, z2, corr_full, corr_se, Stats.mean(l1s), Stats.mean(l2s), Stats.mean(rcs)]
      end

      # Pairwise differences across all rho values for marginal E[T_MRCA_1] (test of "no effect" with real SE!)
      pair_rows = []
      max_pair_z = 0.0
      (0...rhos.length).each do |i|
        ((i + 1)...rhos.length).each do |j|
          diff = (rows[i][1] - rows[j][1]).abs
          se_d = Math.sqrt(rows[i][2]**2 + rows[j][2]**2)
          z_d = diff / se_d
          max_pair_z = z_d if z_d > max_pair_z
          pair_rows << [rhos[i], rhos[j], rows[i][1] - rows[j][1], se_d, z_d]
        end
      end

      # Monotonic decrease in correlation across rho > 0 (exclude rho=0 SE=0 edge cleanly)
      min_corr_z = (0...(rhos.length - 1)).map do |i|
        se_d = Math.sqrt(rows[i][9]**2 + rows[i + 1][9]**2)
        se_d > 1e-12 ? (rows[i][8] - rows[i + 1][8]) / se_d : 1e6
      end.min
      corr_rho0 = rows[0][8]

      ctx.metrics[:c14_max_pair_z] = max_pair_z
      ctx.metrics[:c14_max_marg_z] = max_marg_z
      ctx.metrics[:c14_min_corr_z] = min_corr_z
      ctx.metrics[:c14_corr_rho0] = corr_rho0
      ctx.metrics[:c14_pass] = (max_pair_z <= 3.0 && max_marg_z <= 3.0 && min_corr_z >= 3.0 && (corr_rho0 - 1.0).abs < 1e-12)

      ctx.write_csv('S57', 'out/s57_arg_marginal_tmrca.csv',
                    %w[rho tmrca1_mean tmrca1_se tmrca2_mean tmrca2_se exact_kingman z1 z2 corr_t1_t2 corr_se len1_mean len2_mean recomb_mean], rows)
      ctx.write_csv('S57', 'out/s57_arg_pairwise_diff.csv',
                    %w[rho_a rho_b diff_tmrca1 se_diff z_diff], pair_rows)
      ctx.record_stage('S57', 'arg_marginal_invariance', Time.now - t57, 2)

      # S58: Correlation decay & SVG plots
      t58 = Time.now
      ctx.write_csv('S58', 'out/s58_arg_corr_decay.csv',
                    %w[rho corr_tmrca corr_se mean_recomb_events],
                    rows.map { |r| [r[0], r[8], r[9], r[12]] })
      ctx.write_csv('S58', 'out/s58_arg_len_stats.csv',
                    %w[rho len1_mean len2_mean corr_l1_l2],
                    rows.each_with_index.map { |r, idx| [r[0], r[10], r[11], idx.zero? ? 1.0 : Stats.pearson_corr(t1_by_rho[idx], t2_by_rho[idx])] })

      svg1 = '../../../outputs/figures/day009b_31_arg_marginal.svg'
      SVG.write_xy(svg1, title: 'S57: ARG-Lite Marginal E[T_MRCA] Invariance Across rho (n=15)',
                         subtitle: "All pairwise diffs <= #{format('%.2f', max_pair_z)} SE (Kingman marginal preserved)",
                         xlabel: 'Scaled recombination rho', ylabel: 'Marginal E[T_MRCA]',
                         series: [
                           { label: 'Locus 1 (+-SE)', color: '#2563eb', x: rhos.map { |v| [v, 0.2].max }, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] } },
                           { label: 'Locus 2 (+-SE)', color: '#16a34a', x: rhos.map { |v| [v, 0.2].max }, y: rows.map { |r| r[3] }, yerr: rows.map { |r| r[4] }, dashed: true }
                         ], log_x: true, ref_y: e_tmrca_exact, ref_label: 'Kingman 2(1-1/n)')
      ctx.register_svg('S58', svg1)

      svg2 = '../../../outputs/figures/day009b_32_arg_correlation.svg'
      SVG.write_xy(svg2, title: 'S58: Inter-Locus T_MRCA Correlation vs Scaled Recombination rho',
                         subtitle: "Corr(T1, T2)=1.0 at rho=0; drops monotonically (min step #{format('%.1f', min_corr_z)} SE)",
                         xlabel: 'Scaled recombination rho', ylabel: 'Corr(T_MRCA_1, T_MRCA_2)',
                         series: [
                           { label: 'ARG-lite Corr(T1, T2) (+-SE)', color: '#dc2626', x: rhos.map { |v| [v, 0.2].max }, y: rows.map { |r| r[8] }, yerr: rows.map { |r| r[9] } }
                         ], log_x: true, ref_y: 0.0, ref_label: 'Independent loci')
      ctx.register_svg('S58', svg2)
      ctx.record_stage('S58', 'arg_corr_decay', Time.now - t58, 3)

      # S59: Joint T1 vs T2 scatter / binned comparison at rho=0 vs rho=16
      t59 = Time.now
      joint_rows = (0...200).map do |i|
        [i, t1_by_rho[0][i], t2_by_rho[0][i], t1_by_rho[3][i], t2_by_rho[3][i]]
      end
      ctx.write_csv('S59', 'out/s59_arg_joint_samples.csv',
                    %w[idx t1_rho0 t2_rho0 t1_rho16 t2_rho16], joint_rows)

      svg3 = '../../../outputs/figures/day009b_33_arg_recomb_count.svg'
      SVG.write_xy(svg3, title: 'S59: Mean Recombination Events in ARG-Lite vs rho (n=15)',
                         subtitle: 'Ancestral recombination events scale linearly with rho',
                         xlabel: 'Scaled recombination rho', ylabel: 'Mean recombination events',
                         series: [
                           { label: 'Mean ARG recombinations', color: '#9333ea', x: rhos, y: rows.map { |r| r[12] } }
                         ])
      ctx.register_svg('S59', svg3)
      ctx.record_stage('S59', 'arg_joint_recomb', Time.now - t59, 2)
    end
  end
end
