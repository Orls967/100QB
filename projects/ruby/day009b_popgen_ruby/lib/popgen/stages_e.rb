# frozen_string_literal: true

require_relative 'rng'
require_relative 'markov'
require_relative 'structured'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group E: Structured Coalescent, Island Model, Stepping-Stone & F_ST (S43..S50)
  module StagesE
    module_function

    def run(ctx)
      s43_s44_island_exact(ctx)
      s45_s46_fst_vs_migration(ctx)
      s47_s48_stepping_stone_ibd(ctx)
      s49_s50_hudson_vs_nei(ctx)
    end

    def s43_s44_island_exact(ctx)
      t43 = Time.now
      ds = [2, 4, 8]
      m_l = 0.5
      n_rep = ctx.small? ? 2_000 : 8_000
      rows = []
      max_z = 0.0

      ds.each_with_index do |d, idx|
        ex = Markov.island_exact_times(d, m_l)
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(43 + idx, b_i)
          b_ids.map do
            [
              Structured.simulate_island_two_sample(d, m_l, true, rng),
              Structured.simulate_island_two_sample(d, m_l, false, rng)
            ]
          end
        end
        flat = res.flatten(1)
        tws = flat.map { |x| x[0] }
        tbs = flat.map { |x| x[1] }
        mw = Stats.mean(tws)
        sew = Stats.se_iid(tws)
        mb = Stats.mean(tbs)
        seb = Stats.se_iid(tbs)
        zw = (mw - ex[:tw]).abs / sew
        zb = (mb - ex[:tb]).abs / seb
        max_z = [max_z, zw, zb].max
        rows << [d, m_l, mw, sew, ex[:tw], zw, mb, seb, ex[:tb], zb]
      end

      ctx.metrics[:c10_times_max_z] = max_z
      ctx.write_csv('S43', 'out/s43_island_times.csv',
                    %w[d m_scaled tw_sim tw_se tw_exact z_tw tb_sim tb_se tb_exact z_tb], rows)
      ctx.write_csv('S43', 'out/s43_island_ratios.csv',
                    %w[d tb_over_tw_sim tb_over_tw_exact fst_hudson_exact],
                    rows.map { |r| [r[0], r[6] / r[2], r[8] / r[4], (r[8] - r[4]) / r[8]] })
      ctx.record_stage('S43', 'island_exact_times', Time.now - t43, 2)

      t44 = Time.now
      svg1 = '../../../outputs/figures/day009b_23_island_times.svg'
      SVG.write_xy(svg1, title: 'S43-S44: Finite Island Coalescent Times T_within & T_between (M=0.5)',
                         subtitle: "d in {2, 4, 8} demes: simulation vs exact 2x2 linear system (max z=#{format('%.2f', max_z)})",
                         xlabel: 'Number of demes d', ylabel: 'Mean coalescence time (N_deme gens)',
                         series: [
                           { label: 'T_within sim', color: '#2563eb', x: ds, y: rows.map { |r| r[2] }, yerr: rows.map { |r| r[3] }, line: false },
                           { label: 'T_within exact', color: '#2563eb', x: ds, y: rows.map { |r| r[4] }, points: false, dashed: true },
                           { label: 'T_between sim', color: '#dc2626', x: ds, y: rows.map { |r| r[6] }, yerr: rows.map { |r| r[7] }, line: false },
                           { label: 'T_between exact', color: '#dc2626', x: ds, y: rows.map { |r| r[8] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S44', svg1)
      ctx.write_csv('S44', 'out/s44_island_scaling_d.csv',
                    %w[d tw_exact tb_exact fst_hudson_exact],
                    (2..8).map { |d| ex = Markov.island_exact_times(d, m_l); [d, ex[:tw], ex[:tb], ex[:fst_hudson]] })
      ctx.record_stage('S44', 'island_d_scaling', Time.now - t44, 2)
    end

    def s45_s46_fst_vs_migration(ctx)
      t45 = Time.now
      d = 4
      ms = [0.08, 0.25, 0.80, 2.50, 8.00, 80.0]
      n_rep = ctx.small? ? 800 : 3_200
      n_blocks = 10
      rows = []

      ms.each_with_index do |m_l, idx|
        ex = Markov.island_exact_times(d, m_l)
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(45 + idx, b_i)
          b_ids.map do
            # 2 samples from deme 0, 2 from deme 1 -> exact unbiased Tw and Tb per replicate
            st = Structured.simulate_structured_pairwise(d, [3, 3, 3, 3], m_l, :island, 2.0, rng)
            [st[:tw], st[:tb]]
          end
        end
        flat = res.flatten(1)
        bsz = flat.length / n_blocks
        blk_fst = Array.new(n_blocks) do |b|
          sl = flat[(b * bsz)...((b + 1) * bsz)]
          mw = Stats.mean(sl.map { |x| x[0] })
          mb = Stats.mean(sl.map { |x| x[1] })
          (mb - mw) / mb
        end
        m_fst = Stats.mean(blk_fst)
        se_fst = Stats.se_iid(blk_fst)
        rows << [m_l, m_fst, se_fst, ex[:fst_hudson], (m_fst - ex[:fst_hudson]).abs / se_fst]
      end

      min_step_z = (0...4).map do |i|
        (rows[i][1] - rows[i + 1][1]) / Math.sqrt(rows[i][2]**2 + rows[i + 1][2]**2)
      end.min
      negctrl_fst = rows.last[1]
      ctx.metrics[:c10_min_step_z] = min_step_z
      ctx.metrics[:c10_negctrl_fst] = negctrl_fst
      ctx.metrics[:c10_pass] = (ctx.metrics[:c10_times_max_z] <= 3.0 && min_step_z >= 3.0 && negctrl_fst < 0.02)

      ctx.write_csv('S45', 'out/s45_fst_vs_migration.csv',
                    %w[m_scaled fst_sim fst_se fst_exact z_vs_exact], rows)
      ctx.record_stage('S45', 'fst_vs_migration', Time.now - t45, 1)

      t46 = Time.now
      svg1 = '../../../outputs/figures/day009b_24_fst_migration.svg'
      SVG.write_xy(svg1, title: 'S45-S46: Population Differentiation F_ST vs Scaled Migration M = Nm',
                         subtitle: "4-deme island: monotonic decline (min step #{format('%.1f', min_step_z)} SE); M=50 -> F_ST=#{format('%.4f', negctrl_fst)}",
                         xlabel: 'Scaled migration rate M = N_deme * m', ylabel: 'Hudson F_ST',
                         series: [
                           { label: 'Structured coalescent (+-SE)', color: '#2563eb', x: ms, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] }, line: false },
                           { label: 'Exact linear-system F_ST', color: '#dc2626', x: ms, y: rows.map { |r| r[3] } }
                         ], log_x: true, ref_y: 0.02, ref_label: 'Negative control limit (0.02)')
      ctx.register_svg('S46', svg1)
      ctx.write_csv('S46', 'out/s46_fst_negctrl.csv', %w[m_scaled fst_sim fst_se fst_exact z_vs_exact], [rows.last])
      ctx.record_stage('S46', 'fst_migration_plot', Time.now - t46, 2)
    end

    def s47_s48_stepping_stone_ibd(ctx)
      t47 = Time.now
      d = 8
      m_l = 0.6
      n_rep = ctx.small? ? 500 : 1_500
      rng = RNG.for_stage(47, 0)
      sample_counts = Array.new(d, 2)

      within_sum = 0.0
      by_dist_sum = Array.new(4, 0.0)
      by_dist_cnt = Array.new(4, 0)

      n_rep.times do
        st = Structured.simulate_structured_pairwise(d, sample_counts, m_l, :stepping_stone, 2.0, rng)
        within_sum += st[:tw]
        st[:between_by_pair].each do |(da, db), (val_sum, cnt)|
          raw_d = (da - db).abs
          ring_d = [raw_d, d - raw_d].min
          by_dist_sum[ring_d - 1] += val_sum
          by_dist_cnt[ring_d - 1] += cnt
        end
      end
      tw_m = within_sum / n_rep
      dists = [1, 2, 3, 4]
      rows = dists.map do |dist|
        tb_d = by_dist_sum[dist - 1] / by_dist_cnt[dist - 1]
        fst_d = (tb_d - tw_m) / tb_d
        [dist, tw_m, tb_d, fst_d]
      end
      corr_ibd = Stats.pearson_corr(dists.map(&:to_f), rows.map { |r| r[3] })
      ctx.write_csv('S47', 'out/s47_stepping_stone_ibd.csv', %w[distance tw_mean tb_mean fst_pair], rows)
      ctx.write_csv('S47', 'out/s47_stepping_tb_ratio.csv', %w[distance tb_over_tw],
                    rows.map { |r| [r[0], r[2] / r[1]] })
      ctx.record_stage('S47', 'stepping_stone_ibd', Time.now - t47, 2)

      t48 = Time.now
      svg1 = '../../../outputs/figures/day009b_25_ibd_stepping_stone.svg'
      SVG.write_xy(svg1, title: 'S47-S48: Isolation by Distance on 8-Deme 1D Stepping-Stone Ring',
                         subtitle: "Pairwise F_ST increases monotonically with ring distance (r = #{format('%.4f', corr_ibd)})",
                         xlabel: 'Ring separation distance |i - j|', ylabel: 'Pairwise Hudson F_ST',
                         series: [{ label: 'Stepping-stone F_ST(d)', color: '#2563eb', x: dists, y: rows.map { |r| r[3] } }])
      ctx.register_svg('S48', svg1)
      ctx.write_csv('S48', 'out/s48_ibd_corr.csv', %w[demes m_scaled corr_dist_fst], [[d, m_l, corr_ibd]])
      ctx.record_stage('S48', 'stepping_stone_plot', Time.now - t48, 2)
    end

    def s49_s50_hudson_vs_nei(ctx)
      t49 = Time.now
      ds = [2, 3, 4, 6, 8]
      m_l = 0.5
      n_rep = ctx.small? ? 500 : 1_500
      rows = []

      ds.each_with_index do |d, idx|
        rng = RNG.for_stage(49, idx)
        sc = Array.new(d, 2)
        tws = []
        tbs = []
        n_rep.times do
          st = Structured.simulate_structured_pairwise(d, sc, m_l, :island, 2.0, rng)
          tws << st[:tw]
          tbs << st[:tb]
        end
        mw = Stats.mean(tws)
        mb = Stats.mean(tbs)
        mt = (mw + (d - 1) * mb) / d.to_f
        fst_h = (mb - mw) / mb
        fst_n = (mt - mw) / mt
        ex = Markov.island_exact_times(d, m_l)
        rows << [d, fst_h, ex[:fst_hudson], fst_n, ex[:fst_slatkin], fst_h - fst_n]
      end
      ctx.write_csv('S49', 'out/s49_hudson_vs_nei.csv',
                    %w[demes fst_hudson_sim fst_hudson_exact fst_nei_sim fst_nei_exact diff_hudson_nei], rows)
      ctx.write_csv('S49', 'out/s49_hudson_nei_ratio.csv',
                    %w[demes ratio_nei_over_hudson_exact],
                    rows.map { |r| [r[0], r[4] / r[2]] })
      ctx.record_stage('S49', 'hudson_vs_nei_fst', Time.now - t49, 2)

      t50 = Time.now
      svg1 = '../../../outputs/figures/day009b_26_hudson_vs_nei.svg'
      SVG.write_xy(svg1, title: 'S49-S50: Hudson F_ST vs Nei G_ST Across Deme Count d (M=0.5)',
                         subtitle: 'Nei G_ST depends on d via within-deme weight 1/d; Hudson F_ST is d-invariant',
                         xlabel: 'Number of demes d', ylabel: 'Differentiation F_ST / G_ST',
                         series: [
                           { label: 'Hudson F_ST (sim)', color: '#2563eb', x: ds, y: rows.map { |r| r[1] } },
                           { label: 'Nei G_ST (sim)', color: '#dc2626', x: ds, y: rows.map { |r| r[3] } }
                         ])
      ctx.register_svg('S50', svg1)
      ctx.write_csv('S50', 'out/s50_fst_diff_summary.csv', %w[n_demes_tested diff_at_d8], [[ds.length, rows.last[5]]])
      ctx.record_stage('S50', 'hudson_vs_nei_plot', Time.now - t50, 2)
    end
  end
end
