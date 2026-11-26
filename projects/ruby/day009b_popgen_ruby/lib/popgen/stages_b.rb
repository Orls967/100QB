# frozen_string_literal: true

require_relative 'rng'
require_relative 'wf'
require_relative 'moran'
require_relative 'markov'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group B: Neutral Genetic Drift Stages S11..S20
  module StagesB
    module_function

    def run(ctx)
      s11_s12_het_decay(ctx)
      s13_s14_neutral_fixation(ctx)
      s15_s16_fixation_times(ctx)
      s17_fixation_time_dist(ctx)
      s18_s19_moran_vs_wf(ctx)
      s20_large_n_negctrl(ctx)
    end

    def run_het_ensemble(n_copies, p0, gens, sample_pts, n_rep, stage_id, workers)
      tasks = (0...n_rep).to_a
      # Parallelize across workers in batches to minimize fork overhead
      batch_sz = (n_rep.to_f / workers).ceil
      batches = tasks.each_slice(batch_sz).to_a
      res_batches = RNG.parallel_map(batches, workers: workers) do |b_ids, b_idx|
        rng = RNG.for_stage(stage_id, b_idx)
        b_ids.map do
          cnt = (n_copies * p0).round
          vals = {}
          vals[0] = 2.0 * (cnt.to_f / n_copies) * (1.0 - cnt.to_f / n_copies)
          (1..gens).each do |g|
            cnt = WF.step_biallelic(cnt, n_copies, rng)
            if sample_pts.include?(g)
              p = cnt.to_f / n_copies
              vals[g] = 2.0 * p * (1.0 - p)
            end
          end
          vals
        end
      end
      all_runs = res_batches.flatten(1)
      sample_pts.map do |t|
        xs = all_runs.map { |h| h[t] }
        [t, Stats.mean(xs), Stats.se_iid(xs)]
      end
    end

    def s11_s12_het_decay(ctx)
      t11 = Time.now
      n_rep = ctx.small? ? 600 : 2_000
      # Haploid N=50 -> H_t = H_0 (1 - 1/50)^t; Diploid N=100 (2N=200) -> H_t = H_0 (1 - 1/200)^t
      pts50 = [0, 10, 25, 50, 75, 100]
      res50 = run_het_ensemble(50, 0.5, 100, pts50, n_rep, 11, ctx.workers)
      h0 = 0.5
      rows50 = res50.map do |t, m, se|
        th = h0 * ((1.0 - 1.0 / 50.0)**t)
        z = se > 1e-12 ? (m - th).abs / se : 0.0
        [t, m, se, th, z]
      end
      ctx.write_csv('S11', 'out/s11_het_decay_n50.csv', %w[t h_mean h_se h_theory z_score], rows50)
      ctx.write_csv('S11', 'out/s11_het_ratio_n50.csv', %w[t h_ratio_theory h_ratio_sim],
                    rows50.map { |r| [r[0], (1.0 - 1.0 / 50.0)**r[0], r[1] / 0.5] })
      ctx.record_stage('S11', 'het_decay_haploid_n50', Time.now - t11, 2)

      t12 = Time.now
      pts200 = [0, 20, 50, 100, 150, 200]
      res200 = run_het_ensemble(200, 0.5, 200, pts200, n_rep, 12, ctx.workers)
      rows200 = res200.map do |t, m, se|
        th = h0 * ((1.0 - 1.0 / 200.0)**t)
        z = se > 1e-12 ? (m - th).abs / se : 0.0
        [t, m, se, th, z]
      end
      ctx.write_csv('S12', 'out/s12_het_decay_n200.csv', %w[t h_mean h_se h_theory z_score], rows200)

      # Evaluate c4 on 5 non-zero time points for both N=50 and N=200
      max_z_c4 = (rows50[1..5] + rows200[1..5]).map { |r| r[4] }.max
      ctx.metrics[:c4_max_z] = max_z_c4
      ctx.metrics[:c4_pass] = (max_z_c4 <= 3.0)

      svg1 = '../../../outputs/figures/day009b_05_het_decay.svg'
      SVG.write_xy(svg1, title: 'S11-S12: Heterozygosity Decay Under Neutral Drift',
                         subtitle: "Haploid N=50 vs Diploid 2N=200 (R=#{n_rep} reps, max z=#{format('%.2f', max_z_c4)})",
                         xlabel: 'Generation t', ylabel: 'Mean heterozygosity H_t',
                         series: [
                           { label: 'N=50 sim', color: '#dc2626', x: rows50.map { |r| r[0] }, y: rows50.map { |r| r[1] }, yerr: rows50.map { |r| r[2] }, line: false },
                           { label: 'N=50 theory', color: '#dc2626', x: rows50.map { |r| r[0] }, y: rows50.map { |r| r[3] }, points: false, dashed: true },
                           { label: '2N=200 sim', color: '#2563eb', x: rows200.map { |r| r[0] }, y: rows200.map { |r| r[1] }, yerr: rows200.map { |r| r[2] }, line: false },
                           { label: '2N=200 theory', color: '#2563eb', x: rows200.map { |r| r[0] }, y: rows200.map { |r| r[3] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S12', svg1)
      ctx.record_stage('S12', 'het_decay_diploid_n200', Time.now - t12, 2)
    end

    def s13_s14_neutral_fixation(ctx)
      # c2 requires >= 20,000 replications for p0 in {0.1, 0.5}
      t13 = Time.now
      n_copies = 20
      p_mat = Markov.wf_transition_matrix_float(n_copies)
      u_exact = Markov.fixation_probs(p_mat, rational: false)
      n_rep = ctx.small? ? 4_000 : 20_000

      rows = []
      max_z = 0.0
      [0.1, 0.3, 0.5, 0.7].each_with_index do |p0, idx|
        i0 = (n_copies * p0).round
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(13 + idx, b_i)
          fixed_cnt = 0
          t_sum = 0.0
          b_ids.each do
            out = WF.simulate_to_absorption(n_copies, i0, rng, max_gens: 10_000)
            if out[:fixed]
              fixed_cnt += 1
              t_sum += out[:gens]
            end
          end
          [fixed_cnt, t_sum]
        end
        tot_fix = res.sum { |fc, _| fc }
        p_sim = tot_fix.to_f / n_rep
        se_p = Stats.se_proportion(p_sim, n_rep)
        z = se_p > 1e-12 ? (p_sim - p0).abs / se_p : 0.0
        max_z = z if (p0 == 0.1 || p0 == 0.5) && z > max_z
        rows << [n_copies, p0, n_rep, tot_fix, p_sim, se_p, u_exact[i0], z]
      end
      ctx.metrics[:c2_max_z] = max_z
      ctx.metrics[:c2_pass] = (max_z <= 3.0)

      ctx.write_csv('S13', 'out/s13_neutral_fixation.csv',
                    %w[n_copies p0 n_rep fixed_count p_sim p_se u_markov z_score], rows)
      ctx.write_csv('S13', 'out/s13_fixation_ci.csv',
                    %w[p0 ci95_lo ci95_hi u_markov],
                    rows.map { |r| [r[1], r[4] - 1.96 * r[5], r[4] + 1.96 * r[5], r[6]] })
      ctx.record_stage('S13', 'neutral_fixation_prob', Time.now - t13, 2)

      t14 = Time.now
      svg1 = '../../../outputs/figures/day009b_06_neutral_fixation.svg'
      SVG.write_xy(svg1, title: 'S13-S14: Neutral Fixation Probability vs Initial Frequency p0',
                         subtitle: "N=#{n_copies}, R=#{n_rep} replications per p0 (P_fix = p0)",
                         xlabel: 'Initial frequency p0', ylabel: 'Fixation probability',
                         series: [
                           { label: 'WF simulation (+-SE)', color: '#2563eb', x: rows.map { |r| r[1] }, y: rows.map { |r| r[4] }, yerr: rows.map { |r| r[5] }, line: false },
                           { label: 'Exact Markov u_i = p0', color: '#dc2626', x: rows.map { |r| r[1] }, y: rows.map { |r| r[6] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S14', svg1)
      chain_rows = (1...n_copies).map { |i| [i, i.to_f / n_copies, u_exact[i], (u_exact[i] - i.to_f / n_copies).abs] }
      ctx.write_csv('S14', 'out/s14_fixation_chain_vs_p0.csv',
                    %w[k p0 u_markov abs_err_p0], chain_rows)
      ctx.record_stage('S14', 'fixation_chain_sweep', Time.now - t14, 2)
    end

    def s15_s16_fixation_times(ctx)
      t15 = Time.now
      n_copies = 20
      p_mat = Markov.wf_transition_matrix_float(n_copies)
      t_cond_exact = Markov.mean_conditional_fixation_times(p_mat, rational: false)
      t_abs_exact = Markov.mean_absorption_times(p_mat, rational: false)
      n_rep = ctx.small? ? 4_000 : 16_000

      p0_list = [0.1, 0.25, 0.5, 0.75, 0.9]
      rows = []
      max_z_c5 = 0.0
      max_rel_diff = 0.0

      p0_list.each_with_index do |p0, idx|
        i0 = (n_copies * p0).round
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(15 + idx, b_i)
          fix_t = []
          all_t = []
          b_ids.each do
            out = WF.simulate_to_absorption(n_copies, i0, rng, max_gens: 10_000)
            all_t << out[:gens].to_f
            fix_t << out[:gens].to_f if out[:fixed]
          end
          [fix_t, all_t]
        end
        fix_ts = res.flat_map { |f, _| f }
        all_ts = res.flat_map { |_, a| a }
        m_fix = Stats.mean(fix_ts)
        se_fix = Stats.se_iid(fix_ts)
        m_abs = Stats.mean(all_ts)
        se_abs = Stats.se_iid(all_ts)
        # Kimura-Ohta haploid diffusion approximation: E[T_fix | fix] = -2*N*((1-p0)/p0)*ln(1-p0)
        t_diff = -2.0 * n_copies * ((1.0 - p0) / p0) * Math.log(1.0 - p0)
        z_chain = se_fix > 1e-12 ? (m_fix - t_cond_exact[i0]).abs / se_fix : 0.0
        rel_diff = (m_fix - t_diff).abs / t_diff
        max_z_c5 = z_chain if z_chain > max_z_c5
        max_rel_diff = rel_diff if rel_diff > max_rel_diff

        rows << [p0, m_fix, se_fix, t_cond_exact[i0], t_diff, z_chain, rel_diff, m_abs, se_abs, t_abs_exact[i0]]
      end

      ctx.metrics[:c5_max_z] = max_z_c5
      ctx.metrics[:c5_max_diff_rel] = max_rel_diff
      ctx.metrics[:c5_pass] = (max_z_c5 <= 3.0 && max_rel_diff <= 0.10)

      ctx.write_csv('S15', 'out/s15_fixation_times.csv',
                    %w[p0 t_fix_sim t_fix_se t_fix_markov t_fix_diffusion z_chain rel_diff_diff t_abs_sim t_abs_se t_abs_markov], rows)
      ctx.write_csv('S15', 'out/s15_abs_vs_cond_ratio.csv',
                    %w[p0 t_abs_markov t_fix_cond_markov ratio_cond_abs],
                    rows.map { |r| [r[0], r[9], r[3], r[3] / r[9]] })
      ctx.record_stage('S15', 'fixation_time_vs_markov', Time.now - t15, 2)

      t16 = Time.now
      svg1 = '../../../outputs/figures/day009b_07_fixation_time.svg'
      SVG.write_xy(svg1, title: 'S15-S16: Conditional Fixation Time E[T_fix | fix] (N=20)',
                         subtitle: "Sim vs exact Markov chain (max z=#{format('%.2f', max_z_c5)}) & diffusion approx",
                         xlabel: 'Initial frequency p0', ylabel: 'Mean conditional fixation time (gens)',
                         series: [
                           { label: 'WF simulation (+-SE)', color: '#2563eb', x: p0_list, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] }, line: false },
                           { label: 'Exact Markov chain', color: '#16a34a', x: p0_list, y: rows.map { |r| r[3] }, points: false },
                           { label: 'Diffusion (lit memory)', color: '#dc2626', x: p0_list, y: rows.map { |r| r[4] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S16', svg1)
      diff_rows = (1...n_copies).map do |i|
        p0 = i.to_f / n_copies
        td = -2.0 * n_copies * ((1.0 - p0) / p0) * Math.log(1.0 - p0)
        [p0, t_cond_exact[i], td, (t_cond_exact[i] - td).abs / td]
      end
      ctx.write_csv('S16', 'out/s16_diffusion_vs_markov.csv',
                    %w[p0 t_fix_markov t_fix_diffusion rel_gap], diff_rows)
      ctx.record_stage('S16', 'fixation_time_diffusion', Time.now - t16, 2)
    end

    def s17_fixation_time_dist(ctx)
      t0 = Time.now
      n_copies = 20
      n_rep = ctx.small? ? 2_500 : 8_000
      rng = RNG.for_stage(17, 0)
      fix_t = []
      n_rep.times do
        out = WF.simulate_to_absorption(n_copies, 10, rng, max_gens: 10_000)
        fix_t << out[:gens].to_f if out[:fixed]
      end
      n_bins = 20
      tmax = 120.0
      dt = tmax / n_bins
      hist = Array.new(n_bins, 0.0)
      fix_t.each do |t|
        b = (t / dt).floor
        hist[b] += 1.0 / (fix_t.length * dt) if b >= 0 && b < n_bins
      end
      centers = Array.new(n_bins) { |b| (b + 0.5) * dt }
      rows = centers.each_with_index.map { |c, b| [b, c, hist[b]] }
      ctx.write_csv('S17', 'out/s17_fix_time_hist.csv', %w[bin t_center pdf_fix_time], rows)

      svg1 = '../../../outputs/figures/day009b_08_fix_time_dist.svg'
      SVG.write_hist(svg1, title: 'S17: Conditional Fixation Time Distribution (N=20, p0=0.5)',
                           subtitle: "Skewed right-tailed absorption distribution (mean=#{format('%.2f', Stats.mean(fix_t))} gens)",
                           xlabel: 'Fixation time (generations)', ylabel: 'Probability density',
                           centers: centers, counts: hist)
      ctx.register_svg('S17', svg1)
      ctx.record_stage('S17', 'fixation_time_dist', Time.now - t0, 2)
    end

    def s18_s19_moran_vs_wf(ctx)
      t18 = Time.now
      # Compare neutral Moran (N=20, time unit = N^2/2 steps? Wait: Var(dp) per step = 2p(1-p)/N^2;
      # in haploid WF, Var(dp) per gen = p(1-p)/N => 1 WF gen = N/2 Moran steps, so T_fix in Moran steps = N/2 * T_fix_WF = N^2!)
      # Let's verify both Moran mean conditional fixation time = N(N-1) when i=1? Actually exact Moran at i=N/2:
      # sum_{k=1}^{N/2} ... Let's compare Moran exact Markov chain time / (N/2) vs WF!
      n = 16
      p_moran = Markov.moran_transition_matrix(n, 1.0)
      t_moran_exact = Markov.mean_conditional_fixation_times(p_moran, rational: false)
      p_wf = Markov.wf_transition_matrix_float(n)
      t_wf_exact = Markov.mean_conditional_fixation_times(p_wf, rational: false)

      rows = (1...n).map do |i|
        # Scale Moran steps by 2/N -> generations equivalent
        t_m_scaled = t_moran_exact[i] * (2.0 / n)
        [i, i.to_f / n, t_moran_exact[i], t_m_scaled, t_wf_exact[i], (t_m_scaled - t_wf_exact[i]).abs / t_wf_exact[i]]
      end
      ctx.write_csv('S18', 'out/s18_moran_vs_wf_times.csv',
                    %w[k p0 t_moran_steps t_moran_scaled_gens t_wf_gens rel_diff], rows)
      ctx.record_stage('S18', 'moran_vs_wf_timescale', Time.now - t18, 1)

      t19 = Time.now
      # Heterozygosity decay comparison: Moran N=40 over 1600 steps vs WF N=40 over 80 gens (1 gen = 20 Moran steps = N/2)
      n_rep = ctx.small? ? 300 : 1_000
      rng = RNG.for_stage(19, 0)
      n_m = 40
      steps_m = 1600
      sample_m = 200 # every 200 Moran steps = every 10 WF gens
      m_acc = Array.new(9) { [] }
      n_rep.times do
        hs = Moran.simulate_heterozygosity(n_m, n_m / 2, steps_m, sample_m, rng)
        hs.each_with_index { |v, idx| m_acc[idx] << v }
      end
      h_rows = (0..8).map do |idx|
        gen_eq = idx * 10
        h_th = 0.5 * ((1.0 - 1.0 / n_m)**gen_eq)
        [idx * sample_m, gen_eq, Stats.mean(m_acc[idx]), Stats.se_iid(m_acc[idx]), h_th]
      end
      ctx.write_csv('S19', 'out/s19_moran_het_decay.csv',
                    %w[moran_steps wf_gen_eq h_moran_mean h_moran_se h_wf_theory], h_rows)

      svg1 = '../../../outputs/figures/day009b_09_moran_vs_wf.svg'
      SVG.write_xy(svg1, title: 'S18-S19: Moran vs Wright-Fisher Timescale Equivalence (N=40)',
                         subtitle: '1 WF generation = N/2 Moran birth-death steps (matching drift variance)',
                         xlabel: 'Equivalent WF generation t = 2 * steps / N', ylabel: 'Heterozygosity H_t',
                         series: [
                           { label: 'Moran simulation (+-SE)', color: '#2563eb', x: h_rows.map { |r| r[1] }, y: h_rows.map { |r| r[2] }, yerr: h_rows.map { |r| r[3] }, line: false },
                           { label: 'WF theory (1-1/N)^t / 2', color: '#dc2626', x: h_rows.map { |r| r[1] }, y: h_rows.map { |r| r[4] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S19', svg1)
      ctx.record_stage('S19', 'moran_het_decay', Time.now - t19, 2)
    end

    def s20_large_n_negctrl(ctx)
      t0 = Time.now
      # Negative control: very large population (N = 100,000, p0 = 0.5) never fixes within 200 gens
      n_copies = 100_000
      n_rep = ctx.small? ? 40 : 100
      rng = RNG.for_stage(20, 0)
      rows = []
      fixed_cnt = 0
      n_rep.times do |r|
        cnt = n_copies / 2
        200.times { cnt = WF.step_biallelic(cnt, n_copies, rng) }
        pf = cnt.to_f / n_copies
        fixed_cnt += 1 if cnt.zero? || cnt == n_copies
        rows << [r, n_copies, pf, (pf - 0.5).abs] if r < 30
      end
      ctx.write_csv('S20', 'out/s20_large_n_negctrl.csv', %w[rep n_copies p_after_200g abs_shift], rows)
      ctx.write_csv('S20', 'out/s20_large_n_summary.csv', %w[n_copies n_rep absorbed_count], [[n_copies, n_rep, fixed_cnt]])
      ctx.record_stage('S20', 'large_n_negctrl', Time.now - t0, 2)
    end
  end
end
