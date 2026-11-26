# frozen_string_literal: true

require_relative 'rng'
require_relative 'ld_sweep'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group F: Selective Sweep Hitchhiking & Two-Locus Linkage Disequilibrium (S51..S56)
  module StagesF
    module_function

    def run(ctx)
      s51_s53_hitchhiking(ctx)
      s54_s56_ld_decay_and_eq(ctx)
    end

    def s51_s53_hitchhiking(ctx)
      t51 = Time.now
      n_copies = 400
      s = 0.25
      n_sample = 15
      theta = 4.0
      rs = [0.0, 0.004, 0.015, 0.05, 0.25]
      n_rep = ctx.small? ? 400 : 1_200

      # Pre-generate n_rep conditioned sweep trajectories in parallel
      batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
      trajs = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
        rng = RNG.for_stage(51, b_i)
        b_ids.map { LDSweep.generate_sweep_trajectory(n_copies, s, rng) }
      end.flatten(1)

      # Neutral control (no sweep)
      rng_ctrl = RNG.for_stage(51, 999)
      ctrl_pis = Array.new(n_rep) { LDSweep.simulate_hitchhiking_diversity(n_copies, n_sample, 0.0, 0.0, theta, rng_ctrl)[:pi] }
      ctrl_m = Stats.mean(ctrl_pis)
      ctrl_se = Stats.se_iid(ctrl_pis)

      rows = []
      rs.each_with_index do |r, r_idx|
        t_batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(t_batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(52 + r_idx, b_i)
          b_ids.map { |tid| LDSweep.simulate_hitchhiking_diversity(n_copies, n_sample, s, r, theta, rng, sweep_traj: trajs[tid]) }
        end.flatten(1)
        pis = res.map { |x| x[:pi] }
        ss = res.map { |x| x[:s].to_f }
        m = Stats.mean(pis)
        se = Stats.se_iid(pis)
        z_vs_ctrl = (ctrl_m - m) / Math.sqrt(ctrl_se**2 + se**2)
        rows << [r, m, se, m / ctrl_m, Stats.mean(ss), z_vs_ctrl]
      end

      # Monotonic recovery across r: r[0] < r[1] < r[2] < r[3] < r[4]
      min_rec_z = (0...(rs.length - 1)).map do |i|
        (rows[i + 1][1] - rows[i][1]) / Math.sqrt(rows[i + 1][2]**2 + rows[i][2]**2)
      end.min
      z_drop_r0 = rows[0][5]
      z_negctrl_large_r = rows.last[5].abs

      ctx.metrics[:c12_drop_z] = z_drop_r0
      ctx.metrics[:c12_min_rec_z] = min_rec_z
      ctx.metrics[:c12_large_r_z] = z_negctrl_large_r
      ctx.metrics[:c12_pass] = (z_drop_r0 >= 3.0 && min_rec_z >= 3.0 && z_negctrl_large_r < 3.0)

      ctx.write_csv('S51', 'out/s51_sweep_diversity.csv',
                    %w[r pi_mean pi_se pi_rel_ctrl s_mean z_below_ctrl], rows)
      ctx.write_csv('S51', 'out/s51_sweep_control.csv',
                    %w[n_copies s n_sample theta pi_ctrl_mean pi_ctrl_se],
                    [[n_copies, s, n_sample, theta, ctrl_m, ctrl_se]])
      ctx.write_csv('S51', 'out/s51_sweep_duration.csv',
                    %w[n_copies s mean_fix_gens median_fix_gens],
                    [[n_copies, s, Stats.mean(trajs.map { |tr| tr.length.to_f }), Stats.median_mad(trajs.map { |tr| tr.length.to_f })[0]]])
      ctx.record_stage('S51', 'sweep_hitchhiking', Time.now - t51, 3)

      t52 = Time.now
      svg1 = '../../../outputs/figures/day009b_27_hitchhiking.svg'
      SVG.write_xy(svg1, title: 'S51-S53: Selective Sweep Diversity Trough vs Recombination Distance r',
                         subtitle: "N=#{n_copies}, s=#{s}: r=0 drop #{format('%.1f', z_drop_r0)} SE; recovery min #{format('%.1f', min_rec_z)} SE; r=0.25 matches neutral",
                         xlabel: 'Recombination probability r per gen', ylabel: 'Neutral diversity pi',
                         series: [
                           { label: 'Post-sweep pi (+-SE)', color: '#2563eb', x: rs.map { |v| [v, 5e-4].max }, y: rows.map { |x| x[1] }, yerr: rows.map { |x| x[2] } }
                         ], log_x: true, ref_y: ctrl_m, ref_label: 'Neutral control pi')
      ctx.register_svg('S52', svg1)
      ctx.write_csv('S52', 'out/s52_sweep_traj_sample.csv',
                    %w[gen x_beneficial], trajs[0].each_with_index.map { |x, t| [t, x] })
      ctx.record_stage('S52', 'sweep_recovery_plot', Time.now - t52, 2)

      t53 = Time.now
      svg2 = '../../../outputs/figures/day009b_28_sweep_traj.svg'
      SVG.write_xy(svg2, title: 'S53: Sample Conditioned Beneficial Sweep Trajectory x(t)',
                         subtitle: "N=#{n_copies}, s=#{s} (conditioned on fixation from 1/N to 1)",
                         xlabel: 'Generation t', ylabel: 'Beneficial allele frequency x(t)',
                         series: [
                           { label: 'Sweep #0', color: '#dc2626', x: (0...trajs[0].length).to_a, y: trajs[0], points: false },
                           { label: 'Sweep #1', color: '#2563eb', x: (0...trajs[1].length).to_a, y: trajs[1], points: false }
                         ])
      ctx.register_svg('S53', svg2)
      ctx.write_csv('S53', 'out/s53_sweep_negctrl_large_r.csv',
                    %w[r pi_mean pi_se pi_rel_ctrl s_mean z_below_ctrl], [rows.last])
      ctx.record_stage('S53', 'sweep_large_r_ctrl', Time.now - t53, 2)
    end

    def s54_s56_ld_decay_and_eq(ctx)
      # S54: Deterministic LD decay D_t = D_0 (1 - r)^t vs large finite-N WF (N_copies = 4000, r = 0.05)
      t54 = Time.now
      n_copies = 4_000
      # Start with complete coupling D0 = 0.25: [AB=0.5, Ab=0, aB=0, ab=0.5]
      init_c = [n_copies / 2, 0, 0, n_copies / 2]
      r_ld = 0.05
      gens = 40
      sample_t = [0, 5, 10, 15, 20, 30, 40]
      n_rep = ctx.small? ? 150 : 500

      batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
      res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
        rng = RNG.for_stage(54, b_i)
        b_ids.map { LDSweep.simulate_ld_wf(n_copies, init_c, r_ld, gens, rng, sample_times: sample_t) }
      end.flatten(1)

      d_det = LDSweep.deterministic_decay(0.25, r_ld, gens)
      rows = []
      max_z = 0.0
      sample_t.each do |t|
        ds = res.map { |h| h[t][:d] }
        m = Stats.mean(ds)
        se = Stats.se_iid(ds)
        z = (t.zero? || se < 1e-12) ? 0.0 : (m - d_det[t]).abs / se
        max_z = z if z > max_z
        rows << [t, m, se, d_det[t], z]
      end

      ctx.write_csv('S54', 'out/s54_ld_decay.csv', %w[t d_mean d_se d_det z_score], rows)
      ctx.write_csv('S54', 'out/s54_ld_half_life.csv', %w[r t_half_det d_at_t15],
                    [0.01, 0.02, 0.05, 0.10, 0.20].map { |rv| [rv, Math.log(0.5) / Math.log(1.0 - rv), 0.25 * ((1.0 - rv)**15)] })
      ctx.record_stage('S54', 'ld_decay_det_vs_wf', Time.now - t54, 2)

      # S55: Stationary sigma_d^2 = E[D^2] / E[p_A(1-p_A)p_B(1-p_B)] vs Ohta-Kimura (10+rho)/(22+13*rho+rho^2) & Sved 1/(1+rho)
      t55 = Time.now
      n_eq = 100
      u_mut = 0.001
      rhos = [0.5, 2.0, 6.0, 16.0]
      eq_gens = ctx.small? ? 6_000 : 36_000
      burn = 500
      n_blocks = 10

      eq_rows = RNG.parallel_map((0...rhos.length).to_a, workers: ctx.workers) do |r_idx, _|
        rho = rhos[r_idx]
        rng_eq = RNG.for_stage(55, r_idx)
        r = rho / (2.0 * n_eq)
        counts = [n_eq / 2, 0, 0, n_eq / 2]
        burn.times do
          freqs = counts.map { |c| c.to_f / n_eq }
          fab, fa_b, f_ab, f_a_b = freqs
          pa = fab + fa_b
          pb = fab + f_ab
          d = fab - pa * pb
          nf = [fab - r * d, fa_b + r * d, f_ab + r * d, f_a_b - r * d].map { |v| [v, 0.0].max }
          m0 = 1.0 - u_mut
          m1 = u_mut
          mf = Array.new(4, 0.0)
          4.times do |src|
            sa = (src >> 1) & 1
            sb = src & 1
            4.times do |dst|
              da = (dst >> 1) & 1
              db = dst & 1
              mf[dst] += nf[src] * (sa == da ? m0 : m1) * (sb == db ? m0 : m1)
            end
          end
          s_mf = mf.sum(0.0)
          counts = Samplers.multinomial(n_eq, mf.map { |v| v / s_mf }, rng_eq)
        end

        bsz = eq_gens / n_blocks
        blk_sig = Array.new(n_blocks) do
          d2_sum = 0.0
          den_sum = 0.0
          bsz.times do
            freqs = counts.map { |c| c.to_f / n_eq }
            fab, fa_b, f_ab, f_a_b = freqs
            pa = fab + fa_b
            pb = fab + f_ab
            d = fab - pa * pb
            nf = [fab - r * d, fa_b + r * d, f_ab + r * d, f_a_b - r * d].map { |v| [v, 0.0].max }
            m0 = 1.0 - u_mut
            m1 = u_mut
            mf = Array.new(4, 0.0)
            4.times do |src|
              sa = (src >> 1) & 1
              sb = src & 1
              4.times do |dst|
                da = (dst >> 1) & 1
                db = dst & 1
                mf[dst] += nf[src] * (sa == da ? m0 : m1) * (sb == db ? m0 : m1)
              end
            end
            s_mf = mf.sum(0.0)
            counts = Samplers.multinomial(n_eq, mf.map { |v| v / s_mf }, rng_eq)
            st = LDSweep.ld_stats(counts.map { |c| c.to_f / n_eq })
            d2_sum += st[:d] * st[:d]
            den_sum += st[:pa] * (1.0 - st[:pa]) * st[:pb] * (1.0 - st[:pb])
          end
          den_sum > 1e-12 ? (d2_sum / den_sum) : 0.0
        end

        sig_m = Stats.mean(blk_sig)
        sig_se = Stats.se_iid(blk_sig)
        sved = 1.0 / (1.0 + rho)
        ohta_kimura = (10.0 + rho) / (22.0 + 13.0 * rho + rho * rho)
        rel_ok = (sig_m - ohta_kimura).abs / ohta_kimura
        rel_sv = (sig_m - sved).abs / sved
        [rho, r, sig_m, sig_se, sved, ohta_kimura, rel_sv, rel_ok]
      end

      max_rel_lit = eq_rows.map { |row| row[7] }.max
      ctx.metrics[:c13_decay_max_z] = max_z
      ctx.metrics[:c13_r2_max_rel] = max_rel_lit
      ctx.metrics[:c13_pass] = (max_z <= 3.0 && max_rel_lit <= 0.25)

      ctx.write_csv('S55', 'out/s55_ld_eq_r2.csv',
                    %w[rho r r2_mean r2_se sved_lit ohta_kimura_lit rel_sved rel_ok], eq_rows)
      ctx.write_csv('S55', 'out/s55_ld_ohta_kimura_table.csv',
                    %w[rho sved_1_over_1prho ohta_kimura_sigmad2],
                    [0.2, 0.5, 1.0, 2.0, 4.0, 6.0, 10.0, 16.0, 32.0].map { |rh| [rh, 1.0 / (1.0 + rh), (10.0 + rh) / (22.0 + 13.0 * rh + rh * rh)] })
      ctx.record_stage('S55', 'ld_equilibrium_r2', Time.now - t55, 2)

      t56 = Time.now
      svg1 = '../../../outputs/figures/day009b_29_ld_decay.svg'
      SVG.write_xy(svg1, title: 'S54: Linkage Disequilibrium Decay D_t = D_0 (1 - r)^t (N=4000, r=0.05)',
                         subtitle: "Mean of R=#{n_rep} finite-N WF trajectories vs deterministic curve (max z=#{format('%.2f', max_z)})",
                         xlabel: 'Generation t', ylabel: 'Linkage disequilibrium D_t',
                         series: [
                           { label: 'WF mean (+-SE)', color: '#2563eb', x: sample_t, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] }, line: false },
                           { label: 'Deterministic D0(1-r)^t', color: '#dc2626', x: sample_t, y: rows.map { |r| r[3] } }
                         ])
      ctx.register_svg('S56', svg1)

      svg2 = '../../../outputs/figures/day009b_30_ld_eq_r2.svg'
      SVG.write_xy(svg2, title: 'S55-S56: Equilibrium E[r^2] vs Scaled Recombination rho',
                         subtitle: "Two-locus WF (N=100) vs Sved 1/(1+rho) and Ohta-Kimura (lit memory)",
                         xlabel: 'Scaled recombination rho = 2 N r', ylabel: 'Equilibrium E[r^2]',
                         series: [
                           { label: 'WF E[r^2] (+-SE)', color: '#2563eb', x: rhos, y: eq_rows.map { |r| r[2] }, yerr: eq_rows.map { |r| r[3] } },
                           { label: 'Sved 1/(1+rho)', color: '#dc2626', x: rhos, y: eq_rows.map { |r| r[4] }, dashed: true },
                           { label: 'Ohta-Kimura', color: '#16a34a', x: rhos, y: eq_rows.map { |r| r[5] }, dashed: true }
                         ], log_x: true)
      ctx.register_svg('S56', svg2)
      ctx.write_csv('S56', 'out/s56_ld_summary.csv', %w[decay_max_z eq_r2_max_rel], [[max_z, max_rel_lit]])
      ctx.record_stage('S56', 'ld_plots', Time.now - t56, 3)
    end
  end
end
