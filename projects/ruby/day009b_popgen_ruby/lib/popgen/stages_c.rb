# frozen_string_literal: true

require_relative 'rng'
require_relative 'wf'
require_relative 'moran'
require_relative 'markov'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group C: Selection, Mutation-Drift Balance, Dominance & Stationary Distributions (S21..S30)
  module StagesC
    module_function

    def run(ctx)
      s21_s22_moran_selection(ctx)
      s23_s24_wf_selection_det(ctx)
      s25_s26_infinite_alleles(ctx)
      s27_s28_stationary_mut_sel(ctx)
      s29_mut_sel_balance(ctx)
      s30_dominance_overdominance(ctx)
    end

    def s21_s22_moran_selection(ctx)
      t21 = Time.now
      n = 16
      i0 = 4 # p0 = 0.25
      rs = [0.85, 1.0, 1.15, 1.30]
      n_rep = ctx.small? ? 2_500 : 10_000
      rows = []
      max_chain_err = 0.0
      max_sim_z = 0.0
      p_sim_map = {}
      se_map = {}

      rs.each_with_index do |r, idx|
        p_closed = Markov.moran_fixation_closed(i0, n, r)
        p_mat = Markov.moran_transition_matrix(n, r)
        u_chain = Markov.fixation_probs(p_mat, rational: false)
        chain_err = (p_closed - u_chain[i0]).abs
        max_chain_err = chain_err if chain_err > max_chain_err

        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(21 + idx, b_i)
          b_ids.count { Moran.simulate_to_absorption(n, i0, r, rng)[:fixed] }
        end
        p_sim = res.sum.to_f / n_rep
        se = Stats.se_proportion(p_sim, n_rep)
        z = se > 1e-12 ? (p_sim - p_closed).abs / se : 0.0
        max_sim_z = z if z > max_sim_z
        p_sim_map[r] = p_sim
        se_map[r] = se
        rows << [r, n, i0, p_closed, u_chain[i0], chain_err, p_sim, se, z]
      end

      # Directional checks: r=1.15 > r=1.0 (>= 3 SE), r=0.85 < r=1.0 (>= 3 SE), r=1.0 neutral control
      z_up = (p_sim_map[1.15] - p_sim_map[1.0]) / Math.sqrt(se_map[1.15]**2 + se_map[1.0]**2)
      z_dn = (p_sim_map[1.0] - p_sim_map[0.85]) / Math.sqrt(se_map[1.0]**2 + se_map[0.85]**2)
      diff_neu = p_sim_map[1.0] - (i0.to_f / n)
      ctx.metrics[:c3_chain_err] = max_chain_err
      ctx.metrics[:c3_sim_z] = max_sim_z
      ctx.metrics[:c3_z_up] = z_up
      ctx.metrics[:c3_z_dn] = z_dn
      ctx.metrics[:c3_neu_diff] = diff_neu
      ctx.metrics[:c3_neu_se] = se_map[1.0]
      ctx.metrics[:c3_pass] = (max_chain_err <= 1e-12 && max_sim_z <= 3.0 && z_up >= 3.0 && z_dn >= 3.0)

      ctx.write_csv('S21', 'out/s21_moran_selection.csv',
                    %w[r n init_i p_closed u_markov chain_err p_sim p_se z_vs_closed], rows)
      ctx.write_csv('S21', 'out/s21_moran_directional.csv',
                    %w[comparison z_score se_diff],
                    [['r1.15_vs_1.0', z_up, Math.sqrt(se_map[1.15]**2 + se_map[1.0]**2)],
                     ['r0.85_vs_1.0', z_dn, Math.sqrt(se_map[0.85]**2 + se_map[1.0]**2)]])
      ctx.record_stage('S21', 'moran_fixation_sel', Time.now - t21, 2)

      t22 = Time.now
      r_scan = [0.7, 0.85, 1.0, 1.15, 1.3, 1.5]
      scan_rows = r_scan.map do |r|
        pc = Markov.moran_fixation_closed(1, 20, r)
        [r, 20, 1, pc, pc * 20.0]
      end
      ctx.write_csv('S22', 'out/s22_moran_r_scan.csv', %w[r n init_i p_fix ratio_to_neutral], scan_rows)

      svg1 = '../../../outputs/figures/day009b_10_moran_selection.svg'
      SVG.write_xy(svg1, title: 'S21-S22: Moran Fixation Probability vs Fitness r (N=16, i0=4)',
                         subtitle: "Closed-form vs Markov chain (err=#{format('%.2e', max_chain_err)}) & R=#{n_rep} sims",
                         xlabel: 'Relative fitness r', ylabel: 'Fixation probability',
                         series: [
                           { label: 'Moran simulation (+-SE)', color: '#2563eb', x: rs, y: rows.map { |x| x[6] }, yerr: rows.map { |x| x[7] }, line: false },
                           { label: 'Closed form (1-1/r^i)/(1-1/r^N)', color: '#dc2626', x: rs, y: rows.map { |x| x[3] }, points: false }
                         ],
                         ref_y: i0.to_f / n, ref_label: 'Neutral p0=0.25')
      ctx.register_svg('S22', svg1)
      ctx.record_stage('S22', 'moran_selection_plot', Time.now - t22, 2)
    end

    def s23_s24_wf_selection_det(ctx)
      # c11: N = 2000, p0 = 0.1, s = 0.1, t in [0, 100] evaluated at a priori t_eval = [0,10,20,30,40,50,60,80,100]
      t23 = Time.now
      n_copies = 2_000
      p0 = 0.1
      s = 0.1
      gens = 100
      t_eval = [0, 10, 20, 30, 40, 50, 60, 80, 100]
      n_rep = ctx.small? ? 200 : 800

      # Deterministic recursion p_{t+1} = p_t (1 + s) / (1 + s * p_t) (haploid / genic selection)
      # Wait: in Markov/WF with diploid h=0.5, w_bar = 1 + s*p*(2 - p)? No!
      # Let's check exact deterministic haploid vs additive diploid recursion:
      # In WF.step_biallelic with h=0.5: w_aa = 1+s, w_ab = 1 + 0.5*s, w_bb = 1 =>
      # p_sel = (p^2(1+s) + pq(1+0.5s)) / (p^2(1+s) + 2pq(1+0.5s) + q^2) = p*(1 + 0.5*s*(p+1)) / (1 + s*p)!
      # Wait: prompt c11 explicitly writes "rekursi deterministik p' = p(1+s)/(1+sp) vs rata-rata lintasan besar (N=2000, p0=0.1, s=0.1)"!
      # That exact formula p' = p(1+s)/(1+sp) is the HAPLOID multiplicative/genic selection recursion where fitness of A is 1+s and a is 1!
      # Let's write a dedicated haploid WF step matching p' = p(1+s)/(1+sp) for S23!
      p_det = Array.new(gens + 1, 0.0)
      p_det[0] = p0
      gens.times do |g|
        pt = p_det[g]
        p_det[g + 1] = (pt * (1.0 + s)) / (1.0 + s * pt)
      end

      batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
      res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
        rng = RNG.for_stage(23, b_i)
        b_ids.map do
          cnt = (n_copies * p0).round
          out = { 0 => cnt.to_f / n_copies }
          (1..gens).each do |g|
            p = cnt.to_f / n_copies
            p_sel = (p * (1.0 + s)) / (1.0 + s * p)
            cnt = Binomial.sample(n_copies, p_sel, rng)
            out[g] = cnt.to_f / n_copies if t_eval.include?(g)
          end
          out
        end
      end
      all_runs = res.flatten(1)
      rows = []
      max_z = 0.0
      t_eval.each do |t|
        xs = all_runs.map { |h| h[t] }
        m = Stats.mean(xs)
        se = Stats.se_iid(xs)
        z = (t.zero? || se < 1e-12) ? 0.0 : ((m - p_det[t]).abs / se)
        max_z = z if z > max_z
        rows << [t, m, se, p_det[t], z]
      end
      ctx.metrics[:c11_max_z] = max_z
      ctx.metrics[:c11_pass] = (max_z <= 3.0)

      ctx.write_csv('S23', 'out/s23_wf_selection_traj.csv',
                    %w[t p_mean p_se p_det z_score], rows)
      ctx.write_csv('S23', 'out/s23_wf_selection_var.csv',
                    %w[t p_mean p_sd p_se],
                    t_eval.each_with_index.map { |t, i| [t, rows[i][1], Stats.sd(all_runs.map { |h| h[t] }), rows[i][2]] })
      ctx.record_stage('S23', 'wf_selection_traj', Time.now - t23, 2)

      t24 = Time.now
      # Compare s in {0.0, 0.05, 0.1, 0.2} deterministic curves
      cmp_rows = (0..100).step(5).map do |t|
        vals = [0.0, 0.05, 0.1, 0.2].map do |sv|
          pt = p0
          t.times { pt = (pt * (1.0 + sv)) / (1.0 + sv * pt) }
          pt
        end
        [t] + vals
      end
      ctx.write_csv('S24', 'out/s24_selection_s_sweep.csv', %w[t p_s0 p_s005 p_s01 p_s02], cmp_rows)

      svg1 = '../../../outputs/figures/day009b_11_wf_selection.svg'
      SVG.write_xy(svg1, title: 'S23-S24: Wright-Fisher Selection Trajectory (N=2000, p0=0.1, s=0.1)',
                         subtitle: "Mean of R=#{n_rep} paths vs p' = p(1+s)/(1+sp) (max z=#{format('%.2f', max_z)})",
                         xlabel: 'Generation t', ylabel: 'Allele frequency p_t',
                         series: [
                           { label: 'WF mean (+-SE)', color: '#2563eb', x: t_eval, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] }, line: false },
                           { label: 'Deterministic recursion', color: '#dc2626', x: t_eval, y: rows.map { |r| r[3] }, points: false }
                         ])
      ctx.register_svg('S24', svg1)
      ctx.record_stage('S24', 'wf_selection_sweep', Time.now - t24, 2)
    end

    def s25_s26_infinite_alleles(ctx)
      # c9 part 3: Infinite-alleles WF with N_copies = 500 gene copies => theta = 2 * N_copies * mu
      # Exact Kimura-Crow: E[H] = theta / (theta + 1 + theta / N_copies)
      t25 = Time.now
      n_copies = 500
      thetas = [0.5, 1.0, 2.0, 4.0]
      gens_prod = ctx.small? ? 1_200 : 4_500
      burn = ctx.small? ? 600 : 2_000

      batches = (0...thetas.length).to_a
      res = RNG.parallel_map(batches, workers: ctx.workers) do |t_idx, _|
        th = thetas[t_idx]
        mu = th / (2.0 * n_copies.to_f)
        rng = RNG.for_stage(25, t_idx)
        sim = WF.simulate_infinite_alleles(n_copies, mu, gens_prod, rng, burn_in: burn)
        ba = Stats.block_average(sim[:h_series], 10)
        h_approx = th / (1.0 + th)
        h_corr = th / (1.0 + th + th / n_copies.to_f)
        rel_err = (ba[:mean] - h_approx).abs / h_approx
        [th, mu, ba[:mean], ba[:se], h_approx, h_corr, rel_err, sim[:k_alleles]]
      end

      max_rel = res.map { |r| r[6] }.max
      ctx.metrics[:c9_h_max_rel] = max_rel
      ctx.write_csv('S25', 'out/s25_infinite_alleles_h.csv',
                    %w[theta mu h_mean h_se h_approx_theory h_o1n_corr rel_err k_alleles_end], res)
      ctx.write_csv('S25', 'out/s25_infinite_alleles_k_theory.csv',
                    %w[theta e_k_n500 e_h_exact],
                    thetas.map { |th| [th, th * Math.log(1.0 + n_copies / th), th / (1.0 + th + th / n_copies.to_f)] })
      ctx.record_stage('S25', 'infinite_alleles_drift', Time.now - t25, 2)

      t26 = Time.now
      svg1 = '../../../outputs/figures/day009b_12_infinite_alleles.svg'
      SVG.write_xy(svg1, title: 'S25-S26: Infinite-Alleles Mutation-Drift Equilibrium (N=500)',
                         subtitle: "Expected heterozygosity E[H] vs theta/(1+theta) with O(1/N) correction (max diff=#{format('%.2f%%', 100 * max_rel)})",
                         xlabel: 'Mutation parameter theta', ylabel: 'Equilibrium heterozygosity E[H]',
                         series: [
                           { label: 'WF infinite-alleles (+-SE)', color: '#2563eb', x: thetas, y: res.map { |r| r[2] }, yerr: res.map { |r| r[3] }, line: false },
                           { label: 'Theory theta/(1+theta)', color: '#dc2626', x: thetas, y: res.map { |r| r[4] }, points: false },
                           { label: 'Exact O(1/N) Kimura-Crow', color: '#16a34a', x: thetas, y: res.map { |r| r[5] }, points: false, dashed: true }
                         ])
      ctx.register_svg('S26', svg1)
      kc_rows = thetas.map { |th| [th, th / (1.0 + th), th / (1.0 + th + th / 500.0), (th / (1.0 + th)) - (th / (1.0 + th + th / 500.0))] }
      ctx.write_csv('S26', 'out/s26_kimura_crow_corr.csv', %w[theta h_inf h_n500 o1n_diff], kc_rows)
      ctx.record_stage('S26', 'infinite_alleles_plot', Time.now - t26, 2)
    end

    def s27_s28_stationary_mut_sel(ctx)
      # c16: Stationary distribution of WF Markov chain with mutation + selection (N=20) vs long simulation
      t27 = Time.now
      n_copies = 20
      s = 0.08
      h = 0.5
      u = 0.04
      v = 0.04
      p_mat = Markov.wf_transition_matrix_float(n_copies, s: s, h: h, u: u, v: v)
      pi_exact = Markov.stationary_distribution(p_mat, rational: false)

      n_gens = ctx.small? ? 20_000 : 120_000
      n_blocks = 10
      blk_len = n_gens / n_blocks
      burn = 1_000
      rng = RNG.for_stage(27, 0)
      cnt = n_copies / 2
      burn.times { cnt = WF.step_biallelic(cnt, n_copies, rng, s: s, h: h, u: u, v: v) }

      blk_hist = Array.new(n_blocks) { Array.new(n_copies + 1, 0.0) }
      n_blocks.times do |b|
        blk_len.times do
          cnt = WF.step_biallelic(cnt, n_copies, rng, s: s, h: h, u: u, v: v)
          blk_hist[b][cnt] += 1.0 / blk_len
        end
      end

      sim_pi = Array.new(n_copies + 1, 0.0)
      se_pi = Array.new(n_copies + 1, 0.0)
      rows = []
      (0..n_copies).each do |k|
        vals = blk_hist.map { |bh| bh[k] }
        sim_pi[k] = Stats.mean(vals)
        se_pi[k] = Stats.se_iid(vals)
        rows << [k, k.to_f / n_copies, sim_pi[k], se_pi[k], pi_exact[k], (sim_pi[k] - pi_exact[k]).abs]
      end

      tv = 0.5 * (0..n_copies).sum(0.0) { |k| (sim_pi[k] - pi_exact[k]).abs }
      tv_se_sum = 0.5 * se_pi.sum(0.0)
      tv_thresh = [0.015, 3.0 * tv_se_sum].max
      ctx.metrics[:c16_tv] = tv
      ctx.metrics[:c16_thresh] = tv_thresh
      ctx.metrics[:c16_pass] = (tv <= tv_thresh)

      ctx.write_csv('S27', 'out/s27_stationary_mut_sel.csv',
                    %w[k p_freq pi_sim pi_se pi_markov_exact abs_diff], rows)
      ctx.write_csv('S27', 'out/s27_stationary_moments.csv',
                    %w[mean_k_sim mean_k_exact var_k_sim var_k_exact],
                    [[(0..n_copies).sum(0.0) { |k| k * sim_pi[k] },
                      (0..n_copies).sum(0.0) { |k| k * pi_exact[k] },
                      (0..n_copies).sum(0.0) { |k| k * k * sim_pi[k] } - ((0..n_copies).sum(0.0) { |k| k * sim_pi[k] })**2,
                      (0..n_copies).sum(0.0) { |k| k * k * pi_exact[k] } - ((0..n_copies).sum(0.0) { |k| k * pi_exact[k] })**2]])
      ctx.record_stage('S27', 'stationary_mut_sel', Time.now - t27, 2)

      t28 = Time.now
      ks_arr = (0..n_copies).to_a.map(&:to_f)
      svg1 = '../../../outputs/figures/day009b_13_stationary_dist.svg'
      SVG.write_hist(svg1, title: 'S27-S28: Stationary Distribution (WF Selection + Mutation, N=20)',
                           subtitle: "s=0.08, u=v=0.04, #{n_gens} gens: TV = #{format('%.4f', tv)} (threshold #{format('%.4f', tv_thresh)})",
                           xlabel: 'Allele A copy count k', ylabel: 'Stationary probability pi_k',
                           centers: ks_arr, counts: sim_pi,
                           theory_y: pi_exact, theory_label: 'Exact Markov pi')
      ctx.register_svg('S28', svg1)
      ctx.write_csv('S28', 'out/s28_stationary_summary.csv',
                    %w[n_copies s u v n_gens tv_dist tv_threshold],
                    [[n_copies, s, u, v, n_gens, tv, tv_thresh]])
      ctx.record_stage('S28', 'stationary_dist_plot', Time.now - t28, 2)
    end

    def s29_mut_sel_balance(ctx)
      t0 = Time.now
      # Purifying selection against deleterious recessive (h=0, s=-0.1) vs additive (h=0.5, s=-0.1):
      # Deterministic mutation-selection balance q* ~ sqrt(u / |s|) (recessive) vs u / (h|s|) (additive)
      u = 0.002
      s_del = -0.1
      n_copies = 4_000
      gens = ctx.small? ? 400 : 1_000
      rng = RNG.for_stage(29, 0)

      rows = []
      [[0.0, Math.sqrt(u / s_del.abs), 'recessive_h0'],
       [0.5, u / (0.5 * s_del.abs), 'additive_h05']].each do |h_val, q_theory, lbl|
        # Track deleterious allele q where w_AA=1, w_Aa=1+h*s_del, w_aa=1+s_del
        cnt = (n_copies * q_theory).round
        qs = []
        gens.times do |g|
          q = cnt.to_f / n_copies
          p = 1.0 - q
          w_bar = p * p + 2.0 * p * q * (1.0 + h_val * s_del) + q * q * (1.0 + s_del)
          q_sel = (p * q * (1.0 + h_val * s_del) + q * q * (1.0 + s_del)) / w_bar
          q_next = q_sel + (1.0 - q_sel) * u
          cnt = Binomial.sample(n_copies, q_next, rng)
          qs << (cnt.to_f / n_copies) if g >= gens / 2
        end
        ba = Stats.block_average(qs, 10)
        rows << [lbl, h_val, u, s_del, ba[:mean], ba[:se], q_theory]
      end
      ctx.write_csv('S29', 'out/s29_mut_sel_balance.csv',
                    %w[regime h u s q_mean q_se q_theory], rows)
      ctx.write_csv('S29', 'out/s29_mut_sel_ratio.csv',
                    %w[regime q_mean q_theory ratio_sim_theory],
                    rows.map { |r| [r[0], r[4], r[6], r[4] / r[6]] })
      ctx.record_stage('S29', 'mut_sel_balance', Time.now - t0, 2)
    end

    def s30_dominance_overdominance(ctx)
      t0 = Time.now
      # Overdominance: w_AA = 1 - s1, w_Aa = 1, w_aa = 1 - s2 -> stable equilibrium p* = s2 / (s1 + s2)
      s1 = 0.12
      s2 = 0.24
      p_eq = s2 / (s1 + s2) # 2/3
      n_copies = 1_000
      gens = 120
      rng = RNG.for_stage(30, 0)

      trajs = []
      rows = []
      [0.15, 0.90].each do |p0|
        cnt = (n_copies * p0).round
        tr = [cnt.to_f / n_copies]
        gens.times do
          p = cnt.to_f / n_copies
          q = 1.0 - p
          w_bar = p * p * (1.0 - s1) + 2.0 * p * q + q * q * (1.0 - s2)
          p_next = (p * p * (1.0 - s1) + p * q) / w_bar
          cnt = Binomial.sample(n_copies, p_next, rng)
          tr << (cnt.to_f / n_copies)
        end
        trajs << tr
      end
      (0..gens).step(4).each do |g|
        rows << [g, trajs[0][g], trajs[1][g], p_eq]
      end
      ctx.write_csv('S30', 'out/s30_overdominance.csv', %w[gen p_from_015 p_from_090 p_eq_analytic], rows)
      dom_rows = [0.0, 0.5, 1.0].map do |h_dom|
        p = 0.2
        20.times do
          q = 1.0 - p
          wb = p * p * 1.1 + 2.0 * p * q * (1.0 + 0.1 * h_dom) + q * q
          p = (p * p * 1.1 + p * q * (1.0 + 0.1 * h_dom)) / wb
        end
        [h_dom, 0.2, p]
      end
      ctx.write_csv('S30', 'out/s30_dominance_h_compare.csv', %w[h_dom p0 p_after_20g], dom_rows)

      svg1 = '../../../outputs/figures/day009b_14_overdominance.svg'
      SVG.write_xy(svg1, title: 'S30: Balancing Selection (Overdominance) Polymorphic Equilibrium',
                         subtitle: "s1=#{s1}, s2=#{s2}: both p0=0.15 and p0=0.90 converge to p* = s2/(s1+s2) = 0.667",
                         xlabel: 'Generation t', ylabel: 'Allele A frequency p_t',
                         series: [
                           { label: 'WF from p0=0.15', color: '#2563eb', x: rows.map { |r| r[0] }, y: rows.map { |r| r[1] }, points: false },
                           { label: 'WF from p0=0.90', color: '#16a34a', x: rows.map { |r| r[0] }, y: rows.map { |r| r[2] }, points: false }
                         ],
                         ref_y: p_eq, ref_label: 'Analytic p* = 2/3')
      ctx.register_svg('S30', svg1)
      ctx.record_stage('S30', 'dominance_overdominance', Time.now - t0, 2)
    end
  end
end
