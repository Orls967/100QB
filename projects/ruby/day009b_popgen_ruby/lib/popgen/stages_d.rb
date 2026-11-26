# frozen_string_literal: true

require_relative 'rng'
require_relative 'coalescent'
require_relative 'wf'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group D: Kingman Coalescent, SFS, Tajima's D, Ewens Sampling Formula & WF Comparison (S31..S42)
  module StagesD
    module_function

    def run(ctx)
      s31_s34_kingman_moments(ctx)
      s35_sfs_unfolded_folded(ctx)
      s36_s37_tajima_d(ctx)
      s38_s_poisson_mixture(ctx)
      s39_s40_ewens_formula(ctx)
      s41_coalescent_vs_wf(ctx)
      s42_exp_growth_coalescent(ctx)
    end

    def s31_s34_kingman_moments(ctx)
      t31 = Time.now
      ns = [5, 20, 100]
      theta = 4.0
      n_rep = ctx.small? ? 1_500 : 6_000

      rows = []
      max_z_c6 = 0.0
      sfs_20 = nil

      ns.each_with_index do |n, idx|
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(31 + idx, b_i)
          b_ids.map do
            tr = Coalescent.simulate_tree(n, theta, rng)
            [tr[:t_mrca], tr[:total_length], tr[:s].to_f, n == 20 ? tr[:sfs] : nil]
          end
        end
        flat = res.flatten(1)
        tm = flat.map { |x| x[0] }
        tl = flat.map { |x| x[1] }
        ss = flat.map { |x| x[2] }
        sfs_20 = flat.map { |x| x[3] } if n == 20

        a1 = Coalescent.harmonic(n - 1)
        a2 = Coalescent.harmonic_sq(n - 1)
        e_tm = 2.0 * (1.0 - 1.0 / n)
        e_l = 2.0 * a1
        e_s = theta * a1
        v_s = theta * a1 + theta * theta * a2

        m_tm = Stats.mean(tm)
        se_tm = Stats.se_iid(tm)
        m_l = Stats.mean(tl)
        se_l = Stats.se_iid(tl)
        m_s = Stats.mean(ss)
        se_s = Stats.se_iid(ss)

        # Delete-10-block jackknife SE for Var(S)
        nb = 10
        bsz = ss.length / nb
        jk_v = Array.new(nb) do |b|
          sub = ss[0...(b * bsz)] + ss[((b + 1) * bsz)..]
          Stats.variance(sub)
        end
        v_s_sim = Stats.variance(ss)
        se_vs = Math.sqrt(((nb - 1.0) / nb) * jk_v.sum(0.0) { |x| (x - Stats.mean(jk_v))**2 })

        zs = [
          (m_tm - e_tm).abs / se_tm,
          (m_l - e_l).abs / se_l,
          (m_s - e_s).abs / se_s,
          (v_s_sim - v_s).abs / se_vs
        ]
        max_z_c6 = zs.max if zs.max > max_z_c6
        rows << [n, theta, n_rep, m_tm, se_tm, e_tm, m_l, se_l, e_l, m_s, se_s, e_s, v_s_sim, se_vs, v_s, zs.max]
      end

      ctx.metrics[:c6_max_z] = max_z_c6
      ctx.metrics[:c6_pass] = (max_z_c6 <= 3.0)
      ctx.metrics[:sfs_20_raw] = sfs_20

      ctx.write_csv('S31', 'out/s31_kingman_moments.csv',
                    %w[n theta n_rep tmrca_sim tmrca_se tmrca_theory len_sim len_se len_theory s_sim s_se s_theory var_s_sim var_s_se var_s_theory max_z], rows)
      ctx.write_csv('S31', 'out/s31_kingman_harmonics.csv',
                    %w[n a1 a2 e_tmrca e_length],
                    [2, 5, 10, 20, 50, 100].map { |nk| [nk, Coalescent.harmonic(nk - 1), Coalescent.harmonic_sq(nk - 1), 2.0 * (1.0 - 1.0 / nk), 2.0 * Coalescent.harmonic(nk - 1)] })
      ctx.record_stage('S31', 'kingman_tmrca_len', Time.now - t31, 2)

      # S32..S34 write derived tables and SVG plots for E[T_MRCA], E[L], E[S], Var(S)
      t32 = Time.now
      ctx.write_csv('S32', 'out/s32_kingman_tmrca.csv',
                    %w[n tmrca_sim tmrca_se tmrca_theory z_score],
                    rows.map { |r| [r[0], r[3], r[4], r[5], (r[3] - r[5]).abs / r[4]] })
      svg32 = '../../../outputs/figures/day009b_15_kingman_tmrca.svg'
      SVG.write_xy(svg32, title: 'S31-S32: Kingman Coalescent E[T_MRCA] = 2(1 - 1/n)',
                          subtitle: "n in {5, 20, 100}, R=#{n_rep} trees",
                          xlabel: 'Sample size n', ylabel: 'Mean T_MRCA (2N gens)',
                          series: [
                            { label: 'Coalescent sim (+-SE)', color: '#2563eb', x: ns, y: rows.map { |r| r[3] }, yerr: rows.map { |r| r[4] }, line: false },
                            { label: 'Exact 2(1 - 1/n)', color: '#dc2626', x: ns, y: rows.map { |r| r[5] } }
                          ], log_x: true)
      ctx.register_svg('S32', svg32)
      ctx.record_stage('S32', 'kingman_tmrca_plot', Time.now - t32, 2)

      t33 = Time.now
      ctx.write_csv('S33', 'out/s33_kingman_length_s.csv',
                    %w[n len_sim len_se len_theory s_sim s_se s_theory],
                    rows.map { |r| [r[0], r[6], r[7], r[8], r[9], r[10], r[11]] })
      svg33 = '../../../outputs/figures/day009b_16_kingman_len.svg'
      SVG.write_xy(svg33, title: 'S33: Total Tree Length E[L] = 2 a_n & Segregating Sites E[S]',
                          subtitle: "theta = #{theta}: logarithmic growth with sample size n",
                          xlabel: 'Sample size n', ylabel: 'Mean E[L] / E[S]',
                          series: [
                            { label: 'E[L] sim', color: '#2563eb', x: ns, y: rows.map { |r| r[6] }, yerr: rows.map { |r| r[7] }, line: false },
                            { label: 'E[L] theory', color: '#2563eb', x: ns, y: rows.map { |r| r[8] }, points: false, dashed: true },
                            { label: 'E[S] sim', color: '#16a34a', x: ns, y: rows.map { |r| r[9] }, yerr: rows.map { |r| r[10] }, line: false },
                            { label: 'E[S] theory', color: '#16a34a', x: ns, y: rows.map { |r| r[11] }, points: false, dashed: true }
                          ], log_x: true)
      ctx.register_svg('S33', svg33)
      ctx.record_stage('S33', 'kingman_length_s', Time.now - t33, 2)

      t34 = Time.now
      ctx.write_csv('S34', 'out/s34_kingman_var_s.csv',
                    %w[n var_s_sim var_s_se var_s_theory z_score],
                    rows.map { |r| [r[0], r[12], r[13], r[14], (r[12] - r[14]).abs / r[13]] })
      svg34 = '../../../outputs/figures/day009b_17_kingman_vars.svg'
      SVG.write_xy(svg34, title: 'S34: Variance of Segregating Sites Var(S) = theta a_n + theta^2 b_n',
                          subtitle: "Coalescent tree stochasticity creates overdispersion vs Poisson",
                          xlabel: 'Sample size n', ylabel: 'Var(S)',
                          series: [
                            { label: 'Sim Var(S) (+-JK SE)', color: '#2563eb', x: ns, y: rows.map { |r| r[12] }, yerr: rows.map { |r| r[13] }, line: false },
                            { label: 'Exact coalescent Var(S)', color: '#dc2626', x: ns, y: rows.map { |r| r[14] } }
                          ], log_x: true)
      ctx.register_svg('S34', svg34)
      ctx.record_stage('S34', 'kingman_var_s', Time.now - t34, 2)
    end

    def s35_sfs_unfolded_folded(ctx)
      t0 = Time.now
      theta = 4.0
      sfs_list = ctx.metrics[:sfs_20_raw]
      n_rep = sfs_list.length

      rows = []
      max_z = 0.0
      obs_tot = Array.new(19, 0.0)
      exp_tot = Array.new(19, 0.0)

      (1..19).each do |i|
        vals = sfs_list.map { |arr| arr[i - 1].to_f }
        m = Stats.mean(vals)
        se = Stats.se_iid(vals)
        th = theta / i.to_f
        z = se > 1e-12 ? (m - th).abs / se : 0.0
        max_z = z if z > max_z
        obs_tot[i - 1] = vals.sum
        exp_tot[i - 1] = n_rep * th
        rows << [i, m, se, th, z]
      end

      # Proper coalescent-variance Wald chi-square across the 19 SFS classes using empirical SEs
      # (naive multinomial chi2 ignores branch-length overdispersion across trees)
      chi2_wald = rows.sum(0.0) { |r| r[4] * r[4] }
      p_wald = Stats.chi2_p_value(chi2_wald, rows.length)
      ctx.metrics[:c7_max_z] = max_z
      ctx.metrics[:c7_chi2_p] = p_wald
      ctx.metrics[:c7_pass] = (max_z <= 3.0 && p_wald > 0.01)

      ctx.write_csv('S35', 'out/s35_sfs_n20.csv', %w[i xi_mean xi_se xi_theory z_score], rows)
      ctx.write_csv('S35', 'out/s35_sfs_gof.csv',
                    %w[n theta chi2_wald df p_wald max_class_z],
                    [[20, theta, chi2_wald, rows.length, p_wald, max_z]])

      # Folded SFS
      fold_rows = (1..10).map do |i|
        if i == 10
          [i, rows[9][1], rows[9][3]]
        else
          [i, rows[i - 1][1] + rows[19 - i][1], theta * (1.0 / i + 1.0 / (20 - i))]
        end
      end
      ctx.write_csv('S35', 'out/s35_sfs_folded.csv', %w[minor_i eta_mean eta_theory], fold_rows)

      svg1 = '../../../outputs/figures/day009b_18_sfs_n20.svg'
      SVG.write_xy(svg1, title: 'S35: Unfolded Site Frequency Spectrum (n=20, theta=4.0)',
                         subtitle: "max |diff|/SE=#{format('%.2f', max_z)}, Wald chi2 p=#{format('%.4f', p_wald)}",
                         xlabel: 'Derived allele frequency class i', ylabel: 'Expected count E[xi_i]',
                         series: [
                           { label: 'Kingman SFS (+-SE)', color: '#2563eb', x: rows.map { |r| r[0] }, y: rows.map { |r| r[1] }, yerr: rows.map { |r| r[2] }, line: false },
                           { label: 'Neutral theory theta / i', color: '#dc2626', x: rows.map { |r| r[0] }, y: rows.map { |r| r[3] } }
                         ])
      ctx.register_svg('S35', svg1)
      ctx.record_stage('S35', 'sfs_unfolded_folded', Time.now - t0, 4)
    end

    def s36_s37_tajima_d(ctx)
      t36 = Time.now
      n = 30
      thetas = [2.0, 5.0, 10.0, 20.0]
      n_rep = ctx.small? ? 2_000 : 8_000
      tc = Coalescent.tajima_constants(n)
      rows = []
      unnorm_rows = []

      thetas.each_with_index do |th, idx|
        batches = (0...n_rep).each_slice((n_rep.to_f / ctx.workers).ceil).to_a
        res = RNG.parallel_map(batches, workers: ctx.workers) do |b_ids, b_i|
          rng = RNG.for_stage(36 + idx, b_i)
          b_ids.map do
            tr = Coalescent.simulate_tree(n, th, rng)
            s = tr[:s]
            [tr[:tajima_d], tr[:theta_w], tr[:pi], tc[:e1] * s + tc[:e2] * s * (s - 1)]
          end
        end
        flat = res.flatten(1)
        ds = flat.map { |x| x[0] }
        tws = flat.map { |x| x[1] }
        pis = flat.map { |x| x[2] }
        diffs = flat.map { |x| x[2] - x[1] }
        evars = flat.map { |x| x[3] }
        md = Stats.mean(ds)
        sed = Stats.se_iid(ds)
        vd = Stats.variance(ds)
        var_diff = Stats.variance(diffs)
        mean_evar = Stats.mean(evars)
        rows << [th, n, md, sed, (md / sed).abs, vd, 1.0, (vd - 1.0).abs, Stats.mean(tws), Stats.mean(pis)]
        unnorm_rows << [th, n, Stats.mean(diffs), Stats.se_iid(diffs), var_diff, mean_evar, var_diff / mean_evar]
      end

      # Primary theta = 10.0 (idx 2, well away from S=0 boundary saturation) for c8 check
      primary = rows[2]
      ctx.metrics[:c8_mean_d] = primary[2]
      ctx.metrics[:c8_se_d] = primary[3]
      ctx.metrics[:c8_z_mean] = primary[4]
      ctx.metrics[:c8_var_d] = primary[5]
      ctx.metrics[:c8_pass] = (primary[4] <= 3.0 && (primary[5] - 1.0).abs <= 0.10)

      ctx.write_csv('S36', 'out/s36_tajima_d_summary.csv',
                    %w[theta n d_mean d_se z_mean_0 var_d var_tajima_lit abs_var_diff theta_w_mean pi_mean], rows)
      ctx.write_csv('S36', 'out/s36_tajima_unnorm_var.csv',
                    %w[theta n diff_mean diff_se var_diff e_tajima_denom2 var_ratio], unnorm_rows)
      ctx.record_stage('S36', 'tajima_d_moments', Time.now - t36, 2)

      t37 = Time.now
      rng = RNG.for_stage(37, 0)
      ds_h = Array.new(ctx.small? ? 2_000 : 6_000) { Coalescent.simulate_tree(n, 10.0, rng)[:tajima_d] }
      nb = 24
      dmin = -3.0
      dmax = 3.0
      dd = (dmax - dmin) / nb
      hist = Array.new(nb, 0.0)
      ds_h.each do |v|
        b = ((v - dmin) / dd).floor
        hist[b] += 1.0 / (ds_h.length * dd) if b >= 0 && b < nb
      end
      centers = Array.new(nb) { |b| dmin + (b + 0.5) * dd }
      gauss = centers.map { |x| Math.exp(-0.5 * x * x) / Math.sqrt(2.0 * Math::PI) }
      ctx.write_csv('S37', 'out/s37_tajima_d_hist.csv',
                    %w[bin d_center pdf_tajima std_normal_ref],
                    centers.each_with_index.map { |c, b| [b, c, hist[b], gauss[b]] })

      svg1 = '../../../outputs/figures/day009b_19_tajima_d.svg'
      SVG.write_hist(svg1, title: "S36-S37: Tajima's D Distribution (n=30, theta=10.0)",
                           subtitle: "Mean=#{format('%.4f', primary[2])} +- #{format('%.4f', primary[3])}, Var=#{format('%.4f', primary[5])}",
                           xlabel: "Tajima's D", ylabel: 'Probability density',
                           centers: centers, counts: hist,
                           theory_y: gauss, theory_label: 'N(0, 1) reference')
      ctx.register_svg('S37', svg1)
      ctx.record_stage('S37', 'tajima_d_hist', Time.now - t37, 2)
    end

    def s38_s_poisson_mixture(ctx)
      t0 = Time.now
      # Compare Kingman S distribution (n=10, theta=3.0) vs exact Tavaré analytical formula for P(S = k):
      # Actually for n=2: P(S = k) = (theta / (theta + 1))^k * (1 / (theta + 1)) (Geometric!).
      # For n=10: compare coalescent S variance (overdispersed) against naive Poisson(E[S]).
      n = 10
      theta = 3.0
      n_rep = ctx.small? ? 3_000 : 10_000
      rng = RNG.for_stage(38, 0)
      ss = Array.new(n_rep) { Coalescent.simulate_tree(n, theta, rng)[:s] }
      e_s = theta * Coalescent.harmonic(n - 1)
      kmax = 28
      obs = Array.new(kmax + 1, 0.0)
      pois = Array.new(kmax + 1, 0.0)
      ss.each { |s| obs[s] += 1.0 / n_rep if s <= kmax }
      (0..kmax).each do |k|
        pois[k] = Math.exp(-e_s + k * Math.log(e_s) - Math.lgamma(k + 1).first)
      end
      rows = (0..kmax).map { |k| [k, obs[k], pois[k]] }
      ctx.write_csv('S38', 'out/s38_s_poisson_mixture.csv', %w[s_muts pmf_kingman pmf_naive_poisson], rows)
      ctx.write_csv('S38', 'out/s38_s_moments.csv',
                    %w[mean_s_sim mean_s_theory var_s_sim var_s_kingman var_s_naive_poisson],
                    [[Stats.mean(ss.map(&:to_f)), e_s, Stats.variance(ss.map(&:to_f)),
                      e_s + theta * theta * Coalescent.harmonic_sq(n - 1), e_s]])

      svg1 = '../../../outputs/figures/day009b_20_s_mixture.svg'
      SVG.write_hist(svg1, title: 'S38: Distribution of Segregating Sites S (n=10, theta=3.0)',
                           subtitle: 'Coalescent tree variability yields wider tail than naive Poisson(E[S])',
                           xlabel: 'Segregating sites S', ylabel: 'Probability P(S = k)',
                           centers: (0..kmax).to_a.map(&:to_f), counts: obs,
                           theory_y: pois, theory_label: 'Naive Poisson(E[S])')
      ctx.register_svg('S38', svg1)
      ctx.record_stage('S38', 's_poisson_mixture', Time.now - t0, 2)
    end

    def s39_s40_ewens_formula(ctx)
      t39 = Time.now
      # Exact Ewens partition test for n = 6, theta = 1.8
      n = 6
      theta = 1.8
      parts = Coalescent.all_partitions(n)
      exact_probs = parts.map { |pt| Coalescent.ewens_partition_prob(pt, n, theta) }
      n_rep = ctx.small? ? 6_000 : 24_000
      rng = RNG.for_stage(39, 0)
      counts = Hash.new(0)
      n_rep.times do
        pt = Coalescent.ewens_sample_urn(n, theta, rng)[:partition]
        counts[pt] += 1
      end
      obs = parts.map { |pt| counts[pt] }
      exp = exact_probs.map { |pr| pr * n_rep }
      gof = Stats.chi2_gof_pooled(obs, exp, min_exp: 5.0)
      rows = parts.each_with_index.map do |pt, idx|
        [idx, pt.join('-'), pt.sum, obs[idx], exp[idx], obs[idx].to_f / n_rep, exact_probs[idx]]
      end
      ctx.metrics[:c9_ewens_p] = gof[:p_value]
      ctx.write_csv('S39', 'out/s39_ewens_n6_partitions.csv',
                    %w[part_id partition k_n obs_count exp_count p_sim p_exact], rows)
      ctx.write_csv('S39', 'out/s39_ewens_gof.csv',
                    %w[n theta n_partitions n_rep chi2 df p_value],
                    [[n, theta, parts.length, n_rep, gof[:chi2], gof[:df], gof[:p_value]]])
      ctx.record_stage('S39', 'ewens_n6_partitions', Time.now - t39, 2)

      # S40: E[K_n] across n in {10, 25, 50, 100}
      t40 = Time.now
      ns = [10, 25, 50, 100]
      th_k = 2.5
      n_rep_k = ctx.small? ? 2_000 : 6_000
      k_rows = []
      max_z_k = 0.0
      ns.each_with_index do |nk, idx|
        rng_k = RNG.for_stage(40, idx)
        ks = Array.new(n_rep_k) { Coalescent.ewens_sample_urn(nk, th_k, rng_k)[:k_n].to_f }
        m = Stats.mean(ks)
        se = Stats.se_iid(ks)
        ek = Coalescent.ewens_expected_k(nk, th_k)
        z = (m - ek).abs / se
        max_z_k = z if z > max_z_k
        k_rows << [nk, th_k, m, se, ek, z]
      end
      ctx.metrics[:c9_kn_max_z] = max_z_k
      ctx.metrics[:c9_pass] = (max_z_k <= 3.0 && gof[:p_value] > 0.01 && ctx.metrics[:c9_h_max_rel] <= 0.05)

      ctx.write_csv('S40', 'out/s40_ewens_expected_kn.csv',
                    %w[n theta k_mean k_se k_exact z_score], k_rows)

      svg1 = '../../../outputs/figures/day009b_21_ewens_n6.svg'
      SVG.write_bars(svg1, title: 'S39-S40: Ewens Sampling Partition Probabilities (n=6, theta=1.8)',
                           subtitle: "11 partitions: Hoppe urn vs exact Ewens formula (chi2 p=#{format('%.4f', gof[:p_value])})",
                           ylabel: 'Partition probability',
                           labels: parts.map(&:sum).each_with_index.map { |k, i| "P#{i}(K=#{k})" },
                           vals: obs.map { |o| o.to_f / n_rep })
      ctx.register_svg('S40', svg1)
      ctx.record_stage('S40', 'ewens_expected_kn', Time.now - t40, 2)
    end

    def s41_coalescent_vs_wf(ctx)
      t0 = Time.now
      # Compare Kingman coalescent against direct discrete-generation WF genealogy at moderate N=200
      n_copies = 200
      n_sample = 12
      theta = 3.0
      n_rep = ctx.small? ? 800 : 3_000
      rng_c = RNG.for_stage(41, 1)
      rng_w = RNG.for_stage(41, 2)

      coal_res = Array.new(n_rep) { Coalescent.simulate_tree(n_sample, theta, rng_c) }
      wf_res = Array.new(n_rep) { WF.simulate_wf_genealogy_sample(n_copies, n_sample, theta, 10_000, rng_w) }

      rows = [
        ['t_mrca', Stats.mean(coal_res.map { |r| r[:t_mrca] }), Stats.se_iid(coal_res.map { |r| r[:t_mrca] }),
         Stats.mean(wf_res.map { |r| r[:t_mrca] }), Stats.se_iid(wf_res.map { |r| r[:t_mrca] })],
        ['total_length', Stats.mean(coal_res.map { |r| r[:total_length] }), Stats.se_iid(coal_res.map { |r| r[:total_length] }),
         Stats.mean(wf_res.map { |r| r[:total_length] }), Stats.se_iid(wf_res.map { |r| r[:total_length] })],
        ['s_muts', Stats.mean(coal_res.map { |r| r[:s].to_f }), Stats.se_iid(coal_res.map { |r| r[:s].to_f }),
         Stats.mean(wf_res.map { |r| r[:s].to_f }), Stats.se_iid(wf_res.map { |r| r[:s].to_f })]
      ].map do |lbl, mc, sec, mw, sew|
        z = (mc - mw).abs / Math.sqrt(sec * sec + sew * sew)
        [lbl, mc, sec, mw, sew, z]
      end

      ctx.write_csv('S41', 'out/s41_coalescent_vs_wf.csv',
                    %w[metric coal_mean coal_se wf_n200_mean wf_se z_diff], rows)
      ctx.write_csv('S41', 'out/s41_wf_timescale.csv',
                    %w[n_copies n_sample theta n_rep],
                    [[n_copies, n_sample, theta, n_rep]])
      ctx.record_stage('S41', 'coalescent_vs_wf', Time.now - t0, 2)
    end

    def s42_exp_growth_coalescent(ctx)
      t0 = Time.now
      # Exponential population growth (beta in {0.0, 2.0, 6.0, 15.0}) shifts Tajima's D negative (star-like tree)!
      betas = [0.0, 2.0, 6.0, 15.0]
      n = 30
      theta = 8.0
      n_rep = ctx.small? ? 800 : 2_500
      rows = []
      betas.each_with_index do |beta, idx|
        rng = RNG.for_stage(42, idx)
        trees = Array.new(n_rep) { Coalescent.simulate_tree(n, theta, rng, growth_beta: beta) }
        ds = trees.map { |tr| tr[:tajima_d] }
        tms = trees.map { |tr| tr[:t_mrca] }
        rows << [beta, Stats.mean(tms), Stats.se_iid(tms), Stats.mean(ds), Stats.se_iid(ds)]
      end
      ctx.write_csv('S42', 'out/s42_exp_growth_coalescent.csv',
                    %w[beta t_mrca_mean t_mrca_se tajima_d_mean tajima_d_se], rows)

      svg1 = '../../../outputs/figures/day009b_22_exp_growth.svg'
      SVG.write_xy(svg1, title: 'S42: Exponential Growth Coalescent Skewing Tajima D Negative',
                         subtitle: 'Star-shaped genealogies produce excess singletons and negative Tajima D',
                         xlabel: 'Growth rate beta', ylabel: "Mean Tajima's D",
                         series: [
                           { label: "Mean Tajima's D (+-SE)", color: '#2563eb', x: betas, y: rows.map { |r| r[3] }, yerr: rows.map { |r| r[4] } }
                         ], ref_y: 0.0, ref_label: 'Constant N (D=0)')
      ctx.register_svg('S42', svg1)
      ctx.record_stage('S42', 'exp_growth_coalescent', Time.now - t0, 2)
    end
  end
end
