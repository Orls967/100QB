# frozen_string_literal: true

require_relative 'rng'
require_relative 'binomial'
require_relative 'samplers'
require_relative 'markov'
require_relative 'stats'
require_relative 'svg'

module PopGen
  # Group A: Foundation Stages S01..S10
  module StagesA
    module_function

    def run(ctx)
      s01_rng_streams(ctx)
      s02_binomial_12_pairs(ctx)
      s03_binomial_branch_match(ctx)
      s04_multinomial(ctx)
      s05_poisson(ctx)
      s06_exponential_geom(ctx)
      s07_allele_fixtures(ctx)
      s08_freq_sum_invariants(ctx)
      s09_markov_rational_vs_float(ctx)
      s10_markov_absorption_checks(ctx)
    end

    def s01_rng_streams(ctx)
      t0 = Time.now
      n_draw = ctx.small? ? 5_000 : 20_000
      streams = (0...8).map do |st|
        rng = RNG.for_stage(1, st)
        Array.new(n_draw) { rng.rand }
      end

      rows = []
      8.times do |st|
        xs = streams[st]
        m = Stats.mean(xs)
        v = Stats.variance(xs)
        corr_next = Stats.pearson_corr(xs, streams[(st + 1) % 8])
        rows << [st, RNG.derive_seed(1, st).to_s(16), m, v, corr_next]
      end
      ctx.write_csv('S01', 'out/s01_rng_streams.csv',
                    %w[idx seed_hex u_mean u_var corr_with_next], rows)

      n_bins = 20
      counts = Array.new(n_bins, 0.0)
      streams[0].each do |u|
        b = (u * n_bins).floor
        b = n_bins - 1 if b >= n_bins
        counts[b] += 1.0 / (n_draw.to_f / n_bins)
      end
      centers = Array.new(n_bins) { |b| (b + 0.5) / n_bins }
      theory = Array.new(n_bins, 1.0)
      bin_rows = centers.each_with_index.map { |c, b| [b, c, counts[b], theory[b]] }
      ctx.write_csv('S01', 'out/s01_rng_bins.csv', %w[bin center density theory_density], bin_rows)

      svg1 = '../../../outputs/figures/day009b_01_rng_uniform.svg'
      SVG.write_hist(svg1, title: 'S01: SplitMix64 -> Random Uniformity Audit',
                           subtitle: "Derived stream 0 (N=#{n_draw}, 20 bins)",
                           xlabel: 'u in [0, 1)', ylabel: 'Normalized density',
                           centers: centers, counts: counts,
                           theory_y: theory, theory_label: 'Uniform = 1.0')
      ctx.register_svg('S01', svg1)
      ctx.record_stage('S01', 'rng_streams', Time.now - t0, 3)
    end

    def s02_binomial_12_pairs(ctx)
      t0 = Time.now
      # 12 (n, p) pairs covering both inversion (n*p_min < 10) and BTPE (n*p_min >= 10)
      pairs = [
        [8, 0.25], [15, 0.15], [30, 0.10], [50, 0.08], [100, 0.05], [200, 0.96],
        [40, 0.45], [80, 0.50], [150, 0.30], [300, 0.65], [600, 0.20], [1000, 0.50]
      ]
      n_rep = ctx.small? ? 5_000 : 15_000
      rng = RNG.for_stage(2, 0)

      rows = []
      min_p = 1.0
      all_pass = true
      pvals = []
      labels = []

      pairs.each_with_index do |(n, p), idx|
        r_min = [p, 1.0 - p].min
        branch = (n * r_min < 10.0) ? 'inv' : 'btpe'
        obs = Array.new(n + 1, 0)
        n_rep.times { obs[Binomial.sample(n, p, rng)] += 1 }
        exp = Binomial.pmf_array(n, p).map { |pr| pr * n_rep }
        gof = Stats.chi2_gof_pooled(obs, exp, min_exp: 5.0)
        pv = gof[:p_value]
        min_p = pv if pv < min_p
        all_pass = false if pv <= 0.001
        pvals << pv
        labels << "(#{n},#{p})"

        sample_mean = obs.each_with_index.sum(0.0) { |c, k| c * k } / n_rep
        rows << [idx, n, p, branch, n_rep, sample_mean, n * p, gof[:chi2], gof[:df], pv, pv > 0.001 ? 1 : 0]
      end

      ctx.metrics[:c18_min_p] = min_p
      ctx.metrics[:c18_all_pass] = all_pass
      ctx.write_csv('S02', 'out/s02_binomial_12pairs.csv',
                    %w[idx n p branch n_rep mean_sim mean_theory chi2 df p_value pass_bonf], rows)

      # Also export detailed PMF comparison for pair (40, 0.45)
      obs40 = Array.new(41, 0)
      n_rep.times { obs40[Binomial.sample(40, 0.45, rng)] += 1 }
      pmf_th = Binomial.pmf_array(40, 0.45)
      ks_arr = (0..40).to_a.map(&:to_f)
      emp_pmf = obs40.map { |c| c.to_f / n_rep }
      pmf_rows = (0..40).map { |k| [k, emp_pmf[k], pmf_th[k], emp_pmf[k] - pmf_th[k]] }
      ctx.write_csv('S02', 'out/s02_binomial_pmf40.csv', %w[k pmf_emp pmf_exact diff], pmf_rows)
      ctx.write_csv('S02', 'out/s02_binomial_moments.csv', %w[idx n p var_theory skew_theory],
                    pairs.each_with_index.map { |(n, p), idx| [idx, n, p, n * p * (1.0 - p), (1.0 - 2.0 * p) / Math.sqrt(n * p * (1.0 - p))] })

      svg1 = '../../../outputs/figures/day009b_02_binomial_pmf.svg'
      SVG.write_hist(svg1, title: 'S02: Exact Binomial Sampler Bin(n=40, p=0.45) [BTPE]',
                           subtitle: "N=#{n_rep} draws vs exact analytical PMF",
                           xlabel: 'Success count k', ylabel: 'Probability P(X = k)',
                           centers: ks_arr, counts: emp_pmf,
                           theory_y: pmf_th, theory_label: 'Exact Bin(40, 0.45)')
      ctx.register_svg('S02', svg1)

      svg2 = '../../../outputs/figures/day009b_03_binomial_pvals.svg'
      SVG.write_bars(svg2, title: 'S02: Chi-Square Goodness-of-Fit p-Values (12 Binomial Pairs)',
                           subtitle: 'Covers both inversion (np<10) and BTPE (np>=10); Bonferroni alpha=0.001',
                           ylabel: 'Chi-square p-value', labels: labels, vals: pvals)
      ctx.register_svg('S02', svg2)
      ctx.record_stage('S02', 'binomial_12_pairs', Time.now - t0, 4)
    end

    def s03_binomial_branch_match(ctx)
      t0 = Time.now
      # Compare inv vs btpe right at crossover n*p = 10..20
      cases = [[25, 0.4], [40, 0.3], [50, 0.35]]
      n_rep = ctx.small? ? 5_000 : 12_000
      rng = RNG.for_stage(3, 0)
      rows = []
      cases.each do |n, p|
        obs_inv = Array.new(n + 1, 0)
        obs_btp = Array.new(n + 1, 0)
        n_rep.times do
          obs_inv[Binomial.sample(n, p, rng, force_branch: :inv)] += 1
          obs_btp[Binomial.sample(n, p, rng, force_branch: :btpe)] += 1
        end
        exp = Binomial.pmf_array(n, p).map { |pr| pr * n_rep }
        g_inv = Stats.chi2_gof_pooled(obs_inv, exp)
        g_btp = Stats.chi2_gof_pooled(obs_btp, exp)
        tv = 0.5 * (0..n).sum(0.0) { |k| (obs_inv[k] - obs_btp[k]).abs } / n_rep
        rows << [n, p, g_inv[:p_value], g_btp[:p_value], tv]
      end
      ctx.write_csv('S03', 'out/s03_binomial_branches.csv',
                    %w[n p chi2_p_inv chi2_p_btpe tv_inv_vs_btpe], rows)
      ctx.record_stage('S03', 'binomial_branch_match', Time.now - t0, 1)
    end

    def s04_multinomial(ctx)
      t0 = Time.now
      probs = [0.15, 0.35, 0.30, 0.20]
      n_tot = 60
      n_rep = ctx.small? ? 4_000 : 10_000
      rng = RNG.for_stage(4, 0)
      sums = Array.new(4, 0.0)
      sq_sums = Array.new(4, 0.0)
      n_rep.times do
        c = Samplers.multinomial(n_tot, probs, rng)
        4.times do |i|
          sums[i] += c[i]
          sq_sums[i] += c[i] * c[i]
        end
      end
      rows = []
      4.times do |i|
        m = sums[i] / n_rep
        v = (sq_sums[i] / n_rep - m * m) * n_rep / (n_rep - 1)
        m_th = n_tot * probs[i]
        v_th = n_tot * probs[i] * (1.0 - probs[i])
        rows << [i, probs[i], m, m_th, v, v_th]
      end
      ctx.write_csv('S04', 'out/s04_multinomial.csv',
                    %w[bin prob mean_sim mean_theory var_sim var_theory], rows)
      cov_rows = (0...4).flat_map do |i|
        ((i + 1)...4).map { |j| [i, j, -n_tot * probs[i] * probs[j]] }
      end
      ctx.write_csv('S04', 'out/s04_multinomial_cov.csv', %w[bin_i bin_j cov_theory], cov_rows)
      ctx.record_stage('S04', 'multinomial_sampler', Time.now - t0, 2)
    end

    def s05_poisson(ctx)
      t0 = Time.now
      lams = [0.8, 4.5, 12.0, 35.0]
      n_rep = ctx.small? ? 5_000 : 12_000
      rng = RNG.for_stage(5, 0)
      rows = []
      lams.each do |lam|
        draws = Array.new(n_rep) { Samplers.poisson(lam, rng) }
        m = Stats.mean(draws.map(&:to_f))
        v = Stats.variance(draws.map(&:to_f))
        max_k = [draws.max, (lam + 6.0 * Math.sqrt(lam)).ceil].max
        obs = Array.new(max_k + 1, 0)
        draws.each { |k| obs[k] += 1 }
        exp = Array.new(max_k + 1) do |k|
          n_rep * Math.exp(-lam + k * Math.log(lam) - Math.lgamma(k + 1).first)
        end
        gof = Stats.chi2_gof_pooled(obs, exp)
        rows << [lam, lam < 10.0 ? 'knuth' : 'ptrs', m, v, gof[:p_value]]
      end
      ctx.write_csv('S05', 'out/s05_poisson.csv',
                    %w[theta branch mean_sim var_sim chi2_p], rows)
      pmf5_rows = (0..18).map do |k|
        [k, Math.exp(-4.5 + k * Math.log(4.5) - Math.lgamma(k + 1).first)]
      end
      ctx.write_csv('S05', 'out/s05_poisson_pmf45.csv', %w[k pmf_exact], pmf5_rows)
      ctx.record_stage('S05', 'poisson_sampler', Time.now - t0, 2)
    end

    def s06_exponential_geom(ctx)
      t0 = Time.now
      rates = [0.5, 1.5, 4.0]
      n_rep = ctx.small? ? 4_000 : 10_000
      rng = RNG.for_stage(6, 0)
      rows = []
      rates.each do |lam|
        xs = Array.new(n_rep) { Samplers.exponential(lam, rng) }
        ks = Stats.ks_against_cdf(xs, ->(x) { x <= 0.0 ? 0.0 : (1.0 - Math.exp(-lam * x)) })
        rows << ['exponential', lam, Stats.mean(xs), 1.0 / lam, Stats.variance(xs), 1.0 / (lam * lam), ks[:p_value]]
      end
      [0.2, 0.5, 0.8].each do |p|
        gs = Array.new(n_rep) { Samplers.geometric_trials(p, rng).to_f }
        rows << ['geometric', p, Stats.mean(gs), 1.0 / p, Stats.variance(gs), (1.0 - p) / (p * p), 1.0]
      end
      ctx.write_csv('S06', 'out/s06_exp_geom.csv',
                    %w[dist param mean_sim mean_theory var_sim var_theory ks_p], rows)
      ctx.record_stage('S06', 'exponential_geom', Time.now - t0, 1)
    end

    def s07_allele_fixtures(ctx)
      t0 = Time.now
      # Analytical 2-allele & 3-allele genotype/heterozygosity fixtures
      rows = []
      [[0.3, 0.7], [0.5, 0.5], [0.1, 0.9]].each_with_index do |ps, idx|
        h_exp = 1.0 - ps.sum(0.0) { |x| x * x }
        h_2pq = 2.0 * ps[0] * ps[1]
        rows << ["2allele_#{idx}", ps[0], ps[1], 0.0, h_exp, h_2pq, (h_exp - h_2pq).abs]
      end
      [[0.2, 0.3, 0.5], [1.0 / 3, 1.0 / 3, 1.0 / 3], [0.6, 0.3, 0.1]].each_with_index do |ps, idx|
        h_exp = 1.0 - ps.sum(0.0) { |x| x * x }
        h_pairs = 2.0 * (ps[0] * ps[1] + ps[0] * ps[2] + ps[1] * ps[2])
        rows << ["3allele_#{idx}", ps[0], ps[1], ps[2], h_exp, h_pairs, (h_exp - h_pairs).abs]
      end
      ctx.write_csv('S07', 'out/s07_allele_fixtures.csv',
                    %w[label p1 p2 p3 h_formula h_pairwise abs_err], rows)
      ctx.record_stage('S07', 'allele_fixtures', Time.now - t0, 1)
    end

    def s08_freq_sum_invariants(ctx)
      t0 = Time.now
      rng = RNG.for_stage(8, 0)
      n_copies = 200
      probs = [0.25, 0.35, 0.40]
      counts = probs.map { |pr| (n_copies * pr).round }
      rows = []
      100.times do |g|
        f_cur = counts.map { |c| c.to_f / n_copies }
        counts = Samplers.multinomial(n_copies, f_cur, rng)
        sum_f = counts.sum.to_f / n_copies
        rows << [g + 1, counts[0].to_f / n_copies, counts[1].to_f / n_copies,
                 counts[2].to_f / n_copies, (sum_f - 1.0).abs]
      end
      ctx.write_csv('S08', 'out/s08_freq_sum_invariants.csv',
                    %w[gen f1 f2 f3 sum_err], rows)
      ctx.record_stage('S08', 'freq_sum_invariants', Time.now - t0, 1)
    end

    def s09_markov_rational_vs_float(ctx)
      t0 = Time.now
      # Compare Rational vs Float WF Markov chain for N in {4, 6, 8, 10, 12}
      ns = [4, 6, 8, 10, 12]
      rows = []
      worst_err = 0.0
      ns.each do |n|
        p_rat = Markov.wf_transition_matrix_rational(n)
        p_flt = Markov.wf_transition_matrix_float(n)
        u_rat = Markov.fixation_probs(p_rat, rational: true).map(&:to_f)
        u_flt = Markov.fixation_probs(p_flt, rational: false)
        t_rat = Markov.mean_absorption_times(p_rat, rational: true).map(&:to_f)
        t_flt = Markov.mean_absorption_times(p_flt, rational: false)

        err_p = 0.0
        (0..n).each do |i|
          (0..n).each do |j|
            d = (p_rat[i][j].to_f - p_flt[i][j]).abs
            err_p = d if d > err_p
          end
        end
        err_u = (0..n).map { |i| (u_rat[i] - u_flt[i]).abs }.max
        err_t = (0..n).map { |i| (t_rat[i] - t_flt[i]).abs }.max
        m_err = [err_p, err_u, err_t].max
        worst_err = m_err if m_err > worst_err
        rows << [n, err_p, err_u, err_t, t_rat[n / 2]]
      end
      ctx.metrics[:c1_max_err] = worst_err
      ctx.write_csv('S09', 'out/s09_markov_rational_vs_float.csv',
                    %w[n max_p_err max_u_err max_t_abs_err t_abs_mid], rows)
      p12_r = Markov.wf_transition_matrix_rational(12)
      u12_r = Markov.fixation_probs(p12_r, rational: true)
      t12_r = Markov.mean_absorption_times(p12_r, rational: true)
      ctx.write_csv('S09', 'out/s09_markov_n12_rational.csv', %w[k p0 u_rational t_abs_rational],
                    (0..12).map { |i| [i, i / 12.0, u12_r[i].to_f, t12_r[i].to_f] })

      svg1 = '../../../outputs/figures/day009b_04_rational_vs_float.svg'
      SVG.write_xy(svg1, title: 'S09: Exact Rational vs Float WF Markov Chain (N=4..12)',
                         subtitle: "Max discrepancy <= #{format('%.2e', worst_err)} (threshold 1e-12)",
                         xlabel: 'Gene copies N', ylabel: 'Max |Rational - Float|',
                         series: [
                           { label: 'Fixation prob err', color: '#2563eb', x: ns, y: rows.map { |r| [r[2], 1e-16].max } },
                           { label: 'Absorption time err', color: '#dc2626', x: ns, y: rows.map { |r| [r[3], 1e-16].max } }
                         ],
                         log_y: true, ref_y: 1e-12, ref_label: 'c1 limit (1e-12)')
      ctx.register_svg('S09', svg1)
      ctx.record_stage('S09', 'markov_rational_vs_float', Time.now - t0, 2)
    end

    def s10_markov_absorption_checks(ctx)
      t0 = Time.now
      # N=8 Rational chain: neutral vs selection s=0.25 vs mutation u=v=0.05
      p_neu = Markov.wf_transition_matrix_rational(8)
      u_neu = Markov.fixation_probs(p_neu, rational: true)
      t_abs = Markov.mean_absorption_times(p_neu, rational: true)
      t_fix = Markov.mean_conditional_fixation_times(p_neu, rational: true)

      p_sel = Markov.wf_transition_matrix_rational(8, s: Rational(1, 4), h: Rational(1, 2))
      u_sel = Markov.fixation_probs(p_sel, rational: true)

      rows = (0..8).map do |i|
        [i, i / 8.0, u_neu[i].to_f, t_abs[i].to_f, t_fix[i].to_f, u_sel[i].to_f]
      end
      ctx.write_csv('S10', 'out/s10_markov_n8_exact.csv',
                    %w[k p0 u_neutral t_abs_neutral t_fix_cond_neutral u_sel_s025], rows)

      p_mut = Markov.wf_transition_matrix_rational(8, u: Rational(1, 20), v: Rational(1, 20))
      stat = Markov.stationary_distribution(p_mut, rational: true)
      s_rows = (0..8).map { |i| [i, stat[i].to_f, stat[i].numerator.to_s] }
      ctx.write_csv('S10', 'out/s10_stationary_n8_rational.csv', %w[k pi_exact num_digits_ref], s_rows)
      ctx.record_stage('S10', 'markov_absorption_checks', Time.now - t0, 2)
    end
  end
end
