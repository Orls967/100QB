# frozen_string_literal: true

require_relative 'samplers'
require_relative 'coalescent'

module PopGen
  # Selective sweep (hitchhiking) and two-locus linkage disequilibrium (D, D', r^2).
  module LDSweep
    module_function

    # Compute D, D_prime, and r^2 from 4 haplotype frequencies [f_AB, f_Ab, f_aB, f_ab]
    def ld_stats(freqs)
      fab, fa_b, f_ab, _f_a_b = freqs
      pa = fab + fa_b
      pb = fab + f_ab
      d = fab - pa * pb
      denom_r2 = pa * (1.0 - pa) * pb * (1.0 - pb)
      r2 = denom_r2 > 1e-15 ? (d * d / denom_r2) : 0.0
      dmax = if d >= 0.0
               [pa * (1.0 - pb), (1.0 - pa) * pb].min
             else
               [pa * pb, (1.0 - pa) * (1.0 - pb)].min
             end
      d_prime = dmax > 1e-15 ? (d / dmax) : 0.0
      { d: d, d_prime: d_prime, r2: r2, pa: pa, pb: pb }
    end

    # Deterministic two-locus recursion D_t = D_0 * (1 - r)^t
    def deterministic_decay(d0, r, gens)
      Array.new(gens + 1) { |t| d0 * ((1.0 - r)**t) }
    end

    # Finite-N two-locus Wright-Fisher simulation with recombination r
    # Initial haplotype counts [n_AB, n_Ab, n_aB, n_ab] summing to n_copies
    def simulate_ld_wf(n_copies, init_counts, r, gens, rng, sample_times: nil)
      counts = init_counts.dup
      res = {}
      res[0] = ld_stats(counts.map { |c| c.to_f / n_copies }) if sample_times.nil? || sample_times.include?(0)

      (1..gens).each do |g|
        freqs = counts.map { |c| c.to_f / n_copies }
        fab, fa_b, f_ab, f_a_b = freqs
        pa = fab + fa_b
        pb = fab + f_ab
        d = fab - pa * pb
        # Next generation expected haplotype frequencies after random mating + recombination r
        next_freqs = [
          fab - r * d,
          fa_b + r * d,
          f_ab + r * d,
          f_a_b - r * d
        ]
        # Clamp tiny rounding noise
        next_freqs = next_freqs.map { |v| v < 0.0 ? 0.0 : v }
        sum_f = next_freqs.sum(0.0)
        next_freqs = next_freqs.map { |v| v / sum_f }
        counts = Samplers.multinomial(n_copies, next_freqs, rng)
        if sample_times.nil? || sample_times.include?(g)
          res[g] = ld_stats(counts.map { |c| c.to_f / n_copies })
        end
      end
      res
    end

    # Selective sweep + hitchhiking simulator (structured coalescent conditioned on sweep trajectory):
    # Forward-conditioned beneficial allele trajectory x(t) from 1/N -> 1 under selection s (WF diploid/haploid),
    # then trace n_sample neutral lineages backwards through the sweep at recombination distance r!
    # Also supports neutral control (s = 0, no sweep conditioning -> pure Kingman!).
    def generate_sweep_trajectory(n_copies, s, rng)
      loop do
        traj = [1.0 / n_copies]
        cnt = 1
        while cnt.positive? && cnt < n_copies
          p = cnt.to_f / n_copies
          p_sel = (p * (1.0 + s)) / (1.0 + s * p)
          cnt = Binomial.sample(n_copies, p_sel, rng)
          traj << (cnt.to_f / n_copies)
        end
        return traj if cnt == n_copies
      end
    end

    # Trace n_sample neutral lineages backwards through a fixed sweep trajectory (from t = T_fix down to 0)
    # In each generation backwards with beneficial frequency x:
    # - Lineages in B background (size N*x) coalesce at rate 1/(N*x) or recombine to b at rate r*(1-x)
    # - Lineages in b background (size N*(1-x)) coalesce at rate 1/(N*(1-x)) or recombine to B at rate r*x
    # Once sweep reaches origin (t = 0), remaining lineages coalesce via standard Kingman (N copies).
    def simulate_hitchhiking_diversity(n_copies, n_sample, s, r, theta, rng, sweep_traj: nil)
      if s <= 0.0
        # Neutral control (no sweep)
        res = Coalescent.simulate_tree(n_sample, theta, rng)
        return { pi: res[:pi], s: res[:s], t_mrca: res[:t_mrca] }
      end

      traj = sweep_traj || generate_sweep_trajectory(n_copies, s, rng)
      # All samples taken immediately at fixation => all start on beneficial B background (bg = 1)
      lineages = Array.new(n_sample) { |i| { bg: 1, leaves: [i] } }
      t_pair = Array.new(n_sample) { Array.new(n_sample, 0.0) }
      total_len = 0.0
      t_back = 0

      (traj.length - 1).downto(1) do |idx|
        break if lineages.length <= 1

        x = traj[idx]
        t_back += 1
        total_len += lineages.length.to_f / n_copies

        # 1. Recombination events between B and b backgrounds
        if r > 0.0
          lineages.each do |lin|
            if lin[:bg] == 1
              lin[:bg] = 0 if rng.rand < r * (1.0 - x)
            else
              lin[:bg] = 1 if rng.rand < r * x
            end
          end
        end

        # 2. Coalescence within B background (N*x copies) and within b background (N*(1-x) copies)
        [1, 0].each do |bg|
          sub = lineages.select { |lin| lin[:bg] == bg }
          next if sub.length <= 1

          n_bg = (bg == 1 ? (n_copies * x) : (n_copies * (1.0 - x))).round
          n_bg = 1 if n_bg < 1
          other = lineages.reject { |lin| lin[:bg] == bg }
          groups = Hash.new { |h, k| h[k] = [] }
          sub.each do |lin|
            parent = rng.rand(n_bg)
            groups[parent] << lin
          end
          merged_sub = []
          groups.each_value do |g_list|
            if g_list.length == 1
              merged_sub << g_list[0]
            else
              acc_leaves = g_list[0][:leaves].dup
              g_list[1..].each do |other_lin|
                acc_leaves.each do |la|
                  other_lin[:leaves].each do |lb|
                    tv = t_back.to_f / n_copies
                    t_pair[la][lb] = tv
                    t_pair[lb][la] = tv
                  end
                end
                acc_leaves.concat(other_lin[:leaves])
              end
              merged_sub << { bg: bg, leaves: acc_leaves }
            end
          end
          lineages = other + merged_sub
        end
      end

      # Standard Kingman phase before sweep origin (in coalescent units)
      t_coal = t_back.to_f / n_copies
      k = lineages.length
      while k > 1
        dt = Samplers.exponential(0.5 * k * (k - 1), rng)
        t_coal += dt
        total_len += k * dt
        i = rng.rand(k)
        j = rng.rand(k - 1)
        j += 1 if j >= i
        i, j = j, i if i > j
        lin_b = lineages.delete_at(j)
        lin_a = lineages.delete_at(i)
        lin_a[:leaves].each do |la|
          lin_b[:leaves].each do |lb|
            t_pair[la][lb] = t_coal
            t_pair[lb][la] = t_coal
          end
        end
        lineages << { bg: 0, leaves: lin_a[:leaves] + lin_b[:leaves] }
        k -= 1
      end

      # Compute pairwise pi = theta * mean(T_ij) and S ~ Poisson(0.5 * theta * total_len)
      sum_t = 0.0
      n_pairs = 0.5 * n_sample * (n_sample - 1)
      (0...n_sample).each do |a|
        ((a + 1)...n_sample).each { |b| sum_t += t_pair[a][b] }
      end
      pi_val = theta * (sum_t / n_pairs)
      s_muts = Samplers.poisson(0.5 * theta * total_len, rng)
      { pi: pi_val, s: s_muts, t_mrca: t_coal }
    end
  end
end
