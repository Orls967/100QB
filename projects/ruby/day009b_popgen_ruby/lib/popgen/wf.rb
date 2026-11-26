# frozen_string_literal: true

require_relative 'binomial'
require_relative 'samplers'

module PopGen
  # Wright-Fisher simulators (haploid/diploid, selection, mutation, infinite alleles, V0/V1).
  module WF
    module_function

    # Step 1 generation for biallelic WF with n_copies gene copies (N for haploid, 2N for diploid)
    def step_biallelic(count, n_copies, rng, s: 0.0, h: 0.5, u: 0.0, v: 0.0)
      p = count.to_f / n_copies
      if s != 0.0
        q = 1.0 - p
        w_aa = 1.0 + s
        w_ab = 1.0 + h * s
        w_bar = p * p * w_aa + 2.0 * p * q * w_ab + q * q
        p = (p * p * w_aa + p * q * w_ab) / w_bar
      end
      p = p * (1.0 - u) + (1.0 - p) * v if u != 0.0 || v != 0.0
      Binomial.sample(n_copies, p, rng)
    end

    # V0: Object-per-individual Hash simulation (used in optimization ladder & micro-checks)
    def simulate_v0_hash(n_copies, p0, gens, rng, s: 0.0, u: 0.0, v: 0.0)
      n_a = (n_copies * p0).round
      pop = Array.new(n_copies) { |i| { id: i, allele: (i < n_a ? 1 : 0) } }
      gens.times do
        cnt = pop.count { |ind| ind[:allele] == 1 }
        p = cnt.to_f / n_copies
        p = (p * (1.0 + s)) / (1.0 + s * p) if s != 0.0
        p = p * (1.0 - u) + (1.0 - p) * v if u != 0.0 || v != 0.0
        pop = Array.new(n_copies) { |i| { id: i, allele: (rng.rand < p ? 1 : 0) } }
      end
      pop.count { |ind| ind[:allele] == 1 }.to_f / n_copies
    end

    # V1: Flat allele-frequency Array simulation returning [p_0, p_1, ..., p_gens]
    def simulate_v1_trajectory(n_copies, p0, gens, rng, s: 0.0, h: 0.5, u: 0.0, v: 0.0)
      traj = Array.new(gens + 1, 0.0)
      cnt = (n_copies * p0).round
      traj[0] = cnt.to_f / n_copies
      gens.times do |g|
        cnt = step_biallelic(cnt, n_copies, rng, s: s, h: h, u: u, v: v)
        traj[g + 1] = cnt.to_f / n_copies
      end
      traj
    end

    # Simulate neutral/selected biallelic WF until absorption (or max_gens)
    # Returns { fixed: bool, lost: bool, censored: bool, gens: int }
    def simulate_to_absorption(n_copies, init_count, rng, s: 0.0, h: 0.5, max_gens: 100_000)
      cnt = init_count
      t = 0
      while cnt.positive? && cnt < n_copies && t < max_gens
        cnt = step_biallelic(cnt, n_copies, rng, s: s, h: h)
        t += 1
      end
      { fixed: cnt == n_copies, lost: cnt.zero?, censored: cnt.positive? && cnt < n_copies, gens: t }
    end

    # Finite-N Wright-Fisher infinite-alleles simulator:
    # Tracks allele counts Hash { allele_id => count } over gens generations.
    def simulate_infinite_alleles(n_copies, mu, gens, rng, burn_in: 0)
      counts = { 0 => n_copies }
      next_id = 1
      h_series = []
      total = burn_in + gens

      total.times do |g|
        # Multinomial offspring sampling across existing alleles
        ids = counts.keys
        probs = ids.map { |id| counts[id].to_f / n_copies }
        drawn = Samplers.multinomial(n_copies, probs, rng)
        new_counts = {}
        ids.each_with_index do |id, idx|
          c = drawn[idx]
          new_counts[id] = c if c.positive?
        end

        # Apply Poisson(n_copies * mu) mutations to randomly chosen copies
        n_mut = Samplers.poisson(n_copies * mu, rng)
        if n_mut.positive?
          arr = []
          new_counts.each { |id, c| c.times { arr << id } }
          n_mut.times do
            slot = rng.rand(n_copies)
            arr[slot] = next_id
            next_id += 1
          end
          new_counts = Hash.new(0)
          arr.each { |id| new_counts[id] += 1 }
        end
        counts = new_counts

        if g >= burn_in
          sum_sq = counts.values.sum(0.0) { |c| (c.to_f / n_copies)**2 }
          h_series << (1.0 - sum_sq)
        end
      end
      { counts: counts, h_series: h_series, k_alleles: counts.length }
    end

    # Direct finite-N WF genealogical forward simulation for n_sample lineages
    # Used in S41 to compare direct WF forward genealogy against Kingman coalescent!
    def simulate_wf_genealogy_sample(n_copies, n_sample, theta, gens_max, rng)
      # Trace backwards n_sample lineages through discrete WF generations (exact parent sampling)
      branches = Array.new(n_sample, 0.0)
      # Each lineage accumulates generations; when two coalesce we record subtree
      subtrees = Array.new(n_sample) { |i| { leaves: [i], len: 0.0 } }
      t = 0
      t_mrca = nil

      while subtrees.length > 1 && t < gens_max
        t += 1
        subtrees.each { |st| st[:len] += 1.0 }
        # Assign each active lineage a uniform parent in 0...n_copies
        groups = Hash.new { |h, k| h[k] = [] }
        subtrees.each do |st|
          parent = rng.rand(n_copies)
          groups[parent] << st
        end
        next_subtrees = []
        groups.each_value do |g_list|
          if g_list.length == 1
            next_subtrees << g_list[0]
          else
            # Record lengths of all coalescing children before merging into ancestor
            g_list.each { |child| branches << child[:len] }
            merged_leaves = g_list.flat_map { |c| c[:leaves] }
            next_subtrees << { leaves: merged_leaves, len: 0.0 }
          end
        end
        subtrees = next_subtrees
        t_mrca = t if subtrees.length == 1
      end
      # Collect lengths of branches from coalescence events (exclude initial dummy slots)
      coal_lens = branches[n_sample..] || []
      total_len_coal = coal_lens.sum(0.0) / n_copies.to_f
      t_mrca_coal = (t_mrca || t).to_f / n_copies.to_f
      # Sprinkle mutations at rate theta / 2 per coalescent unit
      s_muts = Samplers.poisson(0.5 * theta * total_len_coal, rng)
      { t_mrca: t_mrca_coal, total_length: total_len_coal, s: s_muts }
    end
  end
end
