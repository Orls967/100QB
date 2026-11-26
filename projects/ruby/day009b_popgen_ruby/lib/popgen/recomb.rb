# frozen_string_literal: true

require_relative 'samplers'

module PopGen
  # Two-locus coalescent with recombination (ARG-lite) and Integer bitset genome representation (V3).
  module Recomb
    module_function

    # Simulate two-locus ARG-lite for n samples with scaled recombination rate rho = 4*N*r (or 2*N*r).
    # Each lineage carries a 2-bit bitmask of ancestral material: bit 0 = locus 1, bit 1 = locus 2.
    # mask == 3 (0b11) means lineage is ancestral at BOTH loci (and eligible to recombine at rate rho/2!).
    # Exact event rates when k lineages remain and k2 of them carry both loci:
    # - Coalescence rate = k * (k - 1) / 2
    # - Recombination rate = k2 * (rho / 2.0)  (consistent with E[T_between] / correlation 1/(1+...)!)
    # Returns { t_mrca_1: Float, t_mrca_2: Float, len_1: Float, len_2: Float, n_recomb: Int }.
    def simulate_arg_lite(n, rho, rng)
      # Each lineage: [mask, leaves_1_count, leaves_2_count]
      lineages = Array.new(n) { [0b11, 1, 1] }
      t = 0.0
      t_mrca_1 = nil
      t_mrca_2 = nil
      len_1 = 0.0
      len_2 = 0.0
      n_recomb = 0

      loop do
        c1 = 0
        c2 = 0
        k2 = 0
        lineages.each do |mask, _, _|
          c1 += 1 unless (mask & 0b01).zero?
          c2 += 1 unless (mask & 0b10).zero?
          k2 += 1 if mask == 0b11
        end
        t_mrca_1 = t if t_mrca_1.nil? && c1 <= 1
        t_mrca_2 = t if t_mrca_2.nil? && c2 <= 1
        break if c1 <= 1 && c2 <= 1

        k = lineages.length
        coal_rate = 0.5 * k * (k - 1)
        rec_rate = 0.5 * rho * k2
        tot_rate = coal_rate + rec_rate
        dt = Samplers.exponential(tot_rate, rng)
        t += dt
        len_1 += c1 * dt if c1 > 1
        len_2 += c2 * dt if c2 > 1

        if rng.rand < (coal_rate / tot_rate)
          i = rng.rand(k)
          j = rng.rand(k - 1)
          j += 1 if j >= i
          i, j = j, i if i > j
          b_lin = lineages.delete_at(j)
          a_lin = lineages.delete_at(i)
          new_mask = a_lin[0] | b_lin[0]
          new_l1 = a_lin[1] + b_lin[1]
          new_l2 = a_lin[2] + b_lin[2]
          # Drop locus from mask once all n samples at that locus have coalesced into this node
          new_mask &= ~0b01 if new_l1 == n
          new_mask &= ~0b10 if new_l2 == n
          lineages << [new_mask, new_l1, new_l2] unless new_mask.zero? && new_l1 == n && new_l2 == n
          if new_l1 == n && t_mrca_1.nil?
            t_mrca_1 = t
          end
          if new_l2 == n && t_mrca_2.nil?
            t_mrca_2 = t
          end
        else
          # Recombine one uniformly chosen two-locus lineage (mask == 0b11)
          two_idxs = []
          lineages.each_with_index { |lin, idx| two_idxs << idx if lin[0] == 0b11 }
          pick = two_idxs[rng.rand(two_idxs.length)]
          _, l1, l2 = lineages.delete_at(pick)
          lineages << [0b01, l1, 0]
          lineages << [0b10, 0, l2]
          n_recomb += 1
        end
      end

      { t_mrca_1: t_mrca_1 || t, t_mrca_2: t_mrca_2 || t, len_1: len_1, len_2: len_2, n_recomb: n_recomb }
    end

    # V3 Optimization Ladder benchmark helper: represent L-locus haplotypes as Ruby Integer bitsets
    # with fast bitwise crossover ((g1 & mask) | (g2 & ~mask)) vs Array-of-bits crossover.
    def benchmark_bitset_recomb(n_ind, n_loci, gens, rng)
      mask_full = (1 << n_loci) - 1
      pop = Array.new(n_ind) { rng.rand(mask_full + 1) }
      gens.times do
        next_pop = Array.new(n_ind, 0)
        n_ind.times do |i|
          g1 = pop[rng.rand(n_ind)]
          g2 = pop[rng.rand(n_ind)]
          cut = rng.rand(n_loci)
          mask_lo = (1 << cut) - 1
          next_pop[i] = (g1 & mask_lo) | (g2 & (mask_full ^ mask_lo))
        end
        pop = next_pop
      end
      pop
    end

    def benchmark_array_recomb(n_ind, n_loci, gens, rng)
      pop = Array.new(n_ind) { Array.new(n_loci) { rng.rand(2) } }
      gens.times do
        next_pop = Array.new(n_ind)
        n_ind.times do |i|
          g1 = pop[rng.rand(n_ind)]
          g2 = pop[rng.rand(n_ind)]
          cut = rng.rand(n_loci)
          child = Array.new(n_loci)
          n_loci.times { |l| child[l] = (l < cut) ? g1[l] : g2[l] }
          next_pop[i] = child
        end
        pop = next_pop
      end
      pop
    end
  end
end
