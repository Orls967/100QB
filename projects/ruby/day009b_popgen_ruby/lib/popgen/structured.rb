# frozen_string_literal: true

require_relative 'samplers'

module PopGen
  # Structured coalescent in d-deme finite island and 1D stepping-stone models.
  module Structured
    module_function

    # Simulate 2-sample or n-sample structured coalescent in d-deme island model
    # Migration rate per lineage in coalescent time (units of N_deme generations) is m_scaled.
    # For 2 samples (within same deme or between two demes):
    #   In same deme: coal rate = 1, migration out of either lineage = 2 * m_scaled.
    #   In different demes: migration rate = 2 * m_scaled; hits same deme with prob 1/(d-1).
    def simulate_island_two_sample(d, m_scaled, same_deme, rng)
      t = 0.0
      in_same = same_deme
      loop do
        if in_same
          tot_rate = 1.0 + 2.0 * m_scaled
          t += Samplers.exponential(tot_rate, rng)
          if rng.rand < (1.0 / tot_rate)
            return t
          else
            in_same = false
          end
        else
          # Only migration that lands the migrating lineage into the other's deme changes state,
          # OR simulate every migration event at rate 2 * m_scaled and switch with prob 1/(d-1)
          tot_rate = 2.0 * m_scaled
          t += Samplers.exponential(tot_rate, rng)
          in_same = true if rng.rand < (1.0 / (d - 1))
        end
      end
    end

    # Multi-sample structured coalescent in d-deme island or 1D stepping-stone (ring) model
    # sample_counts: array of length d giving initial samples per deme.
    # Returns average within-deme T_w, between-deme T_b, Hudson F_ST, and Nei G_ST.
    def simulate_structured_pairwise(d, sample_counts, m_scaled, topology, theta_per_deme, rng)
      # Represent each active lineage as [deme_idx, leaf_set]
      lineages = []
      leaf_deme = []
      leaf_id = 0
      d.times do |dm|
        sample_counts[dm].times do
          lineages << [dm, [leaf_id]]
          leaf_deme << dm
          leaf_id += 1
        end
      end
      n_tot = leaf_id
      t_pair = Array.new(n_tot) { Array.new(n_tot, 0.0) }
      t = 0.0

      while lineages.length > 1
        counts = Array.new(d, 0)
        lineages.each { |dm, _| counts[dm] += 1 }
        coal_rate = counts.sum(0.0) { |c| 0.5 * c * (c - 1) }
        mig_rate = lineages.length * m_scaled
        tot_rate = coal_rate + mig_rate
        t += Samplers.exponential(tot_rate, rng)

        if rng.rand < (coal_rate / tot_rate)
          # Pick deme weighted by c*(c-1)/2
          u = rng.rand * coal_rate
          acc = 0.0
          chosen_d = 0
          d.times do |dm|
            acc += 0.5 * counts[dm] * (counts[dm] - 1)
            if u <= acc
              chosen_d = dm
              break
            end
          end
          idxs = []
          lineages.each_with_index { |lin, idx| idxs << idx if lin[0] == chosen_d }
          i_pos = rng.rand(idxs.length)
          j_pos = rng.rand(idxs.length - 1)
          j_pos += 1 if j_pos >= i_pos
          idx_a = idxs[i_pos]
          idx_b = idxs[j_pos]
          idx_a, idx_b = idx_b, idx_a if idx_a > idx_b
          lin_b = lineages.delete_at(idx_b)
          lin_a = lineages.delete_at(idx_a)

          # Record pairwise coalescence times
          lin_a[1].each do |la|
            lin_b[1].each do |lb|
              t_pair[la][lb] = t
              t_pair[lb][la] = t
            end
          end
          lineages << [chosen_d, lin_a[1] + lin_b[1]]
        else
          # Migrate a uniformly chosen lineage
          l_idx = rng.rand(lineages.length)
          cur_d = lineages[l_idx][0]
          new_d = if topology == :stepping_stone
                    step = rng.rand < 0.5 ? -1 : 1
                    (cur_d + step) % d
                  else
                    nd = rng.rand(d - 1)
                    nd += 1 if nd >= cur_d
                    nd
                  end
          lineages[l_idx][0] = new_d
        end
      end

      # Compute mean within-deme and between-deme pairwise times (and pi = theta * T)
      tw_sum = 0.0
      tw_cnt = 0
      tb_sum = 0.0
      tb_cnt = 0
      within_by_deme = Array.new(d) { [0.0, 0] }
      between_by_pair = Hash.new { |h, k| h[k] = [0.0, 0] }

      (0...n_tot).each do |a|
        ((a + 1)...n_tot).each do |b|
          da = leaf_deme[a]
          db = leaf_deme[b]
          val = t_pair[a][b]
          if da == db
            tw_sum += val
            tw_cnt += 1
            within_by_deme[da][0] += val
            within_by_deme[da][1] += 1
          else
            tb_sum += val
            tb_cnt += 1
            key = da < db ? [da, db] : [db, da]
            between_by_pair[key][0] += val
            between_by_pair[key][1] += 1
          end
        end
      end

      tw = tw_cnt.positive? ? (tw_sum / tw_cnt) : 0.0
      tb = tb_cnt.positive? ? (tb_sum / tb_cnt) : 0.0
      # Total pairwise time T_total (uniform pair across all n_tot samples)
      t_tot_all = (tw_sum + tb_sum) / (0.5 * n_tot * (n_tot - 1))
      fst_hudson = tb > 1e-12 ? ((tb - tw) / tb) : 0.0
      fst_nei = t_tot_all > 1e-12 ? ((t_tot_all - tw) / t_tot_all) : 0.0

      {
        tw: tw,
        tb: tb,
        pi_w: theta_per_deme * tw,
        pi_b: theta_per_deme * tb,
        fst_hudson: fst_hudson,
        fst_nei: fst_nei,
        between_by_pair: between_by_pair
      }
    end
  end
end
