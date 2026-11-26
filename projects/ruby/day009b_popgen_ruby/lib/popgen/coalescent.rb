# frozen_string_literal: true

require_relative 'samplers'

module PopGen
  # Kingman coalescent tree generator, site-frequency spectrum, Watterson/Tajima statistics,
  # and Ewens sampling formula (exact partitions + allele count distribution).
  module Coalescent
    module_function

    def harmonic(n_minus_1)
      return 0.0 if n_minus_1 <= 0

      (1..n_minus_1).sum(0.0) { |i| 1.0 / i }
    end

    def harmonic_sq(n_minus_1)
      return 0.0 if n_minus_1 <= 0

      (1..n_minus_1).sum(0.0) { |i| 1.0 / (i * i) }
    end

    # Tajima (1989) theoretical constants for sample size n
    def tajima_constants(n)
      a1 = harmonic(n - 1)
      a2 = harmonic_sq(n - 1)
      b1 = (n + 1).to_f / (3.0 * (n - 1))
      b2 = 2.0 * (n * n + n + 3).to_f / (9.0 * n * (n - 1))
      c1 = b1 - 1.0 / a1
      c2 = b2 - (n + 2).to_f / (a1 * n) + a2 / (a1 * a1)
      e1 = c1 / a1
      e2 = c2 / (a1 * a1 + a2)
      { a1: a1, a2: a2, b1: b1, b2: b2, c1: c1, c2: c2, e1: e1, e2: e2 }
    end

    # Simulate a Kingman coalescent tree for n samples (with optional exponential growth rate beta >= 0)
    # In standard Kingman (beta = 0): while k lineages remain, waiting time ~ Exp(k*(k-1)/2).
    # Returns tree metrics and unfolded SFS array of length n - 1 (derived counts 1..n-1).
    def simulate_tree(n, theta, rng, growth_beta: 0.0)
      # Active nodes: { size: leaf_count }
      nodes = Array.new(n) { { size: 1 } }
      t_now = 0.0
      total_len = 0.0
      # Branch lengths weighted by descendant count for SFS mutation placement
      branches = []

      k = n
      while k > 1
        rate = 0.5 * k * (k - 1)
        dt = if growth_beta <= 0.0
               Samplers.exponential(rate, rng)
             else
               # Exact time-inhomogeneous coalescent with N(t) = N(0)*exp(-beta*t) backwards:
               # integral_{t_now}^{t_now+dt} rate * exp(beta * u) du = E ~ Exp(1)
               e_draw = Samplers.exponential(1.0, rng)
               arg = 1.0 + (growth_beta * e_draw) / (rate * Math.exp(growth_beta * t_now))
               Math.log(arg) / growth_beta
             end
        t_now += dt
        total_len += k * dt
        nodes.each { |nd| branches << [dt, nd[:size]] }

        # Merge two uniformly chosen lineages
        i = rng.rand(k)
        j = rng.rand(k - 1)
        j += 1 if j >= i
        i, j = j, i if i > j
        n_j = nodes.delete_at(j)
        n_i = nodes.delete_at(i)
        nodes << { size: n_i[:size] + n_j[:size] }
        k -= 1
      end

      # Place mutations along branches as Poisson(0.5 * theta * branch_length)
      sfs = Array.new(n - 1, 0)
      s_tot = 0
      branches.each do |blen, sz|
        next if sz <= 0 || sz >= n

        m = Samplers.poisson(0.5 * theta * blen, rng)
        if m.positive?
          sfs[sz - 1] += m
          s_tot += m
        end
      end

      # Pairwise nucleotide diversity pi from SFS: pi = sum_{i=1}^{n-1} 2*i*(n-i)/(n*(n-1)) * xi_i
      denom_pairs = 0.5 * n * (n - 1)
      pi_val = 0.0
      sfs.each_with_index do |cnt, idx|
        i = idx + 1
        pi_val += cnt * i * (n - i).to_f / denom_pairs if cnt.positive?
      end

      tc = tajima_constants(n)
      theta_w = s_tot.to_f / tc[:a1]
      d_taj = if s_tot.zero?
                0.0
              else
                denom = Math.sqrt(tc[:e1] * s_tot + tc[:e2] * s_tot * (s_tot - 1))
                denom > 1e-15 ? ((pi_val - theta_w) / denom) : 0.0
              end

      {
        t_mrca: t_now,
        total_length: total_len,
        s: s_tot,
        pi: pi_val,
        theta_w: theta_w,
        tajima_d: d_taj,
        sfs: sfs
      }
    end

    # Hoppe urn simulator for Ewens sampling formula (exact allele partition for sample size n)
    # Returns partition vector a = [a_1, a_2, ..., a_n] where a_j = # alleles with j copies, and K_n.
    def ewens_sample_urn(n, theta, rng)
      groups = [1]
      (2..n).each do |i|
        u = rng.rand * (theta + i - 1)
        if u < theta
          groups << 1
        else
          # Choose existing copy uniformly at random
          target = rng.rand(i - 1)
          acc = 0
          groups.each_with_index do |sz, g_idx|
            acc += sz
            if target < acc
              groups[g_idx] += 1
              break
            end
          end
        end
      end
      part = Array.new(n, 0)
      groups.each { |sz| part[sz - 1] += 1 }
      { partition: part, k_n: groups.length }
    end

    # Enumerate all integer partitions of n as arrays [a_1, ..., a_n] (sum j*a_j = n)
    def all_partitions(n)
      res = []
      rec = lambda do |rem, max_j, cur|
        if rem.zero?
          res << cur.dup
          return
        end
        j = [rem, max_j].min
        while j >= 1
          cur[j - 1] += 1
          rec.call(rem - j, j, cur)
          cur[j - 1] -= 1
          j -= 1
        end
      end
      rec.call(n, n, Array.new(n, 0))
      res
    end

    def fact(k)
      return 1 if k <= 1

      (2..k).reduce(1, :*)
    end

    # Exact Ewens Sampling Formula probability P(a_1, ..., a_n; theta)
    def ewens_partition_prob(part, n, theta)
      rising = (0...n).reduce(1.0) { |acc, i| acc * (theta + i) }
      k_tot = part.sum
      den = 1.0
      part.each_with_index do |aj, idx|
        j = idx + 1
        den *= (j**aj) * fact(aj) if aj.positive?
      end
      (fact(n).to_f / rising) * (theta**k_tot) / den
    end

    def ewens_expected_k(n, theta)
      (0...n).sum(0.0) { |i| theta / (theta + i) }
    end
  end
end
