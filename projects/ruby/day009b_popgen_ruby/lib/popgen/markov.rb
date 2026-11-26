# frozen_string_literal: true

require 'rational'

module PopGen
  # Exact Markov chains for Wright-Fisher (Rational for N<=12, Float for N<=40),
  # Moran birth-death chains, stationary distributions, and Island coalescent linear systems.
  module Markov
    module_function

    def choose_int(n, k)
      return 1 if k <= 0 || k >= n

      k = n - k if k > n - k
      num = 1
      den = 1
      k.times do |i|
        num *= (n - i)
        den *= (i + 1)
      end
      num / den
    end

    # Gaussian elimination solving A * x = b for exact Rational or Float matrices
    def solve_linear_rational(a_mat, b_vec)
      n = b_vec.length
      m = a_mat.each_with_index.map { |row, i| row.map(&:to_r) + [b_vec[i].to_r] }

      n.times do |col|
        piv = (col...n).find { |r| !m[r][col].zero? }
        raise 'Singular Rational matrix' unless piv

        m[col], m[piv] = m[piv], m[col] if piv != col
        inv_p = Rational(1, 1) / m[col][col]
        ((col)..n).each { |j| m[col][j] *= inv_p }

        n.times do |r|
          next if r == col || m[r][col].zero?

          fac = m[r][col]
          (col..n).each { |j| m[r][j] -= fac * m[col][j] }
        end
      end
      Array.new(n) { |i| m[i][n] }
    end

    def solve_linear_float(a_mat, b_vec)
      n = b_vec.length
      m = a_mat.each_with_index.map { |row, i| row.map(&:to_f) + [b_vec[i].to_f] }

      n.times do |col|
        piv = (col...n).max_by { |r| m[r][col].abs }
        raise 'Singular Float matrix' if m[piv][col].abs < 1e-15

        m[col], m[piv] = m[piv], m[col] if piv != col
        inv_p = 1.0 / m[col][col]
        (col..n).each { |j| m[col][j] *= inv_p }

        n.times do |r|
          next if r == col || m[r][col].abs < 1e-18

          fac = m[r][col]
          (col..n).each { |j| m[r][j] -= fac * m[col][j] }
        end
      end
      Array.new(n) { |i| m[i][n] }
    end

    # Effective probability p_prime for next generation given count i out of n_copies:
    # Supports selection s, dominance h (additive h=0.5), and bidirectional mutation (u: A->a, v: a->A)
    def wf_next_prob_rational(i, n_copies, s: 0.to_r, h: Rational(1, 2), u: 0.to_r, v: 0.to_r)
      p = Rational(i, n_copies)
      q = 1.to_r - p
      w_aa = 1.to_r + s.to_r
      w_ab = 1.to_r + h.to_r * s.to_r
      w_bb = 1.to_r
      w_bar = p * p * w_aa + 2 * p * q * w_ab + q * q * w_bb
      p_sel = (p * p * w_aa + p * q * w_ab) / w_bar
      p_sel * (1.to_r - u.to_r) + (1.to_r - p_sel) * v.to_r
    end

    def wf_next_prob_float(i, n_copies, s: 0.0, h: 0.5, u: 0.0, v: 0.0)
      p = i.to_f / n_copies
      q = 1.0 - p
      w_aa = 1.0 + s
      w_ab = 1.0 + h * s
      w_bar = p * p * w_aa + 2.0 * p * q * w_ab + q * q
      p_sel = (p * p * w_aa + p * q * w_ab) / w_bar
      p_sel * (1.0 - u) + (1.0 - p_sel) * v
    end

    def wf_transition_matrix_rational(n_copies, s: 0.to_r, h: Rational(1, 2), u: 0.to_r, v: 0.to_r)
      Array.new(n_copies + 1) do |i|
        p = wf_next_prob_rational(i, n_copies, s: s, h: h, u: u, v: v)
        q = 1.to_r - p
        Array.new(n_copies + 1) do |j|
          if p.zero?
            j.zero? ? 1.to_r : 0.to_r
          elsif q.zero?
            j == n_copies ? 1.to_r : 0.to_r
          else
            choose_int(n_copies, j) * (p**j) * (q**(n_copies - j))
          end
        end
      end
    end

    def wf_transition_matrix_float(n_copies, s: 0.0, h: 0.5, u: 0.0, v: 0.0)
      Array.new(n_copies + 1) do |i|
        p = wf_next_prob_float(i, n_copies, s: s, h: h, u: u, v: v)
        Array.new(n_copies + 1) do |j|
          Binomial.pmf(j, n_copies, p)
        end
      end
    end

    # Fixation probability u_i for transient states i=1..n-1 with u_0=0, u_n=1
    def fixation_probs(p_mat, rational: false)
      n = p_mat.length - 1
      m = n - 1
      a = Array.new(m) { Array.new(m, rational ? 0.to_r : 0.0) }
      b = Array.new(m, rational ? 0.to_r : 0.0)
      m.times do |idx|
        i = idx + 1
        m.times do |jdx|
          j = jdx + 1
          a[idx][jdx] = (i == j ? (rational ? 1.to_r : 1.0) : (rational ? 0.to_r : 0.0)) - p_mat[i][j]
        end
        b[idx] = p_mat[i][n]
      end
      sol = rational ? solve_linear_rational(a, b) : solve_linear_float(a, b)
      [rational ? 0.to_r : 0.0] + sol + [rational ? 1.to_r : 1.0]
    end

    # Mean unconditional absorption time t_i (to 0 or N) for i=1..n-1
    def mean_absorption_times(p_mat, rational: false)
      n = p_mat.length - 1
      m = n - 1
      a = Array.new(m) { Array.new(m, rational ? 0.to_r : 0.0) }
      b = Array.new(m, rational ? 1.to_r : 1.0)
      m.times do |idx|
        i = idx + 1
        m.times do |jdx|
          j = jdx + 1
          a[idx][jdx] = (i == j ? (rational ? 1.to_r : 1.0) : (rational ? 0.to_r : 0.0)) - p_mat[i][j]
        end
      end
      sol = rational ? solve_linear_rational(a, b) : solve_linear_float(a, b)
      [rational ? 0.to_r : 0.0] + sol + [rational ? 0.to_r : 0.0]
    end

    # Mean conditional fixation time (conditioned on hitting state N):
    # Solve (I - Q) * (u * t_fix) = u  =>  w_i = u_i * t_fix(i), where (I - Q) w = u
    def mean_conditional_fixation_times(p_mat, rational: false)
      u = fixation_probs(p_mat, rational: rational)
      n = p_mat.length - 1
      m = n - 1
      a = Array.new(m) { Array.new(m, rational ? 0.to_r : 0.0) }
      b = Array.new(m) { |idx| u[idx + 1] }
      m.times do |idx|
        i = idx + 1
        m.times do |jdx|
          j = jdx + 1
          a[idx][jdx] = (i == j ? (rational ? 1.to_r : 1.0) : (rational ? 0.to_r : 0.0)) - p_mat[i][j]
        end
      end
      w = rational ? solve_linear_rational(a, b) : solve_linear_float(a, b)
      res = Array.new(n + 1, rational ? 0.to_r : 0.0)
      m.times { |idx| res[idx + 1] = w[idx] / u[idx + 1] }
      res
    end

    # Stationary distribution of irreducible Markov chain (u, v > 0): (P^T - I) pi = 0, sum(pi) = 1
    def stationary_distribution(p_mat, rational: false)
      d = p_mat.length
      a = Array.new(d) { Array.new(d, rational ? 0.to_r : 0.0) }
      b = Array.new(d, rational ? 0.to_r : 0.0)
      (0...(d - 1)).each do |i|
        d.times do |j|
          a[i][j] = p_mat[j][i] - (i == j ? (rational ? 1.to_r : 1.0) : (rational ? 0.to_r : 0.0))
        end
      end
      d.times { |j| a[d - 1][j] = rational ? 1.to_r : 1.0 }
      b[d - 1] = rational ? 1.to_r : 1.0
      rational ? solve_linear_rational(a, b) : solve_linear_float(a, b)
    end

    # Exact Moran transition matrix with relative fitness r for allele A (birth proportional to fitness, death uniform)
    def moran_transition_matrix(n, r)
      p_mat = Array.new(n + 1) { Array.new(n + 1, 0.0) }
      p_mat[0][0] = 1.0
      p_mat[n][n] = 1.0
      (1...n).each do |i|
        w_tot = r * i + (n - i)
        p_birth_a = (r * i) / w_tot
        p_birth_b = (n - i) / w_tot
        p_death_a = i.to_f / n
        p_death_b = (n - i).to_f / n
        p_up = p_birth_a * p_death_b
        p_dn = p_birth_b * p_death_a
        p_mat[i][i + 1] = p_up
        p_mat[i][i - 1] = p_dn
        p_mat[i][i] = 1.0 - p_up - p_dn
      end
      p_mat
    end

    # Closed-form Moran fixation probability from i=1..n-1 with relative fitness r
    def moran_fixation_closed(i, n, r)
      return i.to_f / n if (r - 1.0).abs < 1e-12

      (1.0 - (1.0 / r)**i) / (1.0 - (1.0 / r)**n)
    end

    # Exact 2-sample coalescent times (T_within, T_between) in d-deme finite island model
    # In coalescent units (1 unit = 2N generations for diploid / N for haploid),
    # migration rate per lineage into each other deme is m_scaled = M / (2*(d-1)), where M = 2*N_deme*m_tot.
    # For discrete generation model with per-lineage migration prob m per generation into a random other deme:
    # Continuous-time rates in units of N_deme generations:
    #   In same deme: coal rate = 1, migration out of either lineage = 2 * N_deme * m = M.
    #     Wait: if M_lineage = N_deme * m per lineage, total rate leaving (within) state is 1 + 2*M_lineage.
    #   From (within): E[T_w] = 1 / (1 + 2*M_l) + (2*M_l / (1 + 2*M_l)) * E[T_b]
    #   From (between): lineages migrate at rate 2*M_l, hitting same deme with prob 1/(d-1):
    #     E[T_b] = (d - 1) / (2 * M_l) + E[T_w]
    def island_exact_times(d, m_lineage)
      # Solve 2x2 linear system:
      # (1 + 2*M_l) * Tw - 2*M_l * Tb = 1
      # -2*M_l/(d-1) * Tw + 2*M_l/(d-1) * Tb = 1
      rate_b = 2.0 * m_lineage / (d - 1)
      a = [
        [1.0 + 2.0 * m_lineage, -2.0 * m_lineage],
        [-rate_b, rate_b]
      ]
      b = [1.0, 1.0]
      tw, tb = solve_linear_float(a, b)
      # Slatkin / Hudson F_ST = 1 - Tw / ((1/d)*Tw + ((d-1)/d)*Tb)  or Hudson F_ST = (Tb - Tw) / Tb
      t_tot = (tw + (d - 1) * tb) / d.to_f
      fst_slatkin = (t_tot - tw) / t_tot
      fst_hudson = (tb - tw) / tb
      { tw: tw, tb: tb, t_total: t_tot, fst_slatkin: fst_slatkin, fst_hudson: fst_hudson }
    end
  end
end
