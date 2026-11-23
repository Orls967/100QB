# frozen_string_literal: true

module CopulaTailDependence
  module Importance
    module_function

    # Anchor 1D kuadratur Gauss-Legendre untuk faktor tunggal Gaussian:
    # Pr(Z1 <= a, ..., Zd <= a) = int Phi((a - sqrt(rho)*f) / sqrt(1-rho))^d phi(f) df.
    def exact_single_factor_prob(d, a, rho)
      return Special.bivariate_norm_cdf(a, a, rho) if d == 2

      nodes, weights = Special.gauss_legendre(64)
      half = 8.0
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)
      sum = 0.0
      nodes.each_with_index do |xi, idx|
        f = half * xi
        u = (a - (sr * f)) / s1r
        p = Special.norm_cdf(u)
        sum += weights[idx] * (p**d) * Special.norm_pdf(f)
      end
      half * sum
    end

    # Pergeseran rerata optimal analitik (saddlepoint / mode log-integran) pada faktor tunggal:
    # mu* ~ d * sqrt(rho) * a / (1 - rho + d * rho), dipoles 5 langkah Newton.
    def optimal_factor_shift(d, a, rho)
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)
      f = (d * sr * a) / (1.0 - rho + (d * rho))
      5.times do
        u = (a - (sr * f)) / s1r
        cdf = Special.norm_cdf(u)
        break if cdf < 1.0e-300

        ratio = Special.norm_pdf(u) / cdf
        g1 = (-(d * sr / s1r) * ratio) - f
        g2 = ((d * rho / (1.0 - rho)) * (-ratio * (u + ratio))) - 1.0
        f -= g1 / g2
      end
      f
    end

    def simulate_is_chunk(len, d, a, sr, s1r, mu_star, rng)
      ns1 = 0.0
      ns2 = 0.0
      is1 = 0.0
      is2 = 0.0
      i = 0
      while i < len
        f_n = rng.next_normal
        all_le = true
        k = 0
        while k < d
          if ((sr * f_n) + (s1r * rng.next_normal)) > a
            all_le = false
            break
          end
          k += 1
        end
        ind = all_le ? 1.0 : 0.0
        ns1 += ind
        ns2 += ind * ind

        f_is = mu_star + rng.next_normal
        lr = Math.exp((-mu_star * f_is) + (0.5 * mu_star * mu_star))
        est = lr * (Special.norm_cdf((a - (sr * f_is)) / s1r)**d)
        is1 += est
        is2 += est * est
        i += 1
      end
      [ns1, ns2, is1, is2]
    end
    private_class_method :simulate_is_chunk

    def compare_naive_vs_is(d:, a:, rho:, n_paths: 100_000, seed: MAIN_SEED, workers: nil)
      exact = exact_single_factor_prob(d, a, rho)
      mu_star = optimal_factor_shift(d, a, rho)
      sr = Math.sqrt(rho)
      s1r = Math.sqrt(1.0 - rho)

      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 10_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 400 + (d * 10) + (rho * 10).round, cidx))
        simulate_is_chunk(len, d, a, sr, s1r, mu_star, rng)
      end

      ns1 = chunks.sum(0.0) { |c| c[0] }
      ns2 = chunks.sum(0.0) { |c| c[1] }
      is1 = chunks.sum(0.0) { |c| c[2] }
      is2 = chunks.sum(0.0) { |c| c[3] }

      n_mean = ns1 / n_paths
      n_var = [(ns2 / n_paths) - (n_mean * n_mean), 0.0].max
      i_mean = is1 / n_paths
      i_var = [(is2 / n_paths) - (i_mean * i_mean), 0.0].max
      i_se = Math.sqrt(i_var / n_paths)

      {
        d: d, a: a, rho: rho, exact: exact, mu_star: mu_star,
        naive_mean: n_mean, naive_se: Math.sqrt(n_var / n_paths), naive_var: n_var,
        is_mean: i_mean, is_se: i_se, is_var: i_var,
        is_z: (i_mean - exact) / [i_se, 1.0e-18].max,
        variance_reduction: i_var.positive? ? n_var / i_var : Float::INFINITY
      }
    end
  end
end
