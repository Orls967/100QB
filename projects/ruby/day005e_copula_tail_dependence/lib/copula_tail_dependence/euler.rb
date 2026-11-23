# frozen_string_literal: true

module CopulaTailDependence
  module Euler
    DEFAULT_WEIGHTS = [0.25, 0.20, 0.20, 0.15, 0.20].freeze
    DEFAULT_SIGMA = [
      [0.0400, 0.0180, 0.0120, 0.0090, 0.0140],
      [0.0180, 0.0529, 0.0160, 0.0110, 0.0150],
      [0.0120, 0.0160, 0.0361, 0.0100, 0.0130],
      [0.0090, 0.0110, 0.0100, 0.0484, 0.0120],
      [0.0140, 0.0150, 0.0130, 0.0120, 0.0441]
    ].freeze

    module_function

    def gaussian_exact(weights = DEFAULT_WEIGHTS, sigma = DEFAULT_SIGMA, alpha: 0.99)
      sw = Linalg.mat_vec(sigma, weights)
      sigma_p = Math.sqrt(Linalg.dot(weights, sw))
      z_a = Special.norm_inv(alpha)
      es_mult = Special.norm_pdf(z_a) / (1.0 - alpha)
      var_tot = sigma_p * z_a
      es_tot = sigma_p * es_mult
      c_var = weights.each_with_index.map { |wi, i| wi * sw[i] * z_a / sigma_p }
      c_es = weights.each_with_index.map { |wi, i| wi * sw[i] * es_mult / sigma_p }
      {
        sigma_p: sigma_p,
        var: var_tot,
        es: es_tot,
        var_contribs: c_var,
        es_contribs: c_es
      }
    end

    def gaussian_es_finite_diff(weights = DEFAULT_WEIGHTS, sigma = DEFAULT_SIGMA, alpha: 0.99, h: 1.0e-5)
      z_a = Special.norm_inv(alpha)
      es_mult = Special.norm_pdf(z_a) / (1.0 - alpha)
      es_fn = lambda do |w_vec|
        sw = Linalg.mat_vec(sigma, w_vec)
        Math.sqrt(Linalg.dot(w_vec, sw)) * es_mult
      end
      weights.each_with_index.map do |wi, i|
        wp = weights.dup
        wm = weights.dup
        wp[i] += h
        wm[i] -= h
        wi * (es_fn.call(wp) - es_fn.call(wm)) / (2.0 * h)
      end
    end

    # Simulasi MC Euler ES: c_i = E[L_i | L >= VaR_alpha] untuk kopula :gaussian, :t, atau :clayton.
    def mc_euler_es(copula:, weights: DEFAULT_WEIGHTS, sigma: DEFAULT_SIGMA, alpha: 0.99,
                    n_paths: 200_000, seed: MAIN_SEED, workers: nil, df: 4.0, theta: 2.0)
      d = weights.length
      chol = Linalg.cholesky(sigma)
      sd_diag = (0...d).map { |i| Math.sqrt(sigma[i][i]) }
      corr_chol = build_corr_chol(sigma, sd_diag)

      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 10_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 200 + copula_stream(copula), cidx))
        draw_euler_chunk(copula, len, d, weights, chol, corr_chol, sd_diag, df, theta, rng)
      end

      losses = Array.new(n_paths, 0.0)
      contribs_paths = Array.new(d) { Array.new(n_paths, 0.0) }
      offset = 0
      chunks.each do |chunk_losses, chunk_li|
        m = chunk_losses.length
        i = 0
        while i < m
          losses[offset + i] = chunk_losses[i]
          j = 0
          while j < d
            contribs_paths[j][offset + i] = chunk_li[j][i]
            j += 1
          end
          i += 1
        end
        offset += m
      end

      compute_tail_euler_stats(losses, contribs_paths, alpha)
    end

    def copula_stream(copula)
      { gaussian: 1, t: 2, clayton: 3 }.fetch(copula)
    end
    private_class_method :copula_stream

    def build_corr_chol(sigma, sd_diag)
      d = sigma.length
      corr = Array.new(d) do |i|
        Array.new(d) { |j| sigma[i][j] / (sd_diag[i] * sd_diag[j]) }
      end
      Linalg.cholesky(corr)
    end
    private_class_method :build_corr_chol

    def draw_euler_chunk(copula, len, d, weights, chol, corr_chol, sd_diag, df, theta, rng)
      chunk_losses = Array.new(len, 0.0)
      chunk_li = Array.new(d) { Array.new(len, 0.0) }
      z = Array.new(d, 0.0)
      i = 0
      while i < len
        k = 0
        while k < d
          z[k] = rng.next_normal
          k += 1
        end
        x_vec = sample_asset_returns(copula, d, z, chol, corr_chol, sd_diag, df, theta, rng)
        tot = 0.0
        r = 0
        while r < d
          li = -weights[r] * x_vec[r]
          chunk_li[r][i] = li
          tot += li
          r += 1
        end
        chunk_losses[i] = tot
        i += 1
      end
      [chunk_losses, chunk_li]
    end
    private_class_method :draw_euler_chunk

    def sample_asset_returns(copula, d, z, chol, corr_chol, sd_diag, df, theta, rng)
      case copula
      when :gaussian
        Linalg.mat_vec(chol, z)
      when :t
        g = Linalg.mat_vec(corr_chol, z)
        scale = Math.sqrt(df / rng.next_chi2(df))
        Array.new(d) { |r| sd_diag[r] * Special.norm_inv(Special.t_cdf(g[r] * scale, df)) }
      when :clayton
        v = rng.next_gamma(1.0 / theta, 1.0)
        Array.new(d) do |r|
          u = (1.0 + (rng.next_exp(1.0) / v))**(-1.0 / theta)
          sd_diag[r] * Special.norm_inv(u.clamp(1.0e-12, 1.0 - 1.0e-12))
        end
      else
        raise ArgumentError, "Kopula tidak dikenal: #{copula}"
      end
    end
    private_class_method :sample_asset_returns

    def compute_tail_euler_stats(losses, contribs_paths, alpha)
      n = losses.length
      d = contribs_paths.length
      var_est = Stats.quantile(losses, alpha)
      tail_idx = []
      i = 0
      while i < n
        tail_idx << i if losses[i] >= var_est
        i += 1
      end
      m = tail_idx.length
      es_est = tail_idx.sum(0.0) { |idx| losses[idx] } / m
      c_est = Array.new(d, 0.0)
      c_se = Array.new(d, 0.0)
      d.times do |j|
        vals = tail_idx.map { |idx| contribs_paths[j][idx] }
        c_est[j] = Stats.mean(vals)
        c_se[j] = Stats.se(vals)
      end
      { var: var_est, es: es_est, contribs: c_est, contrib_ses: c_se, tail_count: m }
    end
    private_class_method :compute_tail_euler_stats
  end
end
