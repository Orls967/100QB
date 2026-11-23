# frozen_string_literal: true

module FactorPortfolioLab
  # Synthetic linear multi-factor return generator with Gaussian or Student-t innovations.
  module FactorWorld
    module_function

    def build_world(n_assets: 30, n_factors: 3, seed: 20_261_123, alpha_scale: 0.0, factor_corr: 0.25)
      rng = Rng::Sampler.new(seed)
      betas = build_betas(n_assets, n_factors, rng)
      sigma_f = build_factor_cov(n_factors, factor_corr)
      mu_f = Array.new(n_factors) { |k| 0.005 - (0.001 * k) }
      alpha = build_alpha(n_assets, alpha_scale, rng)
      idio_var = Array.new(n_assets) { (0.015 + (0.025 * rng.next_float))**2 }

      sigma = compute_total_cov(betas, sigma_f, idio_var)
      mu = compute_expected_returns(alpha, betas, mu_f)

      {
        n_assets: n_assets,
        n_factors: n_factors,
        betas: betas,
        sigma_f: sigma_f,
        chol_f: LinAlg.cholesky(sigma_f),
        mu_f: mu_f,
        alpha: alpha,
        idio_var: idio_var,
        idio_sd: idio_var.map { |v| Math.sqrt(v) },
        sigma: sigma,
        mu: mu
      }
    end

    def simulate(world, t_len:, seed:, df: nil)
      rng = Rng::Sampler.new(seed)
      n_assets = world[:n_assets]
      n_factors = world[:n_factors]
      betas = world[:betas]
      chol_f = world[:chol_f]
      mu_f = world[:mu_f]
      alpha = world[:alpha]
      idio_sd = world[:idio_sd]

      factors = Array.new(t_len)
      returns = Array.new(t_len)
      eps_mat = Array.new(t_len)

      t_len.times do |t|
        z_f = Array.new(n_factors) { draw_innov(rng, df) }
        f_t = LinAlg.mat_vec(chol_f, z_f)
        n_factors.times { |k| f_t[k] += mu_f[k] }

        e_t = Array.new(n_assets) { |i| idio_sd[i] * draw_innov(rng, df) }
        r_t = Array.new(n_assets) do |i|
          val = alpha[i] + e_t[i]
          bi = betas[i]
          n_factors.times { |k| val += bi[k] * f_t[k] }
          val
        end

        factors[t] = f_t
        eps_mat[t] = e_t
        returns[t] = r_t
      end

      { factors: factors, returns: returns, eps: eps_mat }
    end

    def simulate_two_regime(world, t_len:, seed:, split_ratio: 0.5, bear_vol_mult: 2.2, bear_mu_shift: -0.008)
      rng = Rng::Sampler.new(seed)
      split_t = (t_len * split_ratio).to_i
      n_assets = world[:n_assets]
      n_factors = world[:n_factors]
      betas = world[:betas]
      chol_f = world[:chol_f]
      mu_f = world[:mu_f]
      alpha = world[:alpha]
      idio_sd = world[:idio_sd]

      factors = Array.new(t_len)
      returns = Array.new(t_len)

      t_len.times do |t|
        in_bear = t >= split_t
        mult = in_bear ? bear_vol_mult : 1.0
        shift = in_bear ? bear_mu_shift : 0.0

        z_f = Array.new(n_factors) { rng.next_normal }
        f_raw = LinAlg.mat_vec(chol_f, z_f)
        f_t = Array.new(n_factors) { |k| mu_f[k] + (k.zero? ? shift : 0.0) + (mult * f_raw[k]) }

        r_t = Array.new(n_assets) do |i|
          val = alpha[i] + (mult * idio_sd[i] * rng.next_normal)
          bi = betas[i]
          n_factors.times { |k| val += bi[k] * f_t[k] }
          val
        end

        factors[t] = f_t
        returns[t] = r_t
      end

      { factors: factors, returns: returns, split_t: split_t }
    end

    def draw_innov(rng, df)
      df.nil? ? rng.next_normal : rng.next_std_t(df)
    end

    def build_betas(n_assets, n_factors, rng)
      Array.new(n_assets) do |i|
        Array.new(n_factors) do |k|
          if k.zero?
            0.6 + (0.8 * rng.next_float)
          else
            angle = (2.0 * Math::PI * (i + 1) * k) / n_assets
            (0.6 * Math.sin(angle)) + (0.25 * rng.next_normal)
          end
        end
      end
    end

    def build_factor_cov(n_factors, corr)
      vols = Array.new(n_factors) { |k| 0.045 - (0.008 * [k, 3].min) }
      Array.new(n_factors) do |i|
        Array.new(n_factors) do |j|
          c = i == j ? 1.0 : corr
          c * vols[i] * vols[j]
        end
      end
    end

    def build_alpha(n_assets, alpha_scale, rng)
      return Array.new(n_assets, 0.0) if alpha_scale <= 0.0

      raw = Array.new(n_assets) { rng.next_normal }
      mean = raw.sum / n_assets
      raw.map { |v| alpha_scale * (v - mean) }
    end

    def compute_total_cov(betas, sigma_f, idio_var)
      bf = LinAlg.mat_mul(betas, sigma_f)
      bfb = LinAlg.mat_mul(bf, LinAlg.transpose(betas))
      idio_var.each_with_index { |dv, i| bfb[i][i] += dv }
      bfb
    end

    def compute_expected_returns(alpha, betas, mu_f)
      betas.each_with_index.map do |bi, i|
        val = alpha[i]
        mu_f.each_with_index { |mf, k| val += bi[k] * mf }
        val
      end
    end
  end
end
