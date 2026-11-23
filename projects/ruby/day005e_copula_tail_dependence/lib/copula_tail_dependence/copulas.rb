# frozen_string_literal: true

module CopulaTailDependence
  module Copulas
    module_function

    def sample_pair(type, params, rng)
      case type
      when :gaussian
        rho = params.fetch(:rho)
        z1 = rng.next_normal
        z2 = (rho * z1) + (Math.sqrt(1.0 - (rho * rho)) * rng.next_normal)
        [Special.norm_cdf(z1), Special.norm_cdf(z2)]
      when :t
        rho = params.fetch(:rho)
        df = params.fetch(:df)
        z1 = rng.next_normal
        z2 = (rho * z1) + (Math.sqrt(1.0 - (rho * rho)) * rng.next_normal)
        scale = Math.sqrt(df / rng.next_chi2(df))
        [Special.t_cdf(z1 * scale, df), Special.t_cdf(z2 * scale, df)]
      when :clayton
        theta = params.fetch(:theta)
        v = rng.next_gamma(1.0 / theta, 1.0)
        e1 = rng.next_exp(1.0)
        e2 = rng.next_exp(1.0)
        [(1.0 + (e1 / v))**(-1.0 / theta), (1.0 + (e2 / v))**(-1.0 / theta)]
      when :gumbel
        theta = params.fetch(:theta)
        inv_t = 1.0 / theta
        s = rng.next_positive_stable(inv_t)
        e1 = rng.next_exp(1.0)
        e2 = rng.next_exp(1.0)
        [Math.exp(-((e1 / s)**inv_t)), Math.exp(-((e2 / s)**inv_t))]
      else
        raise ArgumentError, "Kopula tidak dikenal: #{type}"
      end
    end

    def exact_kendall_tau(type, params)
      case type
      when :gaussian, :t
        (2.0 / Math::PI) * Math.asin(params.fetch(:rho))
      when :clayton
        th = params.fetch(:theta)
        th / (th + 2.0)
      when :gumbel
        th = params.fetch(:theta)
        1.0 - (1.0 / th)
      else
        raise ArgumentError, "Kopula tidak dikenal: #{type}"
      end
    end

    def exact_spearman_rho(type, params)
      return nil unless type == :gaussian

      (6.0 / Math::PI) * Math.asin(0.5 * params.fetch(:rho))
    end

    def exact_lower_tail(type, params)
      case type
      when :gaussian, :gumbel
        0.0
      when :t
        rho = params.fetch(:rho)
        df = params.fetch(:df)
        arg = -Math.sqrt((df + 1.0) * (1.0 - rho) / (1.0 + rho))
        2.0 * Special.t_cdf(arg, df + 1.0)
      when :clayton
        2.0**(-1.0 / params.fetch(:theta))
      else
        raise ArgumentError, "Kopula tidak dikenal: #{type}"
      end
    end

    def exact_upper_tail(type, params)
      case type
      when :gaussian, :clayton
        0.0
      when :t
        exact_lower_tail(:t, params)
      when :gumbel
        2.0 - (2.0**(1.0 / params.fetch(:theta)))
      else
        raise ArgumentError, "Kopula tidak dikenal: #{type}"
      end
    end

    # Nilai eksak chi(u) = P(U1 > u, U2 > u) / (1 - u) untuk kopula Gaussian pada u berhingga.
    def exact_gaussian_chi(u, rho)
      z_u = Special.norm_inv(u)
      Special.bivariate_norm_cdf(-z_u, -z_u, rho) / (1.0 - u)
    end

    def gaussian_chi_grid(rho, u_grid = [0.9, 0.99, 0.999, 1.0 - 1.0e-4, 1.0 - 1.0e-6, 1.0 - 1.0e-8, 1.0 - 1.0e-10])
      u_grid.map { |u| { u: u, chi: exact_gaussian_chi(u, rho) } }
    end

    # Log-densitas kopula bivariat untuk MLE (Bagian 3).
    def log_density_bivariate(type, u1, u2, params)
      u = u1.clamp(1.0e-12, 1.0 - 1.0e-12)
      v = u2.clamp(1.0e-12, 1.0 - 1.0e-12)
      case type
      when :gaussian
        log_density_gaussian(u, v, params.fetch(:rho))
      when :t
        log_density_t(u, v, params.fetch(:rho), params.fetch(:df))
      when :clayton
        th = params.fetch(:theta)
        Math.log(1.0 + th) - ((1.0 + th) * (Math.log(u) + Math.log(v))) -
          ((2.0 + (1.0 / th)) * Math.log((u**(-th)) + (v**(-th)) - 1.0))
      else
        raise ArgumentError, "Log-densitas belum didukung untuk #{type}"
      end
    end

    def log_density_gaussian(u, v, rho)
      z1 = Special.norm_inv(u)
      z2 = Special.norm_inv(v)
      r2 = rho * rho
      (-0.5 * Math.log(1.0 - r2)) - (((r2 * ((z1 * z1) + (z2 * z2))) - (2.0 * rho * z1 * z2)) / (2.0 * (1.0 - r2)))
    end
    private_class_method :log_density_gaussian

    def log_density_t(u, v, rho, df)
      return log_density_gaussian(u, v, rho) if df.infinite?

      t1 = Special.t_inv(u, df)
      t2 = Special.t_inv(v, df)
      r2 = rho * rho
      q = ((t1 * t1) - (2.0 * rho * t1 * t2) + (t2 * t2)) / (df * (1.0 - r2))
      lg_half = Special.log_gamma(0.5 * (df + 1.0))
      c = Special.log_gamma(0.5 * (df + 2.0)) + Special.log_gamma(0.5 * df) - (2.0 * lg_half)
      c - (0.5 * Math.log(1.0 - r2)) - (0.5 * (df + 2.0) * Math.log(1.0 + q)) +
        (0.5 * (df + 1.0) * (Math.log(1.0 + ((t1 * t1) / df)) + Math.log(1.0 + ((t2 * t2) / df))))
    end
    private_class_method :log_density_t
  end
end
