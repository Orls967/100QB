# frozen_string_literal: true

module FactorPortfolioLab
  # Uniform Strategy interface for portfolio construction in rolling backtests.
  module Strategy
    class Base
      attr_reader :name

      def initialize(name)
        @name = name
      end

      def weights(_returns_window, factors_window: nil)
        raise NotImplementedError
      end
    end

    class EqualWeight < Base
      def initialize = super('EqualWeight')

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        n = returns_window[0].length
        Array.new(n, 1.0 / n)
      end
    end

    class MinVarSample < Base
      def initialize = super('MinVarSample')

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        cov = regularize_diag(LinAlg.cov_matrix(returns_window))
        Markowitz.solve_qp(cov, kappa: 0.0)[:weights]
      end
    end

    class MinVarShrink < Base
      def initialize = super('MinVarShrink')

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        lw = Shrinkage.shrink_cov_lw(returns_window)
        Markowitz.solve_qp(lw[:sigma], kappa: 0.0)[:weights]
      end
    end

    class TangencyShrink < Base
      def initialize(kappa: 0.4)
        super('TangencyShrink')
        @kappa = kappa
      end

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        n = returns_window[0].length
        lw = Shrinkage.shrink_cov_lw(returns_window)
        mu_hat = Pca.column_means(returns_window, n)
        mu_js = Shrinkage.shrink_mu_js(mu_hat, lw[:sigma], returns_window.length)
        Markowitz.solve_qp(lw[:sigma], mu: mu_js, kappa: @kappa)[:weights]
      end
    end

    class RiskParitySample < Base
      def initialize = super('RiskParitySample')

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        lw = Shrinkage.shrink_cov_lw(returns_window)
        RiskParity.solve_newton(lw[:sigma])[:weights]
      end
    end

    class RiskParityFactor < Base
      def initialize(n_factors: 3)
        super('RiskParityFactor')
        @n_factors = n_factors
      end

      def weights(returns_window, factors_window: nil)
        cov = if factors_window
                factor_cov_estimate(returns_window, factors_window)
              else
                pca_denoised_cov(returns_window, @n_factors)
              end
        RiskParity.solve_newton(cov)[:weights]
      end
    end

    class BlackLittermanView < Base
      def initialize(kappa: 0.35)
        super('BlackLittermanView')
        @kappa = kappa
      end

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        n = returns_window[0].length
        lw = Shrinkage.shrink_cov_lw(returns_window)
        sigma = lw[:sigma]
        w_mkt = Array.new(n, 1.0 / n)
        pi = BlackLitterman.implied_equilibrium(sigma, w_mkt, risk_aversion: 2.5)
        mu_hat = Pca.column_means(returns_window, n)
        order = (0...n).sort_by { |i| -mu_hat[i] }
        p_row = Array.new(n, 0.0)
        p_row[order[0]] = 1.0
        p_row[order[-1]] = -1.0
        p_mat = [p_row]
        q_vec = [0.004]
        omega = BlackLitterman.idzorek_omega(sigma, p_mat, [0.35], tau: 0.05)
        post = BlackLitterman.posterior(sigma: sigma, pi: pi, p_mat: p_mat, q_vec: q_vec, omega: omega, tau: 0.05)
        Markowitz.solve_qp(post[:sigma_bl], mu: post[:mu_post], kappa: @kappa)[:weights]
      end
    end

    class PcaFactorMinVar < Base
      def initialize(n_factors: 3)
        super('PcaFactorMinVar')
        @n_factors = n_factors
      end

      def weights(returns_window, factors_window: nil)
        _ = factors_window
        cov = pca_denoised_cov(returns_window, @n_factors)
        Markowitz.solve_qp(cov, kappa: 0.0)[:weights]
      end
    end

    module_function

    def all_strategies
      [
        EqualWeight.new,
        MinVarSample.new,
        MinVarShrink.new,
        TangencyShrink.new,
        RiskParitySample.new,
        RiskParityFactor.new,
        BlackLittermanView.new,
        PcaFactorMinVar.new
      ]
    end

    def regularize_diag(cov, eps: 1.0e-7)
      n = cov.length
      Array.new(n) { |i| Array.new(n) { |j| cov[i][j] + (i == j ? eps : 0.0) } }
    end

    def factor_cov_estimate(returns_window, factors_window)
      ts = TimeSeriesRegression.fit(returns_window, factors_window)
      betas = ts[:betas]
      k = factors_window[0].length
      f_bar = Pca.column_means(factors_window, k)
      sigma_f = Grs.sample_cov_mle(factors_window, f_bar, k)
      idio_var = ts[:sigma_eps_unbiased].each_with_index.map { |row, i| [row[i], 1.0e-6].max }
      FactorWorld.compute_total_cov(betas, sigma_f, idio_var)
    end

    def pca_denoised_cov(returns_window, k_factors)
      fit = Pca.fit(returns_window)
      n = fit[:n_assets]
      k = [k_factors, n - 1].min
      vals = fit[:eigenvalues]
      vecs = fit[:eigenvectors]
      bulk_avg = vals[k...n].sum / [n - k, 1].max
      cov = Array.new(n) { Array.new(n, 0.0) }
      n.times do |i|
        (i...n).each do |j|
          val = (0...k).sum { |m| vals[m] * vecs[i][m] * vecs[j][m] }
          val += bulk_avg if i == j
          cov[i][j] = val
          cov[j][i] = val
        end
      end
      cov
    end
  end

  class Strategy::Base # rubocop:disable Style/ClassAndModuleChildren
    def regularize_diag(cov) = Strategy.regularize_diag(cov)
    def factor_cov_estimate(rw, fw) = Strategy.factor_cov_estimate(rw, fw)
    def pca_denoised_cov(rw, k) = Strategy.pca_denoised_cov(rw, k)
  end
end
