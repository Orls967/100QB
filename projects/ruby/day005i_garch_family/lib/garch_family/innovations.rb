# frozen_string_literal: true

module GarchFamily
  module Innovations
    class Normal
      def name = 'normal'

      def n_dist_params = 0

      def sample(rng, _nu = nil) = rng.next_normal

      def expected_abs_z(_nu = nil) = Special::SQRT_2_OVER_PI

      def log_pdf(z, _nu = nil) = -0.5 * (Math.log(2.0 * Math::PI) + (z * z))

      def obs_loglik(eps, var, _nu = nil)
        -0.5 * (Math.log(2.0 * Math::PI) + Math.log(var) + ((eps * eps) / var))
      end

      def score_wrt_var(eps, var, _nu = nil) = 0.5 * (((eps * eps) / var) - 1.0) / var

      def score_wrt_mu_direct(eps, var, _nu = nil) = eps / var

      def score_wrt_nu(_eps, _var, _nu = nil) = 0.0

      def var_quantile(alpha, _nu = nil) = Special.norm_inv(1.0 - alpha)

      def expected_shortfall(alpha, _nu = nil)
        q = var_quantile(alpha)
        Special.norm_pdf(q) / alpha
      end
    end

    class StudentT
      attr_reader :nu

      def initialize(nu: 5.0)
        @nu = nu.to_f
      end

      def name = "student_t(#{@nu.round(2)})"

      def n_dist_params = 1

      def sample(rng, nu = @nu) = rng.next_std_t(nu)

      def expected_abs_z(nu = @nu)
        2.0 * Math.sqrt((nu - 2.0) / Math::PI) *
          Math.exp(Special.lgamma(0.5 * (nu + 1.0)) - Special.lgamma(0.5 * nu)) / (nu - 1.0)
      end

      def log_const(nu = @nu)
        Special.lgamma(0.5 * (nu + 1.0)) - Special.lgamma(0.5 * nu) - (0.5 * Math.log(Math::PI * (nu - 2.0)))
      end

      def log_pdf(z, nu = @nu)
        log_const(nu) - (0.5 * (nu + 1.0) * Math.log1p((z * z) / (nu - 2.0)))
      end

      def obs_loglik(eps, var, nu = @nu)
        z2 = (eps * eps) / var
        log_const(nu) - (0.5 * Math.log(var)) - (0.5 * (nu + 1.0) * Math.log1p(z2 / (nu - 2.0)))
      end

      def weight(eps, var, nu = @nu) = (nu + 1.0) / ((nu - 2.0) + ((eps * eps) / var))

      def score_wrt_var(eps, var, nu = @nu)
        z2 = (eps * eps) / var
        0.5 * ((weight(eps, var, nu) * z2) - 1.0) / var
      end

      def score_wrt_mu_direct(eps, var, nu = @nu) = weight(eps, var, nu) * eps / var

      def score_wrt_nu(eps, var, nu = @nu)
        z2 = (eps * eps) / var
        dconst = 0.5 * (Special.digamma(0.5 * (nu + 1.0)) - Special.digamma(0.5 * nu) - (1.0 / (nu - 2.0)))
        dlog = -0.5 * Math.log1p(z2 / (nu - 2.0))
        dinside = 0.5 * (nu + 1.0) * z2 / ((nu - 2.0) * ((nu - 2.0) + z2))
        dconst + dlog + dinside
      end

      def var_quantile(alpha, nu = @nu) = -Special.std_t_inv(alpha, nu)

      def expected_shortfall(alpha, nu = @nu)
        t_q = -Special.student_t_inv(alpha, nu)
        Math.sqrt((nu - 2.0) / nu) * (Special.student_t_pdf(t_q, nu) / alpha) * ((nu + (t_q * t_q)) / (nu - 1.0))
      end
    end

    module SkewT
      module_function

      def constants(nu, lam)
        c = Math.exp(Special.lgamma(0.5 * (nu + 1.0)) - Special.lgamma(0.5 * nu)) / Math.sqrt(Math::PI * (nu - 2.0))
        a = 4.0 * lam * c * ((nu - 2.0) / (nu - 1.0))
        b = Math.sqrt(1.0 + (3.0 * lam * lam) - (a * a))
        [a, b, c]
      end

      def log_pdf(z, nu, lam)
        a, b, c = constants(nu, lam)
        arg = (b * z) + a
        s = arg >= 0.0 ? (1.0 + lam) : (1.0 - lam)
        Math.log(b) + Math.log(c) - (0.5 * (nu + 1.0) * Math.log1p(((arg / s)**2) / (nu - 2.0)))
      end

      def inv_cdf(prob, nu, lam)
        a, b, = constants(nu, lam)
        thresh = 0.5 * (1.0 - lam)
        raw = if prob < thresh
                (1.0 - lam) * Special.std_t_inv(prob / (1.0 - lam), nu)
              else
                (1.0 + lam) * Special.std_t_inv(0.5 + ((prob - thresh) / (1.0 + lam)), nu)
              end
        (raw - a) / b
      end
    end
  end
end
