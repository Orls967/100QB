# frozen_string_literal: true

module MlmcHeston
  module HestonCf
    module_function

    # Fungsi karakteristik Heston bentuk "little trap" (Albrecher et al. 2007)
    # untuk x_T = ln(S_T / S_0).
    def char_fn(u_c, params, maturity)
      kappa = params.fetch(:kappa)
      theta = params.fetch(:theta)
      xi = params.fetch(:xi)
      rho = params.fetch(:rho)
      v0 = params.fetch(:v0)
      rate = params.fetch(:rate)

      return bs_char_fn(u_c, params, maturity) if xi.abs < 1.0e-12

      xi2 = xi * xi
      beta = kappa - (Complex(0.0, 1.0) * rho * xi * u_c)
      d = Special.csqrt((beta * beta) + (xi2 * ((Complex(0.0, 1.0) * u_c) + (u_c * u_c))))
      g = (beta - d) / (beta + d)
      exp_dt = Special.cexp(-d * maturity)
      one_minus_g_exp = 1.0 - (g * exp_dt)
      c_term = (Complex(0.0, 1.0) * u_c * rate * maturity) +
               ((kappa * theta / xi2) * (((beta - d) * maturity) - (2.0 * Special.clog(one_minus_g_exp / (1.0 - g)))))
      d_term = ((beta - d) / xi2) * ((1.0 - exp_dt) / one_minus_g_exp)
      Special.cexp(c_term + (d_term * v0))
    end

    def bs_char_fn(u_c, params, maturity)
      kappa = params.fetch(:kappa)
      theta = params.fetch(:theta)
      v0 = params.fetch(:v0)
      rate = params.fetch(:rate)
      v_bar = theta + ((v0 - theta) * (1.0 - Math.exp(-kappa * maturity)) / (kappa * maturity))
      Special.cexp((Complex(0.0, 1.0) * u_c * (rate - (0.5 * v_bar)) * maturity) - (0.5 * v_bar * maturity * u_c * u_c))
    end
    private_class_method :bs_char_fn

    # Harga call Eropa Heston via inversi Gil-Pelaez (Gauss-Legendre 64-titik pada [0, u_max]).
    def call_price(params, strike, maturity, u_max: 140.0, n_gl: 96)
      spot = params.fetch(:spot)
      rate = params.fetch(:rate)
      return [spot - strike, 0.0].max if maturity <= 0.0

      nodes, weights = Special.gauss_legendre(n_gl)
      half = 0.5 * u_max
      k_log = Math.log(strike / spot)
      i_unit = Complex(0.0, 1.0)
      phi_minus_i = char_fn(-i_unit, params, maturity)

      int1 = 0.0
      int2 = 0.0
      nodes.each_with_index do |xi_n, idx|
        u = half * (xi_n + 1.0)
        next if u <= 1.0e-12

        e_iuk = Special.cexp(-i_unit * u * k_log)
        phi_u = char_fn(u, params, maturity)
        phi_u_shift = char_fn(u - i_unit, params, maturity) / phi_minus_i
        int1 += weights[idx] * (e_iuk * phi_u_shift / (i_unit * u)).real
        int2 += weights[idx] * (e_iuk * phi_u / (i_unit * u)).real
      end

      p1 = 0.5 + (half * int1 / Math::PI)
      p2 = 0.5 + (half * int2 / Math::PI)
      (spot * p1) - (strike * Math.exp(-rate * maturity) * p2)
    end

    def put_price(params, strike, maturity, u_max: 140.0, n_gl: 96)
      call_price(params, strike, maturity, u_max: u_max, n_gl: n_gl) -
        params.fetch(:spot) + (strike * Math.exp(-params.fetch(:rate) * maturity))
    end
  end
end
