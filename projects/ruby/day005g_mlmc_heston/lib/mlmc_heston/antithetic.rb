# frozen_string_literal: true

module MlmcHeston
  module Antithetic
    module_function

    # Pasangan antitetik Giles-Szpruch: menukar inkremen pasangan fine (dw_2k, dw_{2k+1})
    # sehingga P_f^{(1)} dan P_f^{(2)} memiliki rata-rata dengan variansi jauh lebih kecil terhadap P_c.
    def heston_antithetic_var(params, strike, maturity, level:, n_paths: 10_000, seed: MAIN_SEED, workers: nil, n0: 4)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 2_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 700 + level, cidx))
        s_std1 = 0.0
        s_std2 = 0.0
        s_ant1 = 0.0
        s_ant2 = 0.0
        len.times do
          y_std, y_ant = draw_antithetic_pair(params, strike, maturity, level, n0, rng)
          s_std1 += y_std
          s_std2 += y_std * y_std
          s_ant1 += y_ant
          s_ant2 += y_ant * y_ant
        end
        [s_std1, s_std2, s_ant1, s_ant2]
      end

      v_std = var_from_sums(chunks.sum(0.0) { |c| c[0] }, chunks.sum(0.0) { |c| c[1] }, n_paths)
      v_ant = var_from_sums(chunks.sum(0.0) { |c| c[2] }, chunks.sum(0.0) { |c| c[3] }, n_paths)
      { level: level, var_standard: v_std, var_antithetic: v_ant, reduction_ratio: v_std / [v_ant, 1.0e-18].max }
    end

    # Diagnostik tambahan: kopel antitetik GBM dengan skema Milstein vs Euler untuk menunjukkan bahwa
    # pertukaran inkremen butuh skema dengan koreksi kuadratik simetris agar variansi level anjlok tajam.
    def gbm_antithetic_check(params, strike, maturity, level:, scheme:, n_paths: 6_000, seed: MAIN_SEED,
                             workers: nil, n0: 4)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 2_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 850 + level + (scheme == :milstein ? 10 : 0), cidx))
        s = [0.0, 0.0, 0.0, 0.0]
        len.times do
          y_std, y_ant = draw_gbm_antithetic_pair(params, strike, maturity, level, n0, rng, scheme)
          s[0] += y_std
          s[1] += y_std * y_std
          s[2] += y_ant
          s[3] += y_ant * y_ant
        end
        s
      end
      v_std = var_from_sums(chunks.sum(0.0) { |c| c[0] }, chunks.sum(0.0) { |c| c[1] }, n_paths)
      v_ant = var_from_sums(chunks.sum(0.0) { |c| c[2] }, chunks.sum(0.0) { |c| c[3] }, n_paths)
      { level: level, scheme: scheme.to_s, var_standard: v_std, var_antithetic: v_ant,
        reduction_ratio: v_std / [v_ant, 1.0e-18].max }
    end

    def cir_antithetic_check(params, maturity, level:, n_paths: 6_000, seed: MAIN_SEED, workers: nil, n0: 4)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 2_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 900 + level, cidx))
        s = [0.0, 0.0, 0.0, 0.0]
        len.times do
          y_std, y_ant = draw_cir_antithetic_pair(params, maturity, level, n0, rng)
          s[0] += y_std
          s[1] += y_std * y_std
          s[2] += y_ant
          s[3] += y_ant * y_ant
        end
        s
      end
      v_std = var_from_sums(chunks.sum(0.0) { |c| c[0] }, chunks.sum(0.0) { |c| c[1] }, n_paths)
      v_ant = var_from_sums(chunks.sum(0.0) { |c| c[2] }, chunks.sum(0.0) { |c| c[3] }, n_paths)
      { level: level, var_standard: v_std, var_antithetic: v_ant, reduction_ratio: v_std / [v_ant, 1.0e-18].max }
    end

    def multi_model_level_table(gbm_params, heston_params, levels: [1, 2, 3, 4, 5, 6], n_paths: 6_000,
                                seed: MAIN_SEED, workers: nil)
      rho = heston_params.fetch(:rho)
      heston_bound = 2.0 / (1.0 + (rho * rho))
      levels.map do |lvl|
        g = gbm_antithetic_check(gbm_params, 100.0, 1.0, level: lvl, scheme: :euler,
                                                         n_paths: n_paths, seed: seed, workers: workers)
        c = cir_antithetic_check(heston_params, 1.0, level: lvl, n_paths: n_paths, seed: seed, workers: workers)
        h = heston_antithetic_var(heston_params, 100.0, 1.0,
                                  level: lvl, n_paths: n_paths, seed: seed + 99, workers: workers)
        { level: lvl, gbm_euler_ratio: g[:reduction_ratio], cir_euler_ratio: c[:reduction_ratio],
          heston_euler_ratio: h[:reduction_ratio], heston_theoretical_ceiling: heston_bound }
      end
    end

    def draw_cir_antithetic_pair(params, maturity, level, n0, rng)
      n_f = n0 * (2**level)
      dt_f = maturity / n_f
      sqrt_dt_f = Math.sqrt(dt_f)
      dw = Array.new(n_f) { sqrt_dt_f * rng.next_normal }
      pf1 = run_cir_path(params, dt_f, dw)
      pf2 = run_cir_path(params, dt_f, Array.new(n_f) { |i| dw[i ^ 1] })
      dw_c = Array.new(n_f / 2) { |j| dw[2 * j] + dw[(2 * j) + 1] }
      pc = run_cir_path(params, 2.0 * dt_f, dw_c)
      [pf1 - pc, (0.5 * (pf1 + pf2)) - pc]
    end
    private_class_method :draw_cir_antithetic_pair

    def run_cir_path(params, dt, dw)
      v = params.fetch(:v0)
      kappa = params.fetch(:kappa)
      theta = params.fetch(:theta)
      xi = params.fetch(:xi)
      dw.each do |d|
        vp = v.positive? ? v : 0.0
        v += (kappa * (theta - vp) * dt) + (xi * Math.sqrt(vp) * d)
      end
      v.positive? ? v : 0.0
    end
    private_class_method :run_cir_path

    def draw_antithetic_pair(params, strike, maturity, level, n0, rng)
      n_f = n0 * (2**level)
      dt_f = maturity / n_f
      sqrt_dt_f = Math.sqrt(dt_f)
      rho = params.fetch(:rho)
      s1r = Math.sqrt(1.0 - (rho * rho))

      dw1 = Array.new(n_f, 0.0)
      dw2 = Array.new(n_f, 0.0)
      n_f.times do |i|
        z1 = rng.next_normal
        z2 = (rho * z1) + (s1r * rng.next_normal)
        dw1[i] = sqrt_dt_f * z1
        dw2[i] = sqrt_dt_f * z2
      end

      pf1 = run_heston_path(params, strike, maturity, dt_f, dw1, dw2)
      # Tukar pasangan (2k, 2k+1) untuk jalur antitetik fine
      dw1_a = Array.new(n_f) { |i| dw1[i ^ 1] }
      dw2_a = Array.new(n_f) { |i| dw2[i ^ 1] }
      pf2 = run_heston_path(params, strike, maturity, dt_f, dw1_a, dw2_a)

      n_c = n_f / 2
      dw1_c = Array.new(n_c) { |j| dw1[2 * j] + dw1[(2 * j) + 1] }
      dw2_c = Array.new(n_c) { |j| dw2[2 * j] + dw2[(2 * j) + 1] }
      pc = run_heston_path(params, strike, maturity, 2.0 * dt_f, dw1_c, dw2_c)

      [pf1 - pc, (0.5 * (pf1 + pf2)) - pc]
    end
    private_class_method :draw_antithetic_pair

    def draw_gbm_antithetic_pair(params, strike, maturity, level, n0, rng, scheme)
      n_f = n0 * (2**level)
      dt_f = maturity / n_f
      sqrt_dt_f = Math.sqrt(dt_f)
      dw = Array.new(n_f) { sqrt_dt_f * rng.next_normal }
      pf1 = run_gbm_path(params, strike, maturity, dt_f, dw, scheme)
      pf2 = run_gbm_path(params, strike, maturity, dt_f, Array.new(n_f) { |i| dw[i ^ 1] }, scheme)
      dw_c = Array.new(n_f / 2) { |j| dw[2 * j] + dw[(2 * j) + 1] }
      pc = run_gbm_path(params, strike, maturity, 2.0 * dt_f, dw_c, scheme)
      [pf1 - pc, (0.5 * (pf1 + pf2)) - pc]
    end
    private_class_method :draw_gbm_antithetic_pair

    def run_gbm_path(params, strike, maturity, dt, dw, scheme)
      s = params.fetch(:spot)
      rate = params.fetch(:rate)
      sigma = params.fetch(:sigma)
      dw.each do |d|
        mil = scheme == :milstein ? 0.5 * sigma * sigma * s * ((d * d) - dt) : 0.0
        s += (rate * s * dt) + (sigma * s * d) + mil
      end
      Math.exp(-rate * maturity) * [s - strike, 0.0].max
    end
    private_class_method :run_gbm_path

    def run_heston_path(params, strike, maturity, dt, dw1, dw2)
      x = Math.log(params.fetch(:spot))
      v = params.fetch(:v0)
      rate = params.fetch(:rate)
      kappa = params.fetch(:kappa)
      theta = params.fetch(:theta)
      xi = params.fetch(:xi)
      i = 0
      n = dw1.length
      while i < n
        vp = v.positive? ? v : 0.0
        sv = Math.sqrt(vp)
        x += ((rate - (0.5 * vp)) * dt) + (sv * dw1[i])
        v += (kappa * (theta - vp) * dt) + (xi * sv * dw2[i])
        i += 1
      end
      Math.exp(-rate * maturity) * [Math.exp(x) - strike, 0.0].max
    end
    private_class_method :run_heston_path

    def var_from_sums(s1, s2, n)
      m = s1 / n
      [(s2 / n) - (m * m), 0.0].max
    end
    private_class_method :var_from_sums
  end
end
