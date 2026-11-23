# frozen_string_literal: true

module MlmcHeston
  module Mlmc
    module_function

    def ols_slope(xs, ys)
      n = xs.length
      mx = xs.sum(0.0) / n
      my = ys.sum(0.0) / n
      num = 0.0
      den = 0.0
      n.times do |i|
        dx = xs[i] - mx
        num += dx * (ys[i] - my)
        den += dx * dx
      end
      den.positive? ? num / den : 0.0
    end

    def level_stats(model:, params:, level:, n_paths:, seed: MAIN_SEED, workers: nil, n0: 4, offset: 0)
      strike = params.fetch(:strike, 100.0)
      maturity = params.fetch(:maturity, 1.0)
      scheme = params.fetch(:scheme, :euler)
      payoff_type = params.fetch(:payoff_type, :european)
      stream = (model == :gbm ? 100 : 200) + (scheme == :milstein ? 50 : 0) + (level * 7) + (offset * 1000)
      chunks = Parallel.map_chunks(n_items: n_paths, workers: workers, chunk_size: 2_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, stream, cidx))
        s_pf = 0.0
        s_y1 = 0.0
        s_y2 = 0.0
        len.times do
          pf, y = if model == :gbm
                    Simulators.gbm_coupled_pair(params, strike, maturity, level, scheme: scheme, rng: rng, n0: n0)
                  else
                    Simulators.heston_coupled_pair(params, strike, maturity, level,
                                                   scheme: scheme, payoff_type: payoff_type, rng: rng, n0: n0)
                  end
          s_pf += pf
          s_y1 += y
          s_y2 += y * y
        end
        [s_pf, s_y1, s_y2]
      end

      sum_pf = chunks.sum(0.0) { |c| c[0] }
      sum_y1 = chunks.sum(0.0) { |c| c[1] }
      sum_y2 = chunks.sum(0.0) { |c| c[2] }
      mean_y = sum_y1 / n_paths
      var_y = [(sum_y2 / n_paths) - (mean_y * mean_y), 1.0e-18].max
      { level: level, n_paths: n_paths, sum_pf: sum_pf, sum_y1: sum_y1, sum_y2: sum_y2,
        mean_y: mean_y, var_y: var_y }
    end

    def estimate_rates(model:, params:, levels: [1, 2, 3, 4, 5], n_probe: 16_000, seed: MAIN_SEED, workers: nil)
      rows = levels.map do |l|
        st = level_stats(model: model, params: params, level: l, n_paths: n_probe, seed: seed, workers: workers)
        { level: l, mean_abs: st[:mean_y].abs, var: st[:var_y], log2_mean: Math.log2([st[:mean_y].abs, 1.0e-16].max),
          log2_var: Math.log2(st[:var_y]) }
      end
      ls = rows.map { |r| r[:level].to_f }
      alpha = -ols_slope(ls, rows.map { |r| r[:log2_mean] })
      beta = -ols_slope(ls, rows.map { |r| r[:log2_var] })
      { rows: rows, alpha_weak: alpha, beta_strong: beta }
    end

    def run_adaptive(model:, params:, eps:, seed: MAIN_SEED, workers: nil, n0: 4, max_level: 7)
      l_max = 2
      n_init = 3_000
      accum = {}
      (0..l_max).each do |l|
        accum[l] = level_stats(model: model, params: params, level: l, n_paths: n_init,
                               seed: seed, workers: workers, n0: n0, offset: 0)
      end

      loop do
        update_optimal_samples(accum, l_max, eps, n0, model, params, seed, workers)
        rem_bias = [accum[l_max][:mean_y].abs, 0.5 * accum[l_max - 1][:mean_y].abs].max
        break if rem_bias <= eps / Math.sqrt(2.0) || l_max >= max_level

        l_max += 1
        accum[l_max] = level_stats(model: model, params: params, level: l_max, n_paths: n_init,
                                   seed: seed, workers: workers, n0: n0, offset: 0)
      end

      build_mlmc_result(accum, l_max, n0, eps, true)
    end

    def update_optimal_samples(accum, l_max, eps, n0, model, params, seed, workers)
      costs = (0..l_max).map { |l| (l.zero? ? 1.0 : 1.5) * n0 * (2**l) }
      vars = (0..l_max).map { |l| accum[l][:var_y] }
      sum_vc = (0..l_max).sum(0.0) { |l| Math.sqrt(vars[l] * costs[l]) }
      mult = 2.0 * sum_vc / (eps * eps)

      (0..l_max).each do |l|
        n_opt = (mult * Math.sqrt(vars[l] / costs[l])).ceil
        extra = n_opt - accum[l][:n_paths]
        next if extra <= 0

        add = level_stats(model: model, params: params, level: l, n_paths: extra,
                          seed: seed, workers: workers, n0: n0, offset: accum[l][:n_paths])
        merge_accum!(accum[l], add)
      end
    end
    private_class_method :update_optimal_samples

    def merge_accum!(acc, add)
      acc[:n_paths] += add[:n_paths]
      acc[:sum_y1] += add[:sum_y1]
      acc[:sum_y2] += add[:sum_y2]
      acc[:mean_y] = acc[:sum_y1] / acc[:n_paths]
      acc[:var_y] = [(acc[:sum_y2] / acc[:n_paths]) - (acc[:mean_y] * acc[:mean_y]), 1.0e-18].max
    end
    private_class_method :merge_accum!

    def build_mlmc_result(accum, l_max, n0, eps, _ok)
      price = (0..l_max).sum(0.0) { |l| accum[l][:mean_y] }
      mc_var = (0..l_max).sum(0.0) { |l| accum[l][:var_y] / accum[l][:n_paths] }
      total_cost = (0..l_max).sum(0.0) { |l| accum[l][:n_paths] * (l.zero? ? 1.0 : 1.5) * n0 * (2**l) }
      {
        eps: eps, price: price, se: Math.sqrt(mc_var), max_level: l_max, total_cost: total_cost,
        n_samples: (0..l_max).map { |l| accum[l][:n_paths] },
        vars: (0..l_max).map { |l| accum[l][:var_y] }
      }
    end
    private_class_method :build_mlmc_result
  end
end
