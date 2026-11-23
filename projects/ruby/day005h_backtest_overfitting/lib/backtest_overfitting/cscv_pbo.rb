# frozen_string_literal: true

module BacktestOverfitting
  module CscvPbo
    S_BLOCKS = 16

    module_function

    def half_masks
      @half_masks ||= begin
        masks = []
        (0...S_BLOCKS).to_a.combination(S_BLOCKS / 2) do |comb|
          lo = 0
          hi = 0
          comb.each { |b| b < 8 ? (lo |= (1 << b)) : (hi |= (1 << (b - 8))) }
          masks << [lo, hi]
        end
        masks.freeze
      end
    end

    def precompute_half_tables(mat, n_blocks: S_BLOCKS)
      n_strats = mat.length
      t_len = (mat.first.length / n_blocks) * n_blocks
      blk_sz = t_len / n_blocks
      b_s1, b_s2, tot_s1, tot_s2 = compute_block_sums(mat, n_strats, n_blocks, blk_sz)
      lo_s1, lo_s2 = build_256_half(b_s1[0, 8], b_s2[0, 8], n_strats)
      hi_s1, hi_s2 = build_256_half(b_s1[8, 8], b_s2[8, 8], n_strats)
      { n_strats: n_strats, half_len: t_len / 2, tot_s1: tot_s1, tot_s2: tot_s2,
        lo_s1: lo_s1, lo_s2: lo_s2, hi_s1: hi_s1, hi_s2: hi_s2 }
    end

    def compute_block_sums(mat, n_strats, n_blocks, blk_sz)
      b_s1 = Array.new(n_blocks) { Array.new(n_strats, 0.0) }
      b_s2 = Array.new(n_blocks) { Array.new(n_strats, 0.0) }
      tot_s1 = Array.new(n_strats, 0.0)
      tot_s2 = Array.new(n_strats, 0.0)
      n_strats.times do |j|
        col = mat[j]
        n_blocks.times do |b|
          s1 = 0.0
          s2 = 0.0
          offset = b * blk_sz
          blk_sz.times do |k|
            v = col[offset + k]
            s1 += v
            s2 += v * v
          end
          b_s1[b][j] = s1
          b_s2[b][j] = s2
          tot_s1[j] += s1
          tot_s2[j] += s2
        end
      end
      [b_s1, b_s2, tot_s1, tot_s2]
    end
    private_class_method :compute_block_sums

    def build_256_half(sub_s1, sub_s2, n_strats)
      out1 = Array.new(256) { Array.new(n_strats, 0.0) }
      out2 = Array.new(256) { Array.new(n_strats, 0.0) }
      (1...256).each do |mask|
        lsb = mask & -mask
        bit = lsb.bit_length - 1
        prev = mask ^ lsb
        p1 = out1[prev]
        p2 = out2[prev]
        add1 = sub_s1[bit]
        add2 = sub_s2[bit]
        c1 = out1[mask]
        c2 = out2[mask]
        n_strats.times do |j|
          c1[j] = p1[j] + add1[j]
          c2[j] = p2[j] + add2[j]
        end
      end
      [out1, out2]
    end
    private_class_method :build_256_half

    def eval_split_mask(tbl, lo_m, hi_m)
      n = tbl[:n_strats]
      hl = tbl[:half_len]
      l1 = tbl[:lo_s1][lo_m]
      l2 = tbl[:lo_s2][lo_m]
      h1 = tbl[:hi_s1][hi_m]
      h2 = tbl[:hi_s2][hi_m]
      ts1 = tbl[:tot_s1]
      ts2 = tbl[:tot_s2]

      best_j = 0
      best_is = -1.0e18
      oos_sr = Array.new(n, 0.0)

      n.times do |j|
        s1 = l1[j] + h1[j]
        s2 = l2[j] + h2[j]
        sr_is = sharpe_from_sums(s1, s2, hl)
        if sr_is > best_is
          best_is = sr_is
          best_j = j
        end
        oos_sr[j] = sharpe_from_sums(ts1[j] - s1, ts2[j] - s2, hl)
      end

      evaluate_oos_outcome(best_j, best_is, oos_sr, n)
    end

    def sharpe_from_sums(s1, s2, len)
      m = s1 / len
      var = (s2 - (len * m * m)) / (len - 1)
      var > 1.0e-24 ? m / Math.sqrt(var) : 0.0
    end
    private_class_method :sharpe_from_sums

    def evaluate_oos_outcome(best_j, best_is, oos_sr, n_strats)
      target = oos_sr[best_j]
      below = 0
      ties = 0
      n_strats.times do |j|
        v = oos_sr[j]
        if v < target
          below += 1
        elsif v == target
          ties += 1
        end
      end
      rank = below + (0.5 * (ties + 1))
      w_bar = rank / (n_strats + 1.0)
      lambda = Math.log(w_bar / (1.0 - w_bar))
      [lambda, best_is, target]
    end
    private_class_method :evaluate_oos_outcome

    def evaluate_cscv(mat, workers: nil, record_splits: false)
      tbl = precompute_half_tables(mat)
      masks = half_masks
      chunks = Parallel.map_chunks(n_items: masks.length, workers: workers, chunk_size: 1_287) do |cidx, len|
        start_i = cidx * 1_287
        Array.new(len) { |k| eval_split_mask(tbl, *masks[start_i + k]) }
      end.flatten(1)

      summarize_cscv_chunks(chunks, masks.length, record_splits)
    end

    def summarize_cscv_chunks(chunks, n_splits, record_splits)
      lambdas = chunks.map { |c| c[0] }
      sr_is = chunks.map { |c| c[1] }
      sr_oos = chunks.map { |c| c[2] }
      pbo = lambdas.count { |l| l <= 0.0 }.to_f / n_splits
      fit = Special.ols_fit(sr_is, sr_oos)
      out = {
        n_splits: n_splits, pbo: pbo, mean_lambda: lambdas.sum(0.0) / n_splits,
        mean_sr_is: sr_is.sum(0.0) / n_splits, mean_sr_oos: sr_oos.sum(0.0) / n_splits,
        degradation_intercept: fit[:intercept], degradation_slope: fit[:slope]
      }
      if record_splits
        out[:lambdas_sample] = lambdas.each_slice(10).map(&:first)
        out[:scatter_sample] = sr_is.zip(sr_oos).each_slice(30).map(&:first)
      end
      out
    end
    private_class_method :summarize_cscv_chunks

    def run_study(seed: MAIN_SEED, workers: nil)
      noise_anchor = noise_exchangeable_pbo_study(n_strats: 20, n_obs: 512, n_mats: 24, seed: seed, workers: workers)
      planted = planted_signal_pbo_study(seed: seed, workers: workers)
      scaling = pbo_scaling_study(seed: seed, workers: workers)
      costs = transaction_cost_study(seed: seed)
      { noise_anchor: noise_anchor, planted_signal: planted, scaling: scaling, cost_study: costs }
    end

    def noise_exchangeable_pbo_study(n_strats:, n_obs:, n_mats:, seed:, workers:)
      pbos = Array.new(n_mats) do |m_idx|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 800, m_idx))
        mat = Array.new(n_strats) { DataGenerator.iid_gaussian(n_obs, sigma: 0.01, rng: rng) }
        evaluate_cscv(mat, workers: workers, record_splits: m_idx.zero?)
      end
      vals = pbos.map { |r| r[:pbo] }
      mean = vals.sum(0.0) / n_mats
      se = Math.sqrt((vals.sum(0.0) { |v| (v - mean)**2 } / (n_mats - 1)) / n_mats)
      { n_strats: n_strats, n_obs: n_obs, n_mats: n_mats, exact_theoretical_pbo: 0.5,
        emp_mean_pbo: mean, emp_se_pbo: se, first_mat: pbos.first }
    end

    def planted_signal_pbo_study(seed:, workers:)
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 830, 0))
      mat = build_planted_matrix(n_strats: 60, n_true: 6, n_obs: 1_024, mu_sig: 0.0020, rng: rng)
      evaluate_cscv(mat, workers: workers, record_splits: true)
    end

    def pbo_scaling_study(seed:, workers:)
      vs_n = [20, 60, 150, 300].map do |n|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 850 + n, 0))
        mat = Array.new(n) { DataGenerator.iid_gaussian(512, sigma: 0.01, rng: rng) }
        res = evaluate_cscv(mat, workers: workers)
        { n_strats: n, pbo: res[:pbo], mean_sr_is: res[:mean_sr_is], mean_sr_oos: res[:mean_sr_oos] }
      end
      vs_t = [256, 512, 1_024].map do |t_len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 880 + t_len, 0))
        mat = build_planted_matrix(n_strats: 60, n_true: 6, n_obs: t_len, mu_sig: 0.0014, rng: rng)
        { n_obs: t_len, pbo: evaluate_cscv(mat, workers: workers)[:pbo] }
      end
      { vs_n_noise: vs_n, vs_t_planted: vs_t }
    end

    def build_planted_matrix(n_strats:, n_true:, n_obs:, mu_sig:, rng:)
      Array.new(n_strats) do |j|
        DataGenerator.iid_gaussian(n_obs, sigma: 0.01, mu: j < n_true ? mu_sig : 0.0, rng: rng)
      end
    end
    private_class_method :build_planted_matrix

    def transaction_cost_study(seed:)
      rng = Rng::Sampler.new(Rng.derive_seed(seed, 900, 0))
      rets = DataGenerator.ar1_gaussian(1_024, phi: 0.14, sigma: 0.01, rng: rng)
      lib = Strategy.build_library
      rows = [0.0, 5.0, 20.0].map do |bps|
        evals = BacktestEngine.evaluate_library(lib, rets, cost_bps: bps)
        hi_to = evals.select { |e| (e[:total_turnover] / rets.length) >= 0.40 }
        lo_to = evals.select { |e| (e[:total_turnover] / rets.length) < 0.08 }
        best = evals.max_by { |e| e[:sharpe] }
        { cost_bps: bps, best_name: best[:name], best_sharpe: best[:sharpe],
          hi_turnover_mean_sharpe: hi_to.sum(0.0) { |e| e[:sharpe] } / hi_to.length,
          lo_turnover_mean_sharpe: lo_to.sum(0.0) { |e| e[:sharpe] } / lo_to.length }
      end
      { n_obs: rets.length, phi: 0.14, rows: rows }
    end
  end
end
