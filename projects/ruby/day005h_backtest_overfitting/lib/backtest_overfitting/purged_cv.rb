# frozen_string_literal: true

module BacktestOverfitting
  module PurgedCv
    module_function

    def generate_overlapping_dataset(n_obs:, horizon:, n_features:, rng:)
      rets = DataGenerator.iid_gaussian(n_obs + horizon, sigma: 0.01, rng: rng)
      indep_z = Array.new(n_features) { DataGenerator.iid_gaussian(n_obs + horizon, sigma: 1.0, rng: rng) }
      y_labels = Array.new(n_obs) { |t| rets[t, horizon].sum(0.0) }
      x_mat = build_rolling_features(n_obs, horizon, n_features, indep_z)
      [x_mat, y_labels]
    end

    def build_rolling_features(n_obs, horizon, n_features, indep_z)
      inv_sh = 1.0 / Math.sqrt(horizon)
      Array.new(n_obs) do |t|
        Array.new(n_features) { |f| indep_z[f][t, horizon].sum(0.0) * inv_sh }
      end
    end
    private_class_method :build_rolling_features

    def split_folds_naive(n_obs, k_folds, rng)
      perm = (0...n_obs).to_a
      (n_obs - 1).downto(1) do |i|
        j = rng.next_int(i + 1)
        perm[i], perm[j] = perm[j], perm[i]
      end
      fold_sz = n_obs / k_folds
      Array.new(k_folds) { |k| perm[(k * fold_sz)...((k + 1) * fold_sz)] }
    end

    def split_folds_contiguous(n_obs, k_folds)
      fold_sz = n_obs / k_folds
      Array.new(k_folds) { |k| ((k * fold_sz)...((k + 1) * fold_sz)).to_a }
    end

    def purged_train_indices(n_obs, test_start, test_end, horizon:, embargo:)
      # Observasi training t memiliki info di [t, t+horizon-1] yang overlap dengan tes jika
      # t >= test_start - horizon dan t < test_end + embargo.
      drop_lo = test_start - horizon + 1
      drop_hi = test_end + horizon + embargo - 1
      (0...n_obs).grep_v(drop_lo..drop_hi)
    end

    def eval_cv_oos_correlation(x_mat, y_labels, train_test_pairs, ridge_lambda: 1.0)
      preds = []
      actuals = []
      n_feat = x_mat.first.length
      train_test_pairs.each do |tr_idx, te_idx|
        w = fit_ridge(x_mat, y_labels, tr_idx, n_feat, ridge_lambda)
        te_idx.each do |t|
          preds << dot(x_mat[t], w)
          actuals << y_labels[t]
        end
      end
      Special.pearson_corr(preds, actuals)
    end

    def fit_ridge(x_mat, y_labels, tr_idx, n_feat, ridge_lambda)
      xtx = Array.new(n_feat) { Array.new(n_feat, 0.0) }
      xty = Array.new(n_feat, 0.0)
      tr_idx.each do |t|
        row = x_mat[t]
        yv = y_labels[t]
        n_feat.times do |i|
          xi = row[i]
          xty[i] += xi * yv
          (i...n_feat).each { |j| xtx[i][j] += xi * row[j] }
        end
      end
      symmetrize_and_reg!(xtx, n_feat, ridge_lambda)
      solve_linear(xtx, xty)
    end
    private_class_method :fit_ridge

    def symmetrize_and_reg!(xtx, n_feat, ridge_lambda)
      n_feat.times do |i|
        xtx[i][i] += ridge_lambda
        (0...i).each { |j| xtx[i][j] = xtx[j][i] }
      end
    end
    private_class_method :symmetrize_and_reg!

    def solve_linear(mat, rhs)
      n = rhs.length
      a = mat.map(&:dup)
      b = rhs.dup
      n.times do |col|
        piv = (col...n).max_by { |r| a[r][col].abs }
        a[col], a[piv] = a[piv], a[col]
        b[col], b[piv] = b[piv], b[col]
        eliminate_below!(a, b, col, n)
      end
      back_substitute(a, b, n)
    end
    private_class_method :solve_linear

    def eliminate_below!(a, b, col, n)
      diag = a[col][col]
      ((col + 1)...n).each do |r|
        fac = a[r][col] / diag
        (col...n).each { |c| a[r][c] -= fac * a[col][c] }
        b[r] -= fac * b[col]
      end
    end
    private_class_method :eliminate_below!

    def back_substitute(a, b, n)
      x = Array.new(n, 0.0)
      (n - 1).downto(0) do |i|
        sum = b[i]
        ((i + 1)...n).each { |j| sum -= a[i][j] * x[j] }
        x[i] = sum / a[i][i]
      end
      x
    end
    private_class_method :back_substitute

    def dot(xs, ys)
      s = 0.0
      xs.length.times { |i| s += xs[i] * ys[i] }
      s
    end
    private_class_method :dot

    def run_study(seed: MAIN_SEED, workers: nil)
      n_reps = 250
      n_obs = 500
      horizon = 10
      embargo = 10
      n_feat = 16
      k_folds = 5

      rows = Parallel.map_chunks(n_items: n_reps, workers: workers, chunk_size: 25) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 950, cidx))
        Array.new(len) { single_purged_rep(n_obs, horizon, embargo, n_feat, k_folds, rng) }
      end.flatten(1)

      summarize_purged_mc(rows, n_reps, n_obs, horizon, embargo, k_folds)
    end

    def single_purged_rep(n_obs, horizon, embargo, n_feat, k_folds, rng)
      x_mat, y_labels = generate_overlapping_dataset(n_obs: n_obs, horizon: horizon, n_features: n_feat, rng: rng)
      c_naive = run_naive_cv(x_mat, y_labels, n_obs, k_folds, rng)
      c_contig = run_contiguous_unpurged_cv(x_mat, y_labels, n_obs, k_folds)
      c_purged = run_purged_embargo_cv(x_mat, y_labels, n_obs, horizon, embargo, k_folds)
      [c_naive, c_contig, c_purged]
    end
    private_class_method :single_purged_rep

    def run_naive_cv(x_mat, y_labels, n_obs, k_folds, rng)
      folds = split_folds_naive(n_obs, k_folds, rng)
      all_idx = (0...n_obs).to_a
      pairs = folds.map { |te| [all_idx - te, te] }
      eval_cv_oos_correlation(x_mat, y_labels, pairs)
    end
    private_class_method :run_naive_cv

    def run_contiguous_unpurged_cv(x_mat, y_labels, n_obs, k_folds)
      folds = split_folds_contiguous(n_obs, k_folds)
      all_idx = (0...n_obs).to_a
      pairs = folds.map { |te| [all_idx - te, te] }
      eval_cv_oos_correlation(x_mat, y_labels, pairs)
    end
    private_class_method :run_contiguous_unpurged_cv

    def run_purged_embargo_cv(x_mat, y_labels, n_obs, horizon, embargo, k_folds)
      fold_sz = n_obs / k_folds
      pairs = Array.new(k_folds) do |k|
        ts = k * fold_sz
        te = ts + fold_sz - 1
        [purged_train_indices(n_obs, ts, te, horizon: horizon, embargo: embargo), (ts..te).to_a]
      end
      eval_cv_oos_correlation(x_mat, y_labels, pairs)
    end
    private_class_method :run_purged_embargo_cv

    def summarize_purged_mc(rows, n_reps, n_obs, horizon, embargo, k_folds)
      naives = rows.map { |r| r[0] }
      nm, nse = mean_and_se(naives, n_reps)
      cm, cse = mean_and_se(rows.map { |r| r[1] }, n_reps)
      purgeds = rows.map { |r| r[2] }
      pm, pse = mean_and_se(purgeds, n_reps)
      { n_reps: n_reps, n_obs: n_obs, horizon: horizon, embargo: embargo, k_folds: k_folds,
        naive_shuffled_corr_mean: nm, naive_shuffled_corr_se: nse,
        contiguous_unpurged_corr_mean: cm, contiguous_unpurged_corr_se: cse,
        purged_embargo_corr_mean: pm, purged_embargo_corr_se: pse,
        purged_zero_z: pm.abs / [pse, 1.0e-12].max,
        naive_sample: naives.first(150), purged_sample: purgeds.first(150) }
    end
    private_class_method :summarize_purged_mc

    def mean_and_se(arr, n_reps)
      m = arr.sum(0.0) / n_reps
      sd = Math.sqrt(arr.sum(0.0) { |x| (x - m)**2 } / (n_reps - 1))
      [m, sd / Math.sqrt(n_reps)]
    end
    private_class_method :mean_and_se
  end
end
