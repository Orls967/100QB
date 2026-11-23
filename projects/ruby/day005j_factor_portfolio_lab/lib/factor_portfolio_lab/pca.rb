# frozen_string_literal: true

module FactorPortfolioLab
  # Principal Component Analysis and Baik-Ben Arous-Peche (BBP) spiked covariance utilities.
  module Pca
    module_function

    def fit(rows, center: true, ddof: 1)
      t_len = rows.length
      raise ArgumentError, 'empty rows' if t_len.zero?

      n_assets = rows[0].length
      means = center ? column_means(rows, n_assets) : Array.new(n_assets, 0.0)
      centered = center ? rows.map { |r| r.each_with_index.map { |v, j| v - means[j] } } : rows
      cov = LinAlg.cov_matrix(centered, ddof: ddof)
      eig = LinAlg.jacobi_sym(cov)
      {
        eigenvalues: eig[:eigenvalues],
        eigenvectors: eig[:eigenvectors],
        means: means,
        cov: cov,
        n_obs: t_len,
        n_assets: n_assets
      }
    end

    def column_means(rows, n_assets)
      sums = Array.new(n_assets, 0.0)
      rows.each do |row|
        n_assets.times { |j| sums[j] += row[j] }
      end
      inv_t = 1.0 / rows.length
      sums.map { |s| s * inv_t }
    end

    def marchenko_pastur_bounds(gamma, sigma2: 1.0)
      root_g = Math.sqrt(gamma)
      {
        lambda_minus: sigma2 * ((1.0 - root_g)**2),
        lambda_plus: sigma2 * ((1.0 + root_g)**2),
        bbp_threshold: sigma2 * root_g
      }
    end

    def bbp_prediction(spike, gamma, sigma2: 1.0)
      ell = spike / sigma2
      thresh = Math.sqrt(gamma)
      if ell <= thresh
        {
          above_threshold: false,
          spike_excess: spike,
          sample_eigenvalue: sigma2 * ((1.0 + thresh)**2),
          alignment_sq: 0.0
        }
      else
        sample_eig = sigma2 * (1.0 + ell) * (1.0 + (gamma / ell))
        num = (ell * ell) - gamma
        den = (ell * ell) + (gamma * ell)
        {
          above_threshold: true,
          spike_excess: spike,
          sample_eigenvalue: sample_eig,
          alignment_sq: num / den
        }
      end
    end

    def eigenvector_alignment_sq(v1, v2)
      n = v1.length
      dot = 0.0
      n1 = 0.0
      n2 = 0.0
      n.times do |i|
        dot += v1[i] * v2[i]
        n1 += v1[i] * v1[i]
        n2 += v2[i] * v2[i]
      end
      denom = n1 * n2
      denom <= 0.0 ? 0.0 : (dot * dot) / denom
    end

    def effective_rank(eigenvalues)
      pos = eigenvalues.select { |e| e > 0.0 }
      return 0.0 if pos.empty?

      total = pos.sum
      probs = pos.map { |e| e / total }
      entropy = -probs.sum { |p| p * Math.log(p) }
      Math.exp(entropy)
    end
  end
end
