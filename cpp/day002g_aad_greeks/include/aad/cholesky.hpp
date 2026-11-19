#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "aad/number.hpp"

namespace aad {

/// Builds d x d correlation matrix from `rho_ij` (i < j, row-major upper triangle)
/// and returns its lower-triangular Cholesky factor `L` (`L * L^T = C`).
template <typename Scalar>
std::vector<std::vector<Scalar>> cholesky_from_rhos(const std::vector<Scalar>& rhos,
                                                      std::size_t d) {
  std::vector<std::vector<Scalar>> corr(d, std::vector<Scalar>(d, Scalar(0.0)));
  std::size_t k = 0;
  for (std::size_t i = 0; i < d; ++i) {
    corr[i][i] = Scalar(1.0);
    for (std::size_t j = i + 1; j < d; ++j) {
      corr[i][j] = rhos[k];
      corr[j][i] = rhos[k];
      ++k;
    }
  }

  std::vector<std::vector<Scalar>> lower(d, std::vector<Scalar>(d, Scalar(0.0)));
  for (std::size_t i = 0; i < d; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      Scalar sum = corr[i][j];
      for (std::size_t m = 0; m < j; ++m) {
        sum = sum - lower[i][m] * lower[j][m];
      }
      if (i == j) {
        if (math::value_of(sum) <= 0.0) {
          throw std::runtime_error("Correlation matrix is not positive-definite");
        }
        lower[i][j] = math::sqrt(sum);
      } else {
        lower[i][j] = sum / lower[j][j];
      }
    }
  }
  return lower;
}

} // namespace aad
