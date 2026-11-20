#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace hw::solvers {

template <typename Fn>
double brent_root(Fn&& f, double a, double b, double tol = 1e-13, int max_iter = 200) {
  double fa = f(a);
  double fb = f(b);
  if (fa * fb > 0.0) {
    throw std::runtime_error("brent_root: root not bracketed");
  }
  if (std::abs(fa) < std::abs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }
  double c = a;
  double fc = fa;
  bool mflag = true;
  double d = 0.0;

  for (int it = 0; it < max_iter; ++it) {
    if (std::abs(fb) < 1e-15 || std::abs(b - a) < tol) {
      return b;
    }
    double s = 0.0;
    if (fa != fc && fb != fc) {
      s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      s = b - fb * (b - a) / (fb - fa);
    }
    const double cond1 = (s < (3.0 * a + b) * 0.25 || s > b) && (b > a);
    const double cond1b = (s > (3.0 * a + b) * 0.25 || s < b) && (b < a);
    const bool out_range = (b > a) ? cond1 : cond1b;
    const bool cond2 = mflag && (std::abs(s - b) >= 0.5 * std::abs(b - c));
    const bool cond3 = !mflag && (std::abs(s - b) >= 0.5 * std::abs(c - d));
    const bool cond4 = mflag && (std::abs(b - c) < tol);
    const bool cond5 = !mflag && (std::abs(c - d) < tol);
    if (out_range || cond2 || cond3 || cond4 || cond5) {
      s = 0.5 * (a + b);
      mflag = true;
    } else {
      mflag = false;
    }
    const double fs = f(s);
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
  }
  return b;
}

inline std::vector<double> solve_linear(std::vector<std::vector<double>> a, std::vector<double> b) {
  const std::size_t n = b.size();
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t piv = col;
    double best = std::abs(a[col][col]);
    for (std::size_t row = col + 1; row < n; ++row) {
      if (std::abs(a[row][col]) > best) {
        best = std::abs(a[row][col]);
        piv = row;
      }
    }
    if (piv != col) {
      std::swap(a[col], a[piv]);
      std::swap(b[col], b[piv]);
    }
    if (std::abs(a[col][col]) < 1e-14) {
      throw std::runtime_error("Singular matrix in solve_linear");
    }
    for (std::size_t row = col + 1; row < n; ++row) {
      const double f = a[row][col] / a[col][col];
      for (std::size_t k = col; k < n; ++k) {
        a[row][k] -= f * a[col][k];
      }
      b[row] -= f * b[col];
    }
  }
  std::vector<double> x(n, 0.0);
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t k = i + 1; k < n; ++k) {
      s -= a[i][k] * x[k];
    }
    x[i] = s / a[i][i];
  }
  return x;
}

template <typename Fn>
std::vector<double> nelder_mead(Fn&& f, const std::vector<double>& x0,
                                const std::vector<double>& step, double tol = 1e-11,
                                int max_iter = 600) {
  const std::size_t n = x0.size();
  struct Vertex {
    std::vector<double> x;
    double val{0.0};
  };
  std::vector<Vertex> simp(n + 1);
  simp[0] = {x0, f(x0)};
  for (std::size_t i = 0; i < n; ++i) {
    auto xi = x0;
    xi[i] += step[i];
    simp[i + 1] = {xi, f(xi)};
  }

  for (int it = 0; it < max_iter; ++it) {
    std::sort(simp.begin(), simp.end(), [](const Vertex& u, const Vertex& v) { return u.val < v.val; });
    double max_span = 0.0;
    for (std::size_t i = 1; i <= n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        max_span = std::max(max_span, std::abs(simp[i].x[j] - simp[0].x[j]));
      }
    }
    if (max_span < tol && std::abs(simp[n].val - simp[0].val) < tol * tol) {
      break;
    }

    std::vector<double> cent(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        cent[j] += simp[i].x[j] / static_cast<double>(n);
      }
    }
    std::vector<double> xr(n);
    for (std::size_t j = 0; j < n; ++j) {
      xr[j] = cent[j] + (cent[j] - simp[n].x[j]);
    }
    const double fr = f(xr);
    if (fr >= simp[0].val && fr < simp[n - 1].val) {
      simp[n] = {xr, fr};
      continue;
    }
    if (fr < simp[0].val) {
      std::vector<double> xe(n);
      for (std::size_t j = 0; j < n; ++j) {
        xe[j] = cent[j] + 2.0 * (xr[j] - cent[j]);
      }
      const double fe = f(xe);
      simp[n] = (fe < fr) ? Vertex{xe, fe} : Vertex{xr, fr};
      continue;
    }
    std::vector<double> xc(n);
    for (std::size_t j = 0; j < n; ++j) {
      xc[j] = cent[j] + 0.5 * (simp[n].x[j] - cent[j]);
    }
    const double fc = f(xc);
    if (fc < simp[n].val) {
      simp[n] = {xc, fc};
      continue;
    }
    for (std::size_t i = 1; i <= n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        simp[i].x[j] = simp[0].x[j] + 0.5 * (simp[i].x[j] - simp[0].x[j]);
      }
      simp[i].val = f(simp[i].x);
    }
  }
  std::sort(simp.begin(), simp.end(), [](const Vertex& u, const Vertex& v) { return u.val < v.val; });
  return simp[0].x;
}

} // namespace hw::solvers
