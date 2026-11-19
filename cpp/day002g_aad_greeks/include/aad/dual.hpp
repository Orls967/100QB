#pragma once

#include <cmath>

#include "aad/number.hpp"

namespace aad {

template <typename T>
struct Dual {
  T val;
  T dot;

  Dual() : val(0.0), dot(0.0) {}
  Dual(double v) : val(v), dot(0.0) {} // NOLINT(google-explicit-constructor)
  Dual(T v, T d) : val(std::move(v)), dot(std::move(d)) {}

  Dual operator-() const { return Dual(-val, -dot); }

  Dual& operator+=(const Dual& rhs) {
    val += rhs.val;
    dot += rhs.dot;
    return *this;
  }
  Dual& operator-=(const Dual& rhs) {
    val -= rhs.val;
    dot -= rhs.dot;
    return *this;
  }
  Dual& operator*=(const Dual& rhs) {
    *this = *this * rhs;
    return *this;
  }
  Dual& operator/=(const Dual& rhs) {
    *this = *this / rhs;
    return *this;
  }

  friend Dual operator+(const Dual& a, const Dual& b) { return Dual(a.val + b.val, a.dot + b.dot); }
  friend Dual operator-(const Dual& a, const Dual& b) { return Dual(a.val - b.val, a.dot - b.dot); }
  friend Dual operator*(const Dual& a, const Dual& b) {
    return Dual(a.val * b.val, a.dot * b.val + a.val * b.dot);
  }
  friend Dual operator/(const Dual& a, const Dual& b) {
    return Dual(a.val / b.val, (a.dot * b.val - a.val * b.dot) / (b.val * b.val));
  }
  friend bool operator<(const Dual& a, const Dual& b) noexcept { return a.val < b.val; }
  friend bool operator>(const Dual& a, const Dual& b) noexcept { return a.val > b.val; }
};

namespace math {

template <typename T>
inline Dual<T> exp(const Dual<T>& x) {
  T e = exp(x.val);
  return Dual<T>(e, e * x.dot);
}

template <typename T>
inline Dual<T> log(const Dual<T>& x) {
  return Dual<T>(log(x.val), x.dot / x.val);
}

template <typename T>
inline Dual<T> sqrt(const Dual<T>& x) {
  T s = sqrt(x.val);
  return Dual<T>(s, x.dot / (T(2.0) * s));
}

template <typename T>
inline Dual<T> pow(const Dual<T>& x, const Dual<T>& y) {
  T p = pow(x.val, y.val);
  T dp = y.val * pow(x.val, y.val - T(1.0)) * x.dot + p * log(x.val) * y.dot;
  return Dual<T>(p, dp);
}

template <typename T>
inline Dual<T> abs(const Dual<T>& x) {
  if (value_of(x.val) >= 0.0) {
    return Dual<T>(x.val, x.dot);
  }
  return Dual<T>(-x.val, -x.dot);
}

template <typename T>
inline Dual<T> max(const Dual<T>& a, const Dual<T>& b) {
  return (value_of(a.val) >= value_of(b.val)) ? a : b;
}

template <typename T>
inline Dual<T> normal_cdf(const Dual<T>& x) {
  constexpr double inv_sqrt_2pi = 0.39894228040143267794;
  T pdf = T(inv_sqrt_2pi) * exp(T(-0.5) * x.val * x.val);
  return Dual<T>(normal_cdf(x.val), pdf * x.dot);
}

} // namespace math

} // namespace aad
