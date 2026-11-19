#pragma once

#include <cmath>
#include <cstddef>

#include "aad/tape.hpp"

namespace aad {

class Number {
public:
  Number() : idx_(active_tape().record_leaf(0.0)) {}
  Number(double val) : idx_(active_tape().record_leaf(val)) {} // NOLINT(google-explicit-constructor)
  static Number from_index(std::size_t idx) noexcept {
    Number n(NoInit{});
    n.idx_ = idx;
    return n;
  }

  [[nodiscard]] double value() const noexcept { return active_tape().value(idx_); }
  [[nodiscard]] double adjoint() const noexcept { return active_tape().adjoint(idx_); }
  [[nodiscard]] std::size_t index() const noexcept { return idx_; }

  Number operator-() const {
    return from_index(active_tape().record_unary(-value(), idx_, -1.0));
  }

  Number& operator+=(const Number& rhs) {
    *this = *this + rhs;
    return *this;
  }
  Number& operator-=(const Number& rhs) {
    *this = *this - rhs;
    return *this;
  }
  Number& operator*=(const Number& rhs) {
    *this = *this * rhs;
    return *this;
  }
  Number& operator/=(const Number& rhs) {
    *this = *this / rhs;
    return *this;
  }

  friend Number operator+(const Number& lhs, const Number& rhs) {
    return from_index(
        active_tape().record_binary(lhs.value() + rhs.value(), lhs.idx_, 1.0, rhs.idx_, 1.0));
  }
  friend Number operator-(const Number& lhs, const Number& rhs) {
    return from_index(
        active_tape().record_binary(lhs.value() - rhs.value(), lhs.idx_, 1.0, rhs.idx_, -1.0));
  }
  friend Number operator*(const Number& lhs, const Number& rhs) {
    const double lv = lhs.value();
    const double rv = rhs.value();
    return from_index(active_tape().record_binary(lv * rv, lhs.idx_, rv, rhs.idx_, lv));
  }
  friend Number operator/(const Number& lhs, const Number& rhs) {
    const double lv = lhs.value();
    const double rv = rhs.value();
    const double inv = 1.0 / rv;
    return from_index(
        active_tape().record_binary(lv * inv, lhs.idx_, inv, rhs.idx_, -lv * inv * inv));
  }

  friend bool operator<(const Number& lhs, const Number& rhs) noexcept {
    return lhs.value() < rhs.value();
  }
  friend bool operator<=(const Number& lhs, const Number& rhs) noexcept {
    return lhs.value() <= rhs.value();
  }
  friend bool operator>(const Number& lhs, const Number& rhs) noexcept {
    return lhs.value() > rhs.value();
  }
  friend bool operator>=(const Number& lhs, const Number& rhs) noexcept {
    return lhs.value() >= rhs.value();
  }

private:
  struct NoInit {};
  explicit Number(NoInit) noexcept : idx_(0) {}
  std::size_t idx_;
};

namespace math {

inline double exp(double x) noexcept { return std::exp(x); }
inline double log(double x) noexcept { return std::log(x); }
inline double sqrt(double x) noexcept { return std::sqrt(x); }
inline double pow(double x, double y) noexcept { return std::pow(x, y); }
inline double abs(double x) noexcept { return std::abs(x); }
inline double max(double a, double b) noexcept { return a > b ? a : b; }
inline double min(double a, double b) noexcept { return a < b ? a : b; }
inline double normal_pdf(double x) noexcept {
  constexpr double inv_sqrt_2pi = 0.39894228040143267794;
  return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}
inline double normal_cdf(double x) noexcept {
  constexpr double inv_sqrt2 = 0.70710678118654752440;
  return 0.5 * std::erfc(-x * inv_sqrt2);
}

inline Number exp(const Number& x) {
  const double v = std::exp(x.value());
  return Number::from_index(active_tape().record_unary(v, x.index(), v));
}

inline Number log(const Number& x) {
  const double xv = x.value();
  return Number::from_index(active_tape().record_unary(std::log(xv), x.index(), 1.0 / xv));
}

inline Number sqrt(const Number& x) {
  const double v = std::sqrt(x.value());
  const double deriv = (v > 0.0) ? (0.5 / v) : 0.0;
  return Number::from_index(active_tape().record_unary(v, x.index(), deriv));
}

inline Number pow(const Number& x, const Number& y) {
  const double xv = x.value();
  const double yv = y.value();
  const double v = std::pow(xv, yv);
  const double dx = (xv != 0.0) ? (yv * std::pow(xv, yv - 1.0)) : 0.0;
  const double dy = (xv > 0.0) ? (v * std::log(xv)) : 0.0;
  return Number::from_index(active_tape().record_binary(v, x.index(), dx, y.index(), dy));
}

inline Number pow(const Number& x, double y) {
  const double xv = x.value();
  const double v = std::pow(xv, y);
  const double dx = (xv != 0.0) ? (y * std::pow(xv, y - 1.0)) : 0.0;
  return Number::from_index(active_tape().record_unary(v, x.index(), dx));
}

inline Number abs(const Number& x) {
  const double xv = x.value();
  const double deriv = (xv > 0.0) ? 1.0 : ((xv < 0.0) ? -1.0 : 0.0);
  return Number::from_index(active_tape().record_unary(std::abs(xv), x.index(), deriv));
}

inline Number max(const Number& a, const Number& b) {
  if (a.value() >= b.value()) {
    return Number::from_index(active_tape().record_unary(a.value(), a.index(), 1.0));
  }
  return Number::from_index(active_tape().record_unary(b.value(), b.index(), 1.0));
}

inline Number normal_cdf(const Number& x) {
  const double xv = x.value();
  return Number::from_index(active_tape().record_unary(normal_cdf(xv), x.index(), normal_pdf(xv)));
}

inline double value_of(double x) noexcept { return x; }
inline double value_of(const Number& x) noexcept { return x.value(); }

} // namespace math

} // namespace aad
