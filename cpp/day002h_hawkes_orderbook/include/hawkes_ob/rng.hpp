#pragma once

#include <bit>
#include <cmath>
#include <cstdint>

namespace hob {

class SplitMix64 {
public:
  explicit constexpr SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}

  constexpr std::uint64_t next() noexcept {
    state_ += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state_;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

private:
  std::uint64_t state_;
};

[[nodiscard]] constexpr std::uint64_t derive_seed(std::uint64_t seed, std::size_t idx) noexcept {
  return seed ^ (static_cast<std::uint64_t>(idx + 1) * 0x9e3779b97f4a7c15ULL);
}

class Xoshiro256StarStar {
public:
  explicit constexpr Xoshiro256StarStar(std::uint64_t seed) noexcept {
    SplitMix64 sm(seed);
    s_[0] = sm.next();
    s_[1] = sm.next();
    s_[2] = sm.next();
    s_[3] = sm.next();
    if ((s_[0] | s_[1] | s_[2] | s_[3]) == 0ULL) {
      s_[0] = 1ULL;
    }
  }

  [[nodiscard]] constexpr std::uint64_t next_u64() noexcept {
    const std::uint64_t result = std::rotl(s_[1] * 5ULL, 7) * 9ULL;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = std::rotl(s_[3], 45);
    return result;
  }

  [[nodiscard]] constexpr double next_uniform() noexcept {
    constexpr double inv_2_53 = 1.0 / 9007199254740992.0;
    return static_cast<double>(next_u64() >> 11) * inv_2_53;
  }

  [[nodiscard]] double next_exp(double rate) noexcept {
    const double u = 1.0 - next_uniform();
    return -std::log(u > 1e-16 ? u : 1e-16) / rate;
  }

  /// Geometric distribution on {1, 2, ...} with mean `1 / p`.
  [[nodiscard]] std::int64_t next_geometric(double p) noexcept {
    const double u = 1.0 - next_uniform();
    const double val = std::floor(std::log(u > 1e-16 ? u : 1e-16) / std::log(1.0 - p)) + 1.0;
    return static_cast<std::int64_t>(val >= 1.0 ? val : 1.0);
  }

private:
  std::uint64_t s_[4]{};
};

} // namespace hob
