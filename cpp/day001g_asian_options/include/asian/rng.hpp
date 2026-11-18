#pragma once

#include <bit>
#include <cmath>
#include <cstdint>

namespace asian {

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

[[nodiscard]] constexpr std::uint64_t derive_chunk_seed(std::uint64_t seed,
                                                        std::size_t chunk_idx) noexcept {
  return seed ^ (static_cast<std::uint64_t>(chunk_idx + 1) * 0x9e3779b97f4a7c15ULL);
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

  [[nodiscard]] double next_normal() noexcept {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    double u = 0.0;
    double v = 0.0;
    double s = 0.0;
    do {
      u = 2.0 * next_uniform() - 1.0;
      v = 2.0 * next_uniform() - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    const double mul = std::sqrt(-2.0 * std::log(s) / s);
    spare_ = v * mul;
    has_spare_ = true;
    return u * mul;
  }

private:
  std::uint64_t s_[4]{};
  bool has_spare_{false};
  double spare_{0.0};
};

} // namespace asian
