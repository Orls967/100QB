#pragma once

#include <array>
#include <cstdint>

namespace lbm {

inline constexpr std::uint64_t MASTER_SEED = 20261127ULL;

struct SplitMix64 {
    std::uint64_t state;

    explicit constexpr SplitMix64(std::uint64_t seed) noexcept : state(seed) {}

    constexpr std::uint64_t next() noexcept {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

struct Xoshiro256pp {
    std::array<std::uint64_t, 4> s{};

    explicit constexpr Xoshiro256pp(std::uint64_t seed) noexcept {
        SplitMix64 sm(seed);
        for (int i = 0; i < 4; ++i) {
            s[static_cast<std::size_t>(i)] = sm.next();
        }
    }

    static constexpr std::uint64_t rotl(std::uint64_t x, int k) noexcept {
        return (x << k) | (x >> (64 - k));
    }

    constexpr std::uint64_t next_u64() noexcept {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    // Uniform double in [0, 1)
    constexpr double uniform01() noexcept {
        return static_cast<double>(next_u64() >> 11) * 0x1.0p-53;
    }

    // Uniform double in [-1, 1)
    constexpr double uniform_sym() noexcept {
        return 2.0 * uniform01() - 1.0;
    }
};

inline constexpr std::uint64_t derive_seed(int stage_id, int task_id = 0,
                                           std::uint64_t master = MASTER_SEED) noexcept {
    SplitMix64 sm(master ^ (static_cast<std::uint64_t>(stage_id + 1) * 0x9E3779B97F4A7C15ULL));
    const std::uint64_t h = sm.next();
    SplitMix64 sm2(h ^ (static_cast<std::uint64_t>(task_id + 1) * 0xBF58476D1CE4E5B9ULL));
    return sm2.next();
}

}  // namespace lbm
