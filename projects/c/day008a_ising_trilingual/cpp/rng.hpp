#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace ising {

inline constexpr uint64_t MASTER_SEED = 20261126ULL;

inline uint64_t splitmix64_mix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

inline uint64_t splitmix64_next(uint64_t &state) {
    state += 0x9E3779B97F4A7C15ULL;
    return splitmix64_mix(state);
}

struct Xoshiro256pp {
    uint64_t s[4]{};

    explicit Xoshiro256pp(uint64_t seed = MASTER_SEED) { reseed(seed); }

    void reseed(uint64_t seed) {
        uint64_t st = seed;
        s[0] = splitmix64_next(st);
        s[1] = splitmix64_next(st);
        s[2] = splitmix64_next(st);
        s[3] = splitmix64_next(st);
    }

    static inline uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    inline uint64_t next_u64() {
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    inline double u01() {
        return static_cast<double>(next_u64() >> 11) * 0x1.0p-53;
    }

    inline int rand_int(int n) {
        double u = u01();
        int idx = static_cast<int>(std::floor(u * static_cast<double>(n)));
        if (idx >= n) idx = n - 1;
        return idx;
    }
};

inline uint64_t compute_task_seed(int L, int t_index, int algo_id, int replica) {
    uint64_t h = 0x517CC1B727220A95ULL;
    h = splitmix64_mix(h ^ static_cast<uint64_t>(L));
    h = splitmix64_mix(h ^ static_cast<uint64_t>(t_index));
    h = splitmix64_mix(h ^ static_cast<uint64_t>(algo_id));
    h = splitmix64_mix(h ^ static_cast<uint64_t>(replica));
    return splitmix64_mix(MASTER_SEED ^ h);
}

inline uint64_t fnv1a_64(const int8_t *data, size_t len) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= static_cast<uint8_t>(data[i]);
        h *= 0x100000001B3ULL;
    }
    return h;
}

}  // namespace ising
