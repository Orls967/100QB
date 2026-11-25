#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace vicsek {

inline constexpr uint64_t MASTER_SEED = 20261126ULL;
inline constexpr double PI = 3.14159265358979323846;
inline constexpr double TWO_PI = 2.0 * PI;

inline uint64_t splitmix64_mix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

struct Xoshiro256pp {
    uint64_t s[4]{};
    explicit Xoshiro256pp(uint64_t seed = MASTER_SEED) {
        uint64_t st = seed;
        for (int i = 0; i < 4; ++i) {
            st += 0x9E3779B97F4A7C15ULL;
            s[i] = splitmix64_mix(st);
        }
    }
    inline uint64_t next_u64() {
        auto rotl = [](uint64_t x, int k) { return (x << k) | (x >> (64 - k)); };
        uint64_t res = rotl(s[0] + s[3], 23) + s[0];
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t; s[3] = rotl(s[3], 45);
        return res;
    }
    inline double u01() { return static_cast<double>(next_u64() >> 11) * 0x1.0p-53; }
};

inline double wrap_pos(double x, double L) {
    x = std::fmod(x, L);
    return (x < 0.0) ? (x + L) : x;
}

inline double min_img(double d, double L) {
    double h = 0.5 * L;
    if (d > h) d -= L;
    else if (d < -h) d += L;
    return d;
}

struct Sim {
    int N;
    double L, v0 = 0.03, r = 1.0, eta;
    std::vector<double> x, y, theta, cos_t, sin_t, new_t;
    int M;
    double cell_w;
    std::vector<int> head, next, nbr_buf;
    Xoshiro256pp rng;

    Sim(int n, double l, double noise, uint64_t seed, bool uniform_theta = false)
        : N(n), L(l), eta(noise), x(n), y(n), theta(n), cos_t(n), sin_t(n), new_t(n),
          M(std::max(1, static_cast<int>(std::floor(l / 1.0)))), cell_w(l / M),
          head(M * M, -1), next(n, -1), rng(seed) {
        nbr_buf.reserve(64);
        for (int i = 0; i < N; ++i) {
            x[i] = rng.u01() * L;
            y[i] = rng.u01() * L;
            theta[i] = uniform_theta ? 0.0 : ((rng.u01() - 0.5) * TWO_PI);
        }
    }

    double polarization() const {
        double sc = 0.0, ss = 0.0;
        for (int i = 0; i < N; ++i) {
            sc += std::cos(theta[i]);
            ss += std::sin(theta[i]);
        }
        return std::sqrt(sc * sc + ss * ss) / static_cast<double>(N);
    }

    void precompute_trig() {
        for (int i = 0; i < N; ++i) {
            cos_t[i] = std::cos(theta[i]);
            sin_t[i] = std::sin(theta[i]);
        }
    }

    void step_bruteforce() {
        precompute_trig();
        const double r2 = r * r;
        for (int i = 0; i < N; ++i) {
            double sc = 0.0, ss = 0.0;
            for (int j = 0; j < N; ++j) {
                double dx = min_img(x[j] - x[i], L);
                double dy = min_img(y[j] - y[i], L);
                if (dx * dx + dy * dy < r2) {
                    sc += cos_t[j];
                    ss += sin_t[j];
                }
            }
            new_t[i] = std::atan2(ss, sc) + eta * (rng.u01() - 0.5);
        }
        advance();
    }

    void step_cell_list(bool sort_neighbors = true) {
        precompute_trig();
        std::fill(head.begin(), head.end(), -1);
        for (int i = 0; i < N; ++i) {
            int cx = std::min(M - 1, static_cast<int>(x[i] / cell_w));
            int cy = std::min(M - 1, static_cast<int>(y[i] / cell_w));
            int c = cy * M + cx;
            next[i] = head[c];
            head[c] = i;
        }
        const double r2 = r * r;
        const int span = (M < 3) ? 0 : 1;
        for (int i = 0; i < N; ++i) {
            int cx = std::min(M - 1, static_cast<int>(x[i] / cell_w));
            int cy = std::min(M - 1, static_cast<int>(y[i] / cell_w));
            double sc = 0.0, ss = 0.0;
            if (sort_neighbors) {
                nbr_buf.clear();
                if (M < 3) {
                    for (int c = 0; c < M * M; ++c)
                        for (int j = head[c]; j != -1; j = next[j]) nbr_buf.push_back(j);
                } else {
                    for (int dy = -span; dy <= span; ++dy) {
                        int ny = (cy + dy + M) % M;
                        for (int dx = -span; dx <= span; ++dx) {
                            int nx = (cx + dx + M) % M;
                            for (int j = head[ny * M + nx]; j != -1; j = next[j]) {
                                double ddx = min_img(x[j] - x[i], L);
                                double ddy = min_img(y[j] - y[i], L);
                                if (ddx * ddx + ddy * ddy < r2) nbr_buf.push_back(j);
                            }
                        }
                    }
                }
                std::sort(nbr_buf.begin(), nbr_buf.end());
                for (int j : nbr_buf) {
                    if (M < 3) {
                        double ddx = min_img(x[j] - x[i], L);
                        double ddy = min_img(y[j] - y[i], L);
                        if (ddx * ddx + ddy * ddy >= r2) continue;
                    }
                    sc += cos_t[j];
                    ss += sin_t[j];
                }
            } else {
                for (int dy = -span; dy <= span; ++dy) {
                    int ny = (cy + dy + M) % M;
                    for (int dx = -span; dx <= span; ++dx) {
                        int nx = (cx + dx + M) % M;
                        for (int j = head[ny * M + nx]; j != -1; j = next[j]) {
                            double ddx = min_img(x[j] - x[i], L);
                            double ddy = min_img(y[j] - y[i], L);
                            if (ddx * ddx + ddy * ddy < r2) {
                                sc += cos_t[j];
                                ss += sin_t[j];
                            }
                        }
                    }
                }
            }
            new_t[i] = std::atan2(ss, sc) + eta * (rng.u01() - 0.5);
        }
        advance();
    }

    void advance() {
        for (int i = 0; i < N; ++i) {
            theta[i] = new_t[i];
            x[i] = wrap_pos(x[i] + v0 * std::cos(theta[i]), L);
            y[i] = wrap_pos(y[i] + v0 * std::sin(theta[i]), L);
        }
    }
};

}  // namespace vicsek
