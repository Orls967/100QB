#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <thread>
#include <vector>

namespace wator {

inline constexpr uint64_t MASTER_SEED = 20261126ULL;

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
    inline int rand_int(int n) {
        int idx = static_cast<int>(std::floor(u01() * static_cast<double>(n)));
        return (idx >= n) ? (n - 1) : idx;
    }
};

enum Species : uint8_t { EMPTY = 0, FISH = 1, SHARK = 2 };

struct Cell {
    Species kind = EMPTY;
    int16_t age = 0;
    int16_t energy = 0;
};

struct Params {
    int W = 200, H = 200;
    int fish_breed = 3, shark_breed = 10;
    int shark_init_energy = 5, fish_energy = 3;
    double init_fish = 0.30, init_shark = 0.05;
};

struct World {
    Params p;
    std::vector<Cell> grid;
    Xoshiro256pp rng;
    int n_fish = 0, n_sharks = 0;
    std::vector<int> order;

    World(Params params, uint64_t seed) : p(params), grid(static_cast<size_t>(p.W * p.H)), rng(seed) {
        int N = p.W * p.H;
        std::vector<int> idx(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) idx[static_cast<size_t>(i)] = i;
        for (int i = N - 1; i > 0; --i) std::swap(idx[static_cast<size_t>(i)], idx[static_cast<size_t>(rng.rand_int(i + 1))]);
        int nf = static_cast<int>(std::round(p.init_fish * N));
        int ns = static_cast<int>(std::round(p.init_shark * N));
        for (int k = 0; k < nf; ++k) grid[static_cast<size_t>(idx[static_cast<size_t>(k)])] = {FISH, static_cast<int16_t>(rng.rand_int(p.fish_breed)), 0};
        for (int k = nf; k < nf + ns; ++k) grid[static_cast<size_t>(idx[static_cast<size_t>(k)])] = {SHARK, static_cast<int16_t>(rng.rand_int(p.shark_breed)), static_cast<int16_t>(p.shark_init_energy)};
        n_fish = nf; n_sharks = ns;
    }

    inline int wrap_nbr(int i, int dir) const {
        int x = i % p.W, y = i / p.W;
        if (dir == 0) x = (x + 1 == p.W) ? 0 : x + 1;
        else if (dir == 1) y = (y + 1 == p.H) ? 0 : y + 1;
        else if (dir == 2) x = (x == 0) ? p.W - 1 : x - 1;
        else y = (y == 0) ? p.H - 1 : y - 1;
        return y * p.W + x;
    }

    bool check_invariants() const {
        int cf = 0, cs = 0;
        for (const auto &c : grid) {
            if (c.kind == FISH) ++cf;
            else if (c.kind == SHARK) {
                if (c.energy < 1) return false;
                ++cs;
            }
        }
        return (cf == n_fish) && (cs == n_sharks);
    }

    void shuffle_species(Species sp) {
        order.clear();
        for (int i = 0; i < static_cast<int>(grid.size()); ++i)
            if (grid[static_cast<size_t>(i)].kind == sp) order.push_back(i);
        for (int i = static_cast<int>(order.size()) - 1; i > 0; --i)
            std::swap(order[static_cast<size_t>(i)], order[static_cast<size_t>(rng.rand_int(i + 1))]);
    }

    void step() {
        int nbrs[4];
        shuffle_species(FISH);
        for (int i : order) {
            if (grid[static_cast<size_t>(i)].kind != FISH) continue;
            grid[static_cast<size_t>(i)].age++;
            int k = 0;
            for (int d = 0; d < 4; ++d) {
                int nb = wrap_nbr(i, d);
                if (grid[static_cast<size_t>(nb)].kind == EMPTY) nbrs[k++] = nb;
            }
            if (k > 0) {
                int dst = nbrs[rng.rand_int(k)];
                Cell c = grid[static_cast<size_t>(i)];
                if (c.age >= p.fish_breed) {
                    grid[static_cast<size_t>(dst)] = {FISH, 0, 0};
                    grid[static_cast<size_t>(i)] = {FISH, 0, 0};
                    ++n_fish;
                } else {
                    grid[static_cast<size_t>(dst)] = c;
                    grid[static_cast<size_t>(i)] = {EMPTY, 0, 0};
                }
            }
        }
        shuffle_species(SHARK);
        for (int i : order) {
            if (grid[static_cast<size_t>(i)].kind != SHARK) continue;
            Cell c = grid[static_cast<size_t>(i)];
            c.age++; c.energy--;
            if (c.energy <= 0) {
                grid[static_cast<size_t>(i)] = {EMPTY, 0, 0};
                --n_sharks;
                continue;
            }
            int kf = 0, ke = 0, fnbr[4], enbr[4];
            for (int d = 0; d < 4; ++d) {
                int nb = wrap_nbr(i, d);
                if (grid[static_cast<size_t>(nb)].kind == FISH) fnbr[kf++] = nb;
                else if (grid[static_cast<size_t>(nb)].kind == EMPTY) enbr[ke++] = nb;
            }
            int dst = -1;
            if (kf > 0) {
                dst = fnbr[rng.rand_int(kf)];
                c.energy += static_cast<int16_t>(p.fish_energy);
                --n_fish;
            } else if (ke > 0) {
                dst = enbr[rng.rand_int(ke)];
            }
            if (dst >= 0) {
                if (c.age >= p.shark_breed) {
                    grid[static_cast<size_t>(dst)] = {SHARK, 0, c.energy};
                    grid[static_cast<size_t>(i)] = {SHARK, 0, static_cast<int16_t>(p.shark_init_energy)};
                    ++n_sharks;
                } else {
                    grid[static_cast<size_t>(dst)] = c;
                    grid[static_cast<size_t>(i)] = {EMPTY, 0, 0};
                }
            } else {
                grid[static_cast<size_t>(i)] = c;
            }
        }
    }
};

struct LVParams { double a, b, c, d; };

inline LVParams fit_lv_ols(const std::vector<int> &F, const std::vector<int> &S, int t0) {
    auto ols_line = [](const std::vector<double> &x, const std::vector<double> &y) {
        double mx = 0, my = 0;
        for (size_t i = 0; i < x.size(); ++i) { mx += x[i]; my += y[i]; }
        mx /= static_cast<double>(x.size()); my /= static_cast<double>(y.size());
        double num = 0, den = 0;
        for (size_t i = 0; i < x.size(); ++i) { num += (x[i] - mx) * (y[i] - my); den += (x[i] - mx) * (x[i] - mx); }
        double slope = (den > 1e-18) ? (num / den) : 0.0;
        return std::pair<double, double>{my - slope * mx, slope};
    };
    std::vector<double> xf, yf, xs, ys;
    for (size_t t = static_cast<size_t>(t0); t + 1 < F.size(); ++t) {
        if (F[t] > 0 && S[t] > 0) {
            xf.push_back(static_cast<double>(S[t]));
            yf.push_back(static_cast<double>(F[t + 1] - F[t]) / static_cast<double>(F[t]));
            xs.push_back(static_cast<double>(F[t]));
            ys.push_back(static_cast<double>(S[t + 1] - S[t]) / static_cast<double>(S[t]));
        }
    }
    auto [a, neg_b] = ols_line(xf, yf);
    auto [neg_c, d] = ols_line(xs, ys);
    return {a, -neg_b, -neg_c, d};
}

inline void rk4_lv_step(double &F, double &S, double dt, const LVParams &lv) {
    auto fF = [&](double f, double s) { return (lv.a - lv.b * s) * f; };
    auto fS = [&](double f, double s) { return (-lv.c + lv.d * f) * s; };
    double k1F = fF(F, S), k1S = fS(F, S);
    double k2F = fF(F + 0.5 * dt * k1F, S + 0.5 * dt * k1S), k2S = fS(F + 0.5 * dt * k1F, S + 0.5 * dt * k1S);
    double k3F = fF(F + 0.5 * dt * k2F, S + 0.5 * dt * k2S), k3S = fS(F + 0.5 * dt * k2F, S + 0.5 * dt * k2S);
    double k4F = fF(F + dt * k3F, S + dt * k3S), k4S = fS(F + dt * k3F, S + dt * k3S);
    F += (dt / 6.0) * (k1F + 2.0 * k2F + 2.0 * k3F + k4F);
    S += (dt / 6.0) * (k1S + 2.0 * k2S + 2.0 * k3S + k4S);
}

inline double lv_first_integral(double F, double S, const LVParams &lv) {
    return lv.d * F - lv.c * std::log(F) + lv.b * S - lv.a * std::log(S);
}

}  // namespace wator
