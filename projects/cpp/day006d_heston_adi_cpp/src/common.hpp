#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <new>
#include <vector>

struct Params {
    double S0 = 100.0, K = 100.0, T = 1.0, r = 0.03;
    double kappa = 1.5, theta = 0.04, sigma = 0.5, rho = -0.7, v0 = 0.04;
    static Params set_a() { return {}; }
    static Params set_b() { Params p; p.kappa = 3.0; p.sigma = 0.3; return p; }
};

struct Grid {
    int nx, nv, m, ld, i0, j0;
    double dx, dv, dt, x_min;
    Grid(int nx_, int nv_, int m_, const Params& p)
        : nx(nx_), nv(nv_), m(m_), ld((nx_ + 15) & ~15), i0((nx_ - 1) / 2),
          dx(6.0 / (nx_ - 1)), dv(1.0 / (nv_ - 1)), dt(p.T / m_), x_min(std::log(p.K) - 3.0) {
        j0 = static_cast<int>(std::llround(p.v0 / dv));
    }
};

struct AlignedBuf {
    double* ptr = nullptr;
    AlignedBuf() = default;
    explicit AlignedBuf(size_t n) {
        size_t sz = (n * sizeof(double) + 127) & ~size_t(127);
        ptr = static_cast<double*>(::operator new[](sz, std::align_val_t(128)));
        std::memset(ptr, 0, sz);
    }
    ~AlignedBuf() { if (ptr) ::operator delete[](ptr, std::align_val_t(128)); }
    AlignedBuf(const AlignedBuf&) = delete;
    AlignedBuf& operator=(const AlignedBuf&) = delete;
    double* get() { return ptr; }
    const double* get() const { return ptr; }
};

struct Problem {
    Params p; Grid g;
    std::vector<double> ax_l, ax_d, ax_u, av_l, av_d, av_u, m0, xl, xd, xu, vl, vd, vu, vinv, vcp;
    double ax_b0_d, ax_b0_u, ax_bn_l, ax_bn_d, av_e0, xb0_d, xb0_u, xbn_l, xbn_d, f0;
    AlignedBuf xinv_buf, xcp_buf;
    Problem(const Params& p_, int nx, int nv, int m) : p(p_), g(nx, nv, m, p_),
        ax_l(nv), ax_d(nv), ax_u(nv), av_l(nv), av_d(nv), av_u(nv), m0(nv),
        xl(nv), xd(nv), xu(nv), vl(nv), vd(nv), vu(nv), vinv(nv), vcp(nv),
        xinv_buf(size_t(nv) * g.ld), xcp_buf(size_t(nv) * g.ld) {
        double dx = g.dx, dv = g.dv, dt = g.dt;
        ax_b0_d = dt * (-p.r / dx - 0.5 * p.r); ax_b0_u = dt * (p.r / dx);
        ax_bn_l = dt * (-p.r / dx);             ax_bn_d = dt * (p.r / dx - 0.5 * p.r);
        xb0_d = 1.0 - 0.5 * ax_b0_d; xb0_u = -0.5 * ax_b0_u;
        xbn_l = -0.5 * ax_bn_l;      xbn_d = 1.0 - 0.5 * ax_bn_d;
        double *xinv = xinv_buf.get(), *xcp = xcp_buf.get();
        for (int j = 0; j < nv; ++j) {
            double v = j * dv, al = 0.5 * v / (dx * dx), be = (p.r - 0.5 * v) / (2.0 * dx);
            ax_l[j] = dt * (al - be); ax_d[j] = dt * (-2.0 * al - 0.5 * p.r); ax_u[j] = dt * (al + be);
            xl[j] = -0.5 * ax_l[j];   xd[j] = 1.0 - 0.5 * ax_d[j];           xu[j] = -0.5 * ax_u[j];
            xinv[j * g.ld] = 1.0 / xb0_d; xcp[j * g.ld] = xb0_u * xinv[j * g.ld];
            for (int i = 1; i < nx - 1; ++i) {
                double inv = 1.0 / (xd[j] - xl[j] * xcp[j * g.ld + i - 1]);
                xinv[j * g.ld + i] = inv; xcp[j * g.ld + i] = xu[j] * inv;
            }
            xinv[j * g.ld + nx - 1] = 1.0 / (xbn_d - xbn_l * xcp[j * g.ld + nx - 2]);
            m0[j] = (j == 0 || j == nv - 1) ? 0.0 : dt * p.rho * p.sigma * v / (4.0 * dx * dv);
            double gam = 0.5 * p.sigma * p.sigma * v / (dv * dv), eta = p.kappa * (p.theta - v) / (2.0 * dv);
            if (j == 0) {
                av_d[0] = dt * (-1.5 * p.kappa * p.theta / dv - 0.5 * p.r);
                av_u[0] = dt * (2.0 * p.kappa * p.theta / dv); av_e0 = dt * (-0.5 * p.kappa * p.theta / dv);
            } else if (j == nv - 1) {
                av_l[j] = dt * (2.0 * gam); av_d[j] = dt * (-2.0 * gam - 0.5 * p.r);
            } else {
                av_l[j] = dt * (gam - eta); av_d[j] = dt * (-2.0 * gam - 0.5 * p.r); av_u[j] = dt * (gam + eta);
            }
            vl[j] = -0.5 * av_l[j]; vd[j] = 1.0 - 0.5 * av_d[j]; vu[j] = -0.5 * av_u[j];
        }
        f0 = (-0.5 * av_e0) / vu[1]; vd[0] -= f0 * vl[1]; vu[0] -= f0 * vd[1];
        vinv[0] = 1.0 / vd[0]; vcp[0] = vu[0] * vinv[0];
        for (int j = 1; j < nv - 1; ++j) { vinv[j] = 1.0 / (vd[j] - vl[j] * vcp[j - 1]); vcp[j] = vu[j] * vinv[j]; }
        vinv[nv - 1] = 1.0 / (vd[nv - 1] - vl[nv - 1] * vcp[nv - 2]);
    }
};

inline double tri(double c0, double c1, double c2, double v0, double v1, double v2) { return (c0 * v0 + c1 * v1) + c2 * v2; }
inline double bi(double c0, double c1, double v0, double v1) { return c0 * v0 + c1 * v1; }
inline double mixed(double m, double ne, double nw, double se, double sw) { return m * ((ne - nw) - (se - sw)); }

inline void init_payoff_strided(double* u, const Grid& g, const Params& p) {
    for (int j = 0; j < g.nv; ++j) {
        double* row = u + size_t(j) * g.ld;
        for (int i = 0; i < g.nx; ++i)
            row[i] = (i == g.i0) ? (p.K * std::exp(0.5 * g.dx) - p.K - 0.5 * p.K * g.dx) / g.dx
                                 : std::max(std::exp(g.x_min + i * g.dx) - p.K, 0.0);
    }
}
inline void pack_dense(const double* src, double* dst, const Grid& g) {
    for (int j = 0; j < g.nv; ++j) std::memcpy(dst + size_t(j) * g.nx, src + size_t(j) * g.ld, size_t(g.nx) * sizeof(double));
}
inline uint64_t fnv1a_hash(const double* data, size_t n) {
    uint64_t h = 14695981039346656037ULL;
    const auto* b = reinterpret_cast<const uint8_t*>(data);
    for (size_t i = 0; i < n * sizeof(double); ++i) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}
inline double price_at(const double* dense, const Grid& g) { return dense[size_t(g.j0) * g.nx + g.i0]; }
