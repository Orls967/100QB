#include "variants.hpp"

namespace {
inline void x_solve_div(int nx, double xl, double xd, double xu,
                        double b0d, double b0u, double bnl, double bnd,
                        double* __restrict row, double* __restrict cp) {
    cp[0] = b0u / b0d;
    row[0] = row[0] / b0d;
    for (int i = 1; i < nx - 1; ++i) {
        double den = xd - xl * cp[i - 1];
        cp[i] = xu / den;
        row[i] = (row[i] - xl * row[i - 1]) / den;
    }
    double den = bnd - bnl * cp[nx - 2];
    row[nx - 1] = (row[nx - 1] - bnl * row[nx - 2]) / den;
    for (int i = nx - 2; i >= 0; --i) row[i] -= cp[i] * row[i + 1];
}

inline void v_solve_div(int nv, int ld, int i, const double* __restrict vl,
                        const double* __restrict vd, const double* __restrict vu,
                        double f0, double* __restrict grid, double* __restrict cp) {
    grid[i] -= f0 * grid[ld + i];
    cp[0] = vu[0] / vd[0];
    grid[i] /= vd[0];
    for (int j = 1; j < nv - 1; ++j) {
        double den = vd[j] - vl[j] * cp[j - 1];
        cp[j] = vu[j] / den;
        grid[size_t(j) * ld + i] = (grid[size_t(j) * ld + i] - vl[j] * grid[size_t(j - 1) * ld + i]) / den;
    }
    double den = vd[nv - 1] - vl[nv - 1] * cp[nv - 2];
    grid[size_t(nv - 1) * ld + i] = (grid[size_t(nv - 1) * ld + i] - vl[nv - 1] * grid[size_t(nv - 2) * ld + i]) / den;
    for (int j = nv - 2; j >= 0; --j) grid[size_t(j) * ld + i] -= cp[j] * grid[size_t(j + 1) * ld + i];
}
} // namespace

void solve_v1(const Problem& prob, std::vector<double>& out_dense) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld, m = prob.g.m;
    const size_t tot = size_t(nv) * ld;
    AlignedBuf u_b(tot), y_b(tot), w_b(tot), a1_b(tot), a2_b(tot), cp_b(std::max(nx, nv));
    double* __restrict U = u_b.get();  double* __restrict Y = y_b.get();
    double* __restrict W = w_b.get();  double* __restrict A1 = a1_b.get();
    double* __restrict A2 = a2_b.get(); double* __restrict cp = cp_b.get();
    init_payoff_strided(U, prob.g, prob.p);

    for (int step = 0; step < m; ++step) {
        for (int j = 0; j < nv; ++j) {
            const double* __restrict u = U + size_t(j) * ld;
            double* __restrict y = Y + size_t(j) * ld;
            double* __restrict w = W + size_t(j) * ld;
            double* __restrict a1r = A1 + size_t(j) * ld;
            double* __restrict a2r = A2 + size_t(j) * ld;
            const double* __restrict um = (j > 0) ? u - ld : u;
            const double* __restrict up = (j + 1 < nv) ? u + ld : u;
            for (int i = 0; i < nx; ++i) {
                double a0 = (i > 0 && i + 1 < nx && j > 0 && j + 1 < nv)
                    ? mixed(prob.m0[j], up[i+1], up[i-1], um[i+1], um[i-1]) : 0.0;
                double a1 = (i == 0) ? bi(prob.ax_b0_d, prob.ax_b0_u, u[0], u[1])
                          : (i == nx - 1) ? bi(prob.ax_bn_l, prob.ax_bn_d, u[nx-2], u[nx-1])
                          : tri(prob.ax_l[j], prob.ax_d[j], prob.ax_u[j], u[i-1], u[i], u[i+1]);
                double a2 = (j == 0) ? tri(prob.av_d[0], prob.av_u[0], prob.av_e0, u[i], u[ld+i], u[2*ld+i])
                          : (j == nv - 1) ? bi(prob.av_l[j], prob.av_d[j], um[i], u[i])
                          : tri(prob.av_l[j], prob.av_d[j], prob.av_u[j], um[i], u[i], up[i]);
                a1r[i] = a1; a2r[i] = a2;
                double y0 = ((u[i] + a0) + a1) + a2;
                w[i] = y0 - 0.5 * a0;
                y[i] = y0 - 0.5 * a1;
            }
            x_solve_div(nx, prob.xl[j], prob.xd[j], prob.xu[j], prob.xb0_d, prob.xb0_u, prob.xbn_l, prob.xbn_d, y, cp);
            for (int i = 0; i < nx; ++i) y[i] -= 0.5 * a2r[i];
        }
        for (int i = 0; i < nx; ++i)
            v_solve_div(nv, ld, i, prob.vl.data(), prob.vd.data(), prob.vu.data(), prob.f0, Y, cp);
        for (int j = 0; j < nv; ++j) {
            const double* __restrict y = Y + size_t(j) * ld;
            const double* __restrict ym = (j > 0) ? y - ld : y;
            const double* __restrict yp = (j + 1 < nv) ? y + ld : y;
            const double* __restrict w = W + size_t(j) * ld;
            const double* __restrict a1r = A1 + size_t(j) * ld;
            const double* __restrict a2r = A2 + size_t(j) * ld;
            double* __restrict u = U + size_t(j) * ld;
            for (int i = 0; i < nx; ++i) {
                double b0 = (i > 0 && i + 1 < nx && j > 0 && j + 1 < nv)
                    ? mixed(prob.m0[j], yp[i+1], yp[i-1], ym[i+1], ym[i-1]) : 0.0;
                u[i] = (w[i] + 0.5 * b0) - 0.5 * a1r[i];
            }
            x_solve_div(nx, prob.xl[j], prob.xd[j], prob.xu[j], prob.xb0_d, prob.xb0_u, prob.xbn_l, prob.xbn_d, u, cp);
            for (int i = 0; i < nx; ++i) u[i] -= 0.5 * a2r[i];
        }
        for (int i = 0; i < nx; ++i)
            v_solve_div(nv, ld, i, prob.vl.data(), prob.vd.data(), prob.vu.data(), prob.f0, U, cp);
    }
    out_dense.resize(size_t(nv) * nx);
    pack_dense(U, out_dense.data(), prob.g);
}
