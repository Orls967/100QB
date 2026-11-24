#include "variants.hpp"

namespace {
void thomas_div(int n, const double* lo, const double* di, const double* up,
                const double* r, double* x) {
    std::vector<double> cp(n), y(n);
    cp[0] = up[0] / di[0];
    y[0] = r[0] / di[0];
    for (int k = 1; k < n; ++k) {
        double den = di[k] - lo[k] * cp[k - 1];
        cp[k] = (k + 1 < n) ? (up[k] / den) : 0.0;
        y[k] = (r[k] - lo[k] * y[k - 1]) / den;
    }
    x[n - 1] = y[n - 1];
    for (int k = n - 2; k >= 0; --k) x[k] = y[k] - cp[k] * x[k + 1];
}
} // namespace

void solve_v0(const Problem& prob, std::vector<double>& out_dense) {
    const int nx = prob.g.nx, nv = prob.g.nv, m = prob.g.m;
    std::vector<std::vector<double>> U(nx, std::vector<double>(nv, 0.0));
    for (int i = 0; i < nx; ++i) {
        double val = (i == prob.g.i0)
            ? (prob.p.K * std::exp(0.5 * prob.g.dx) - prob.p.K - 0.5 * prob.p.K * prob.g.dx) / prob.g.dx
            : std::max(std::exp(prob.g.x_min + i * prob.g.dx) - prob.p.K, 0.0);
        for (int j = 0; j < nv; ++j) U[i][j] = val;
    }

    for (int step = 0; step < m; ++step) {
        std::vector<std::vector<double>> A1(nx, std::vector<double>(nv, 0.0));
        std::vector<std::vector<double>> A2(nx, std::vector<double>(nv, 0.0));
        std::vector<std::vector<double>> W(nx, std::vector<double>(nv, 0.0));
        std::vector<std::vector<double>> Y(nx, std::vector<double>(nv, 0.0));

        for (int j = 0; j < nv; ++j) {
            for (int i = 0; i < nx; ++i) {
                double a0 = (i > 0 && i + 1 < nx && j > 0 && j + 1 < nv)
                    ? mixed(prob.m0[j], U[i+1][j+1], U[i-1][j+1], U[i+1][j-1], U[i-1][j-1]) : 0.0;
                double a1 = (i == 0) ? bi(prob.ax_b0_d, prob.ax_b0_u, U[0][j], U[1][j])
                          : (i == nx - 1) ? bi(prob.ax_bn_l, prob.ax_bn_d, U[nx-2][j], U[nx-1][j])
                          : tri(prob.ax_l[j], prob.ax_d[j], prob.ax_u[j], U[i-1][j], U[i][j], U[i+1][j]);
                double a2 = (j == 0) ? tri(prob.av_d[0], prob.av_u[0], prob.av_e0, U[i][0], U[i][1], U[i][2])
                          : (j == nv - 1) ? bi(prob.av_l[j], prob.av_d[j], U[i][j-1], U[i][j])
                          : tri(prob.av_l[j], prob.av_d[j], prob.av_u[j], U[i][j-1], U[i][j], U[i][j+1]);
                A1[i][j] = a1; A2[i][j] = a2;
                double y0 = ((U[i][j] + a0) + a1) + a2;
                W[i][j] = y0 - 0.5 * a0;
                Y[i][j] = y0 - 0.5 * a1;
            }
        }
        for (int j = 0; j < nv; ++j) {
            std::vector<double> lo(nx, prob.xl[j]), di(nx, prob.xd[j]), up(nx, prob.xu[j]), rhs(nx), sol(nx);
            di[0] = prob.xb0_d; up[0] = prob.xb0_u; lo[nx-1] = prob.xbn_l; di[nx-1] = prob.xbn_d;
            for (int i = 0; i < nx; ++i) rhs[i] = Y[i][j];
            thomas_div(nx, lo.data(), di.data(), up.data(), rhs.data(), sol.data());
            for (int i = 0; i < nx; ++i) U[i][j] = sol[i] - 0.5 * A2[i][j];
        }
        for (int i = 0; i < nx; ++i) {
            std::vector<double> rhs(nv), sol(nv);
            for (int j = 0; j < nv; ++j) rhs[j] = U[i][j];
            rhs[0] -= prob.f0 * rhs[1];
            thomas_div(nv, prob.vl.data(), prob.vd.data(), prob.vu.data(), rhs.data(), sol.data());
            for (int j = 0; j < nv; ++j) Y[i][j] = sol[j];
        }
        for (int j = 0; j < nv; ++j) {
            std::vector<double> lo(nx, prob.xl[j]), di(nx, prob.xd[j]), up(nx, prob.xu[j]), rhs(nx), sol(nx);
            di[0] = prob.xb0_d; up[0] = prob.xb0_u; lo[nx-1] = prob.xbn_l; di[nx-1] = prob.xbn_d;
            for (int i = 0; i < nx; ++i) {
                double b0 = (i > 0 && i + 1 < nx && j > 0 && j + 1 < nv)
                    ? mixed(prob.m0[j], Y[i+1][j+1], Y[i-1][j+1], Y[i+1][j-1], Y[i-1][j-1]) : 0.0;
                rhs[i] = (W[i][j] + 0.5 * b0) - 0.5 * A1[i][j];
            }
            thomas_div(nx, lo.data(), di.data(), up.data(), rhs.data(), sol.data());
            for (int i = 0; i < nx; ++i) U[i][j] = sol[i] - 0.5 * A2[i][j];
        }
        for (int i = 0; i < nx; ++i) {
            std::vector<double> rhs(nv), sol(nv);
            for (int j = 0; j < nv; ++j) rhs[j] = U[i][j];
            rhs[0] -= prob.f0 * rhs[1];
            thomas_div(nv, prob.vl.data(), prob.vd.data(), prob.vu.data(), rhs.data(), sol.data());
            for (int j = 0; j < nv; ++j) U[i][j] = sol[j];
        }
    }
    out_dense.resize(size_t(nv) * nx);
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nx; ++i) out_dense[size_t(j) * nx + i] = U[i][j];
}
