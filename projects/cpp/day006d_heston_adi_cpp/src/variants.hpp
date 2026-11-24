#pragma once
#include "common.hpp"
#include <algorithm>
#include <vector>

struct Opts { int B = 16, tile = 16, threads = 1; };

inline void x_solve_lu_row(int nx, double xl, double bnl, const double* __restrict xinv,
                           const double* __restrict xcp, double* __restrict row) {
    row[0] *= xinv[0];
    for (int i = 1; i < nx - 1; ++i) row[i] = (row[i] - xl * row[i - 1]) * xinv[i];
    row[nx - 1] = (row[nx - 1] - bnl * row[nx - 2]) * xinv[nx - 1];
    for (int i = nx - 2; i >= 0; --i) row[i] -= xcp[i] * row[i + 1];
}

inline void eval_p1_row(const Problem& prob, int j, const double* __restrict U, double* __restrict Y,
                        double* __restrict W, double* __restrict A1, double* __restrict A2) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld;
    const double* __restrict u = U + size_t(j) * ld;
    const double* __restrict um = (j > 0) ? u - ld : u;
    const double* __restrict up = (j + 1 < nv) ? u + ld : u;
    double *__restrict y = Y + size_t(j) * ld, *__restrict w = W + size_t(j) * ld;
    double *__restrict a1r = A1 + size_t(j) * ld, *__restrict a2r = A2 + size_t(j) * ld;
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
        w[i] = y0 - 0.5 * a0; y[i] = y0 - 0.5 * a1;
    }
    x_solve_lu_row(nx, prob.xl[j], prob.xbn_l, prob.xinv_buf.get() + size_t(j) * ld, prob.xcp_buf.get() + size_t(j) * ld, y);
}

inline void eval_p3_row(const Problem& prob, int j, const double* __restrict Y,
                        const double* __restrict W, const double* __restrict A1, double* __restrict U) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld;
    const double* __restrict y = Y + size_t(j) * ld;
    const double* __restrict ym = (j > 0) ? y - ld : y;
    const double* __restrict yp = (j + 1 < nv) ? y + ld : y;
    const double *__restrict w = W + size_t(j) * ld, *__restrict a1r = A1 + size_t(j) * ld;
    double* __restrict u = U + size_t(j) * ld;
    for (int i = 0; i < nx; ++i) {
        double b0 = (i > 0 && i + 1 < nx && j > 0 && j + 1 < nv)
            ? mixed(prob.m0[j], yp[i+1], yp[i-1], ym[i+1], ym[i-1]) : 0.0;
        u[i] = (w[i] + 0.5 * b0) - 0.5 * a1r[i];
    }
    x_solve_lu_row(nx, prob.xl[j], prob.xbn_l, prob.xinv_buf.get() + size_t(j) * ld, prob.xcp_buf.get() + size_t(j) * ld, u);
}

template <bool FuseA2>
inline void v_solve_batched(int nv, int ld, int i_beg, int i_end, int B, const double* __restrict vl,
                            const double* __restrict vinv, const double* __restrict vcp, double f0,
                            double* __restrict grid, const double* __restrict a2) {
    for (int ib = i_beg; ib < i_end; ib += B) {
        const int ie = std::min(i_end, ib + B);
        double *__restrict g0 = grid, *__restrict g1 = grid + ld, inv0 = vinv[0];
        if constexpr (FuseA2) {
            const double *__restrict a0 = a2, *__restrict a1 = a2 + ld;
            for (int i = ib; i < ie; ++i) g0[i] = ((g0[i] - 0.5 * a0[i]) - f0 * (g1[i] - 0.5 * a1[i])) * inv0;
        } else {
            for (int i = ib; i < ie; ++i) g0[i] = (g0[i] - f0 * g1[i]) * inv0;
        }
        for (int j = 1; j < nv; ++j) {
            double* __restrict gj = grid + size_t(j) * ld;
            const double* __restrict gp = grid + size_t(j - 1) * ld;
            const double l = vl[j], inv = vinv[j];
            if constexpr (FuseA2) {
                const double* __restrict aj = a2 + size_t(j) * ld;
                for (int i = ib; i < ie; ++i) gj[i] = ((gj[i] - 0.5 * aj[i]) - l * gp[i]) * inv;
            } else {
                for (int i = ib; i < ie; ++i) gj[i] = (gj[i] - l * gp[i]) * inv;
            }
        }
        for (int j = nv - 2; j >= 0; --j) {
            double* __restrict gj = grid + size_t(j) * ld;
            const double* __restrict gn = grid + size_t(j + 1) * ld;
            const double c = vcp[j];
            for (int i = ib; i < ie; ++i) gj[i] -= c * gn[i];
        }
    }
}

void solve_v0(const Problem& prob, std::vector<double>& out_dense);
void solve_v1(const Problem& prob, std::vector<double>& out_dense);
void solve_v2(const Problem& prob, std::vector<double>& out_dense);
void solve_v3(const Problem& prob, std::vector<double>& out_dense, const Opts& opt);
void solve_v4(const Problem& prob, std::vector<double>& out_dense, const Opts& opt);
void solve_v5(const Problem& prob, std::vector<double>& out_dense, const Opts& opt);
