#include "variants.hpp"

namespace {
inline void v_solve_lu(int nv, int ld, int i, const double* __restrict vl,
                       const double* __restrict vinv, const double* __restrict vcp,
                       double f0, double* __restrict grid) {
    grid[i] -= f0 * grid[ld + i];
    grid[i] *= vinv[0];
    for (int j = 1; j < nv; ++j)
        grid[size_t(j) * ld + i] = (grid[size_t(j) * ld + i] - vl[j] * grid[size_t(j - 1) * ld + i]) * vinv[j];
    for (int j = nv - 2; j >= 0; --j)
        grid[size_t(j) * ld + i] -= vcp[j] * grid[size_t(j + 1) * ld + i];
}
} // namespace

void solve_v2(const Problem& prob, std::vector<double>& out_dense) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld, m = prob.g.m;
    const size_t tot = size_t(nv) * ld;
    AlignedBuf u_b(tot), y_b(tot), w_b(tot), a1_b(tot), a2_b(tot);
    double *U = u_b.get(), *Y = y_b.get(), *W = w_b.get(), *A1 = a1_b.get(), *A2 = a2_b.get();
    init_payoff_strided(U, prob.g, prob.p);
    for (int step = 0; step < m; ++step) {
        for (int j = 0; j < nv; ++j) {
            eval_p1_row(prob, j, U, Y, W, A1, A2);
            for (int i = 0; i < nx; ++i) Y[size_t(j) * ld + i] -= 0.5 * A2[size_t(j) * ld + i];
        }
        for (int i = 0; i < nx; ++i) v_solve_lu(nv, ld, i, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, Y);
        for (int j = 0; j < nv; ++j) {
            eval_p3_row(prob, j, Y, W, A1, U);
            for (int i = 0; i < nx; ++i) U[size_t(j) * ld + i] -= 0.5 * A2[size_t(j) * ld + i];
        }
        for (int i = 0; i < nx; ++i) v_solve_lu(nv, ld, i, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, U);
    }
    out_dense.resize(size_t(nv) * nx);
    pack_dense(U, out_dense.data(), prob.g);
}
