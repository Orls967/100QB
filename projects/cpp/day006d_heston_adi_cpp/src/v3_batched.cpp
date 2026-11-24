#include "variants.hpp"

void solve_v3(const Problem& prob, std::vector<double>& out_dense, const Opts& opt) {
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
        v_solve_batched<false>(nv, ld, 0, nx, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, Y, nullptr);
        for (int j = 0; j < nv; ++j) {
            eval_p3_row(prob, j, Y, W, A1, U);
            for (int i = 0; i < nx; ++i) U[size_t(j) * ld + i] -= 0.5 * A2[size_t(j) * ld + i];
        }
        v_solve_batched<false>(nv, ld, 0, nx, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, U, nullptr);
    }
    out_dense.resize(size_t(nv) * nx);
    pack_dense(U, out_dense.data(), prob.g);
}
