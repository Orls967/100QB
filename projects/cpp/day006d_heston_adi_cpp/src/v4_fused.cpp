#include "variants.hpp"

void solve_v4(const Problem& prob, std::vector<double>& out_dense, const Opts& opt) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld, m = prob.g.m;
    const size_t tot = size_t(nv) * ld;
    AlignedBuf u_b(tot), y_b(tot), w_b(tot), a1_b(tot), a2_b(tot);
    double *U = u_b.get(), *Y = y_b.get(), *W = w_b.get(), *A1 = a1_b.get(), *A2 = a2_b.get();
    init_payoff_strided(U, prob.g, prob.p);
    for (int step = 0; step < m; ++step) {
        for (int jb = 0; jb < nv; jb += opt.tile)
            for (int j = jb; j < std::min(nv, jb + opt.tile); ++j) eval_p1_row(prob, j, U, Y, W, A1, A2);
        v_solve_batched<true>(nv, ld, 0, nx, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, Y, A2);
        for (int jb = 0; jb < nv; jb += opt.tile)
            for (int j = jb; j < std::min(nv, jb + opt.tile); ++j) eval_p3_row(prob, j, Y, W, A1, U);
        v_solve_batched<true>(nv, ld, 0, nx, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, U, A2);
    }
    out_dense.resize(size_t(nv) * nx);
    pack_dense(U, out_dense.data(), prob.g);
}
