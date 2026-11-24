#include "variants.hpp"
#include <barrier>
#include <thread>

void solve_v5(const Problem& prob, std::vector<double>& out_dense, const Opts& opt) {
    const int nx = prob.g.nx, nv = prob.g.nv, ld = prob.g.ld, m = prob.g.m, T = std::max(1, opt.threads);
    const size_t tot = size_t(nv) * ld;
    AlignedBuf u_b(tot), y_b(tot), w_b(tot), a1_b(tot), a2_b(tot);
    double *U = u_b.get(), *Y = y_b.get(), *W = w_b.get(), *A1 = a1_b.get(), *A2 = a2_b.get();
    init_payoff_strided(U, prob.g, prob.p);
    std::barrier sync(T);
    const int n_chunks = (nx + 15) / 16;
    auto worker = [&](int tid) {
        const int j_beg = (tid * nv) / T, j_end = ((tid + 1) * nv) / T;
        const int i_beg = ((tid * n_chunks) / T) * 16, i_end = std::min(nx, (((tid + 1) * n_chunks) / T) * 16);
        for (int step = 0; step < m; ++step) {
            for (int j = j_beg; j < j_end; ++j) eval_p1_row(prob, j, U, Y, W, A1, A2);
            sync.arrive_and_wait();
            if (i_beg < i_end)
                v_solve_batched<true>(nv, ld, i_beg, i_end, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, Y, A2);
            sync.arrive_and_wait();
            for (int j = j_beg; j < j_end; ++j) eval_p3_row(prob, j, Y, W, A1, U);
            sync.arrive_and_wait();
            if (i_beg < i_end)
                v_solve_batched<true>(nv, ld, i_beg, i_end, opt.B, prob.vl.data(), prob.vinv.data(), prob.vcp.data(), prob.f0, U, A2);
            sync.arrive_and_wait();
        }
    };
    std::vector<std::thread> ts;
    for (int t = 1; t < T; ++t) ts.emplace_back(worker, t);
    worker(0);
    for (auto& th : ts) th.join();
    out_dense.resize(size_t(nv) * nx);
    pack_dense(U, out_dense.data(), prob.g);
}
