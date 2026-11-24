#include "check.hpp"
#include "reference.hpp"
#include "variants.hpp"
#include <cmath>
#include <cstdio>

namespace {
std::vector<double> gauss_dense(std::vector<double> A, std::vector<double> b, int n) {
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int r = k + 1; r < n; ++r) if (std::abs(A[r * n + k]) > std::abs(A[piv * n + k])) piv = r;
        for (int c = k; c < n; ++c) std::swap(A[k * n + c], A[piv * n + c]);
        std::swap(b[k], b[piv]);
        for (int r = k + 1; r < n; ++r) {
            double f = A[r * n + k] / A[k * n + k];
            for (int c = k; c < n; ++c) A[r * n + c] -= f * A[k * n + c];
            b[r] -= f * b[k];
        }
    }
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i]; for (int c = i + 1; c < n; ++c) s -= A[i * n + c] * x[c];
        x[i] = s / A[i * n + i];
    }
    return x;
}
} // namespace

TEST_CASE(gl16_exact_poly31) {
    GL16 gl; double sum = 0.0;
    for (int q = 0; q < 16; ++q) sum += gl.w[q] * (std::pow(gl.x[q], 30) + std::pow(gl.x[q], 31) + 1.0);
    CHECK(std::abs(sum - (2.0 / 31.0 + 2.0)) <= 1e-14);
}
TEST_CASE(bs_quadrature_match) {
    GL16 gl; BsPsi psi{100.0, 0.03, 1.0, 0.2}; double bs = bs_call(100.0, 100.0, 0.03, 1.0, 0.2);
    CHECK(std::abs(gp_call(psi, 100.0, 100.0, 0.03, 1.0, gl) - bs) <= 1e-10);
    CHECK(std::abs(lewis_call(psi, 100.0, 100.0, 0.03, 1.0, gl) - bs) <= 1e-10);
}
TEST_CASE(heston_martingale_and_parity) {
    GL16 gl;
    for (const Params& p : {Params::set_a(), Params::set_b()}) {
        HestonPsi psi{p, p.S0, p.T, p.v0};
        CHECK(std::abs(psi(std::complex<double>(0.0, -1.0)) / (p.S0 * std::exp(p.r * p.T)) - 1.0) <= 1e-12);
        double c_lew = lewis_call(psi, p.S0, p.K, p.r, p.T, gl), p_gp = gp_put(psi, p.S0, p.K, p.r, p.T, gl);
        CHECK(std::abs(c_lew - (p_gp + p.S0 - p.K * std::exp(-p.r * p.T))) <= 1e-8);
    }
}
TEST_CASE(heston_sigma_zero_limit) {
    GL16 gl; double bs = bs_call(100.0, 100.0, 0.03, 1.0, 0.2), prev = 1e9;
    for (double sig : {0.02, 0.01, 0.005}) {
        Params p = Params::set_a(); p.rho = 0.0; p.sigma = sig;
        double err = std::abs(gp_call(HestonPsi{p, p.S0, p.T, p.v0}, p.S0, p.K, p.r, p.T, gl) - bs);
        CHECK(err < prev); prev = err;
    }
    CHECK(prev <= 1e-4);
}
TEST_CASE(thomas_vs_dense_gauss) {
    const int n = 20; uint64_t rng = 20261124ULL;
    auto next_u = [&]() { rng = rng * 6364136223846793005ULL + 1442695040888963407ULL; return ((rng >> 33) + 1) / double(1ULL << 31); };
    std::vector<double> A(n * n, 0.0), r(n), lo(n, 0.0), di(n), up(n, 0.0), cp(n), y(n), xt(n);
    for (int i = 0; i < n; ++i) {
        if (i > 0) { lo[i] = -0.2 * next_u(); A[i * n + i - 1] = lo[i]; }
        if (i + 1 < n) { up[i] = -0.2 * next_u(); A[i * n + i + 1] = up[i]; }
        di[i] = 2.0 + next_u(); A[i * n + i] = di[i]; r[i] = next_u() - 0.5;
    }
    auto xd = gauss_dense(A, r, n);
    cp[0] = up[0] / di[0]; y[0] = r[0] / di[0];
    for (int i = 1; i < n; ++i) { double d = di[i] - lo[i] * cp[i-1]; cp[i] = up[i] / d; y[i] = (r[i] - lo[i] * y[i-1]) / d; }
    xt[n-1] = y[n-1]; for (int i = n - 2; i >= 0; --i) xt[i] = y[i] - cp[i] * xt[i+1];
    for (int i = 0; i < n; ++i) CHECK(std::abs(xt[i] - xd[i]) <= 1e-12);
    Problem prob(Params::set_a(), 15, 15, 20); int nv = 15;
    std::vector<double> Av(nv * nv, 0.0), rv(nv), xlu(nv);
    Av[0] = 1.0 - 0.5 * prob.av_d[0]; Av[1] = -0.5 * prob.av_u[0]; Av[2] = -0.5 * prob.av_e0;
    for (int j = 1; j < nv; ++j) { Av[j * nv + j - 1] = prob.vl[j]; Av[j * nv + j] = prob.vd[j]; if (j + 1 < nv) Av[j * nv + j + 1] = prob.vu[j]; }
    for (int j = 0; j < nv; ++j) rv[j] = next_u();
    auto xdv = gauss_dense(Av, rv, nv);
    rv[0] -= prob.f0 * rv[1]; xlu[0] = rv[0] * prob.vinv[0];
    for (int j = 1; j < nv; ++j) xlu[j] = (rv[j] - prob.vl[j] * xlu[j - 1]) * prob.vinv[j];
    for (int j = nv - 2; j >= 0; --j) xlu[j] -= prob.vcp[j] * xlu[j + 1];
    for (int j = 0; j < nv; ++j) CHECK(std::abs(xlu[j] - xdv[j]) <= 1e-12);
}
TEST_CASE(v0_v1_v2_match_and_pde_sanity) {
    Problem prob(Params::set_a(), 37, 21, 10); std::vector<double> u0, u1, u2;
    solve_v0(prob, u0); solve_v1(prob, u1); solve_v2(prob, u2);
    double mx = 0.0, d01 = 0.0, d02 = 0.0;
    for (size_t k = 0; k < u0.size(); ++k) { mx = std::max(mx, std::abs(u0[k])); d01 = std::max(d01, std::abs(u1[k] - u0[k])); d02 = std::max(d02, std::abs(u2[k] - u0[k])); }
    CHECK(d01 / mx <= 1e-12); CHECK(d02 / mx <= 1e-12);
}
int main() {
    for (auto& t : TestReg::all()) t.fn();
    std::printf("Tests run: %zu, failures: %d\n", TestReg::all().size(), TestReg::fails());
    return TestReg::fails() ? 1 : 0;
}
