#include "../src/wator.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

static void test_rng_python_reference() {
    wator::Xoshiro256pp rng(wator::MASTER_SEED);
    assert(rng.next_u64() == 0xA9920E0E7DAB3EF0ULL);
    assert(rng.next_u64() == 0x033B9CEB4B3BA7A0ULL);
}

static void test_torus_wrap() {
    wator::Params p; p.W = 10; p.H = 10;
    wator::World w(p, wator::MASTER_SEED);
    assert(w.wrap_nbr(0, 2) == 9);    /* x=0 left -> x=9 */
    assert(w.wrap_nbr(9, 0) == 0);    /* x=9 right -> x=0 */
    assert(w.wrap_nbr(0, 3) == 90);   /* y=0 down -> y=9 */
    assert(w.wrap_nbr(90, 1) == 0);   /* y=9 up -> y=0 */
}

static void test_invariants_and_determinism() {
    wator::Params p; p.W = 40; p.H = 40;
    wator::World w1(p, wator::MASTER_SEED);
    wator::World w2(p, wator::MASTER_SEED);
    for (int t = 0; t < 200; ++t) {
        w1.step();
        w2.step();
        assert(w1.check_invariants());
        assert(w1.n_fish == w2.n_fish && w1.n_sharks == w2.n_sharks);
    }
    for (size_t i = 0; i < w1.grid.size(); ++i) {
        assert(w1.grid[i].kind == w2.grid[i].kind);
    }
}

static void test_lv_rk4_conservation() {
    wator::LVParams lv{0.1, 0.0002, 0.1, 0.00005};
    double F = 2500.0, S = 600.0;
    double H0 = wator::lv_first_integral(F, S, lv);
    for (int k = 0; k < 5000; ++k) {
        wator::rk4_lv_step(F, S, 0.02, lv);
    }
    double H1 = wator::lv_first_integral(F, S, lv);
    assert(std::fabs(H1 - H0) < 1e-6);
}

int main() {
    test_rng_python_reference();
    test_torus_wrap();
    test_invariants_and_determinism();
    test_lv_rk4_conservation();
    std::printf("Wa-Tor C++ unit tests PASSED (4/4)\n");
    return 0;
}
