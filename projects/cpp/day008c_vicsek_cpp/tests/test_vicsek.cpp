#include "../src/vicsek.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

static void test_rng_python_reference() {
    vicsek::Xoshiro256pp rng(vicsek::MASTER_SEED);
    assert(rng.next_u64() == 0xA9920E0E7DAB3EF0ULL);
    assert(rng.next_u64() == 0x033B9CEB4B3BA7A0ULL);
}

static void test_minimum_image_edge() {
    double L = 16.0;
    assert(std::fabs(vicsek::min_img(15.8 - 0.1, L) - (-0.3)) < 1e-12);
    assert(std::fabs(vicsek::min_img(0.1 - 15.8, L) - 0.3) < 1e-12);
}

static void test_cell_list_equals_bruteforce() {
    vicsek::Sim s_bf(200, 8.0, 1.5, vicsek::MASTER_SEED);
    vicsek::Sim s_cl(200, 8.0, 1.5, vicsek::MASTER_SEED);
    for (int t = 0; t < 50; ++t) {
        s_bf.step_bruteforce();
        s_cl.step_cell_list(true);
    }
    for (int i = 0; i < 200; ++i) {
        assert(s_bf.x[i] == s_cl.x[i]);
        assert(s_bf.y[i] == s_cl.y[i]);
        assert(s_bf.theta[i] == s_cl.theta[i]);
    }
}

static void test_zero_noise_aligned_stays_one() {
    vicsek::Sim s(100, 10.0, 0.0, vicsek::MASTER_SEED, true);
    assert(std::fabs(s.polarization() - 1.0) < 1e-12);
    for (int t = 0; t < 50; ++t) s.step_cell_list(true);
    assert(std::fabs(s.polarization() - 1.0) < 1e-12);
}

static void test_circular_mean_hand_case() {
    double a1 = 0.25, a2 = -0.25;
    double avg = std::atan2(std::sin(a1) + std::sin(a2), std::cos(a1) + std::cos(a2));
    assert(std::fabs(avg) < 1e-12);
    double b1 = 3.0, b2 = -3.0; /* across +/- pi */
    double avg2 = std::atan2(std::sin(b1) + std::sin(b2), std::cos(b1) + std::cos(b2));
    assert(std::fabs(std::fabs(avg2) - vicsek::PI) < 1e-12);
}

int main() {
    test_rng_python_reference();
    test_minimum_image_edge();
    test_cell_list_equals_bruteforce();
    test_zero_noise_aligned_stays_one();
    test_circular_mean_hand_case();
    std::printf("Vicsek C++ unit tests PASSED (5/5)\n");
    return 0;
}
