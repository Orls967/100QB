#include "check.hpp"
#include "variants.hpp"

TEST_CASE(v3_v4_bit_identical_to_v2) {
    Problem prob(Params::set_a(), 51, 25, 15);
    std::vector<double> u2, u3a, u3b, u3c, u4;
    solve_v2(prob, u2);
    uint64_t h2 = fnv1a_hash(u2.data(), u2.size());
    solve_v3(prob, u3a, Opts{1, 16, 1});
    solve_v3(prob, u3b, Opts{4, 16, 1});
    solve_v3(prob, u3c, Opts{16, 16, 1});
    solve_v4(prob, u4, Opts{16, 8, 1});
    CHECK(fnv1a_hash(u3a.data(), u3a.size()) == h2);
    CHECK(fnv1a_hash(u3b.data(), u3b.size()) == h2);
    CHECK(fnv1a_hash(u3c.data(), u3c.size()) == h2);
    CHECK(fnv1a_hash(u4.data(), u4.size()) == h2);
}

TEST_CASE(v5_thread_invariance) {
    Problem prob(Params::set_a(), 51, 25, 15);
    std::vector<double> u4, u5;
    solve_v4(prob, u4, Opts{16, 8, 1});
    uint64_t h4 = fnv1a_hash(u4.data(), u4.size());
    for (int t : {1, 2, 3, 4}) {
        solve_v5(prob, u5, Opts{16, 8, t});
        CHECK(fnv1a_hash(u5.data(), u5.size()) == h4);
    }
}
