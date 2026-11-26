#pragma once

#include <array>
#include <concepts>
#include <cstddef>

namespace lbm {

template <typename T>
concept LatticeModel = requires {
    { T::D } -> std::convertible_to<int>;
    { T::Q } -> std::convertible_to<int>;
    { T::CS2 } -> std::convertible_to<double>;
    { T::INV_CS2 } -> std::convertible_to<double>;
    { T::INV_CS4 } -> std::convertible_to<double>;
    { T::cx[0] } -> std::convertible_to<int>;
    { T::cy[0] } -> std::convertible_to<int>;
    { T::cz[0] } -> std::convertible_to<int>;
    { T::w[0] } -> std::convertible_to<double>;
    { T::opp[0] } -> std::convertible_to<int>;
};

struct D2Q9 {
    static constexpr int D = 2;
    static constexpr int Q = 9;
    static constexpr double CS2 = 1.0 / 3.0;
    static constexpr double CS4 = 1.0 / 9.0;
    static constexpr double INV_CS2 = 3.0;
    static constexpr double INV_CS4 = 9.0;
    static constexpr const char* NAME = "D2Q9";

    // Ordering: 0=center, 1=E, 2=N, 3=W, 4=S, 5=NE, 6=NW, 7=SW, 8=SE
    static constexpr std::array<int, Q> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr std::array<int, Q> cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr std::array<int, Q> cz = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    static constexpr std::array<int, Q> opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    static constexpr std::array<double, Q> w = {
        4.0 / 9.0,
        1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
};

struct D3Q19 {
    static constexpr int D = 3;
    static constexpr int Q = 19;
    static constexpr double CS2 = 1.0 / 3.0;
    static constexpr double CS4 = 1.0 / 9.0;
    static constexpr double INV_CS2 = 3.0;
    static constexpr double INV_CS4 = 9.0;
    static constexpr const char* NAME = "D3Q19";

    // 0: (0,0,0)
    // 1..6: (+-1,0,0), (0,+-1,0), (0,0,+-1) paired as (1,2), (3,4), (5,6)
    // 7..18: face diagonals paired sequentially (7,8), (9,10), (11,12), (13,14), (15,16), (17,18)
    static constexpr std::array<int, Q> cx = {
        0,
        1, -1, 0, 0, 0, 0,
        1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0
    };
    static constexpr std::array<int, Q> cy = {
        0,
        0, 0, 1, -1, 0, 0,
        1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1
    };
    static constexpr std::array<int, Q> cz = {
        0,
        0, 0, 0, 0, 1, -1,
        0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1
    };
    static constexpr std::array<int, Q> opp = {
        0,
        2, 1, 4, 3, 6, 5,
        8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17
    };
    static constexpr std::array<double, Q> w = {
        1.0 / 3.0,
        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };
};

namespace detail {

constexpr std::array<int, 27> make_d3q27_c(int axis) {
    std::array<int, 27> arr{};
    int idx = 0;
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                arr[static_cast<std::size_t>(idx++)] = (axis == 0 ? x : (axis == 1 ? y : z));
            }
        }
    }
    return arr;
}

constexpr std::array<int, 27> make_d3q27_opp() {
    std::array<int, 27> arr{};
    for (int i = 0; i < 27; ++i) arr[static_cast<std::size_t>(i)] = 26 - i;
    return arr;
}

constexpr std::array<double, 27> make_d3q27_w() {
    constexpr auto cx = make_d3q27_c(0);
    constexpr auto cy = make_d3q27_c(1);
    constexpr auto cz = make_d3q27_c(2);
    std::array<double, 27> arr{};
    for (int i = 0; i < 27; ++i) {
        const auto ii = static_cast<std::size_t>(i);
        const int s = (cx[ii] != 0 ? 1 : 0) + (cy[ii] != 0 ? 1 : 0) + (cz[ii] != 0 ? 1 : 0);
        arr[ii] = (s == 0) ? (8.0 / 27.0) : (s == 1) ? (2.0 / 27.0) : (s == 2) ? (1.0 / 54.0) : (1.0 / 216.0);
    }
    return arr;
}

}  // namespace detail

struct D3Q27 {
    static constexpr int D = 3;
    static constexpr int Q = 27;
    static constexpr double CS2 = 1.0 / 3.0;
    static constexpr double CS4 = 1.0 / 9.0;
    static constexpr double INV_CS2 = 3.0;
    static constexpr double INV_CS4 = 9.0;
    static constexpr const char* NAME = "D3Q27";

    static constexpr std::array<int, Q> cx = detail::make_d3q27_c(0);
    static constexpr std::array<int, Q> cy = detail::make_d3q27_c(1);
    static constexpr std::array<int, Q> cz = detail::make_d3q27_c(2);
    static constexpr std::array<int, Q> opp = detail::make_d3q27_opp();
    static constexpr std::array<double, Q> w = detail::make_d3q27_w();
};

static_assert(LatticeModel<D2Q9>);
static_assert(LatticeModel<D3Q19>);
static_assert(LatticeModel<D3Q27>);

struct IsotropyReport {
    double weight_sum_err = 0.0;
    double order1_err = 0.0;
    double order2_err = 0.0;
    double order3_err = 0.0;
    double order4_err = 0.0;
    double order5_err = 0.0;
    double opp_check_err = 0.0;
};

IsotropyReport check_d2q9_isotropy();
IsotropyReport check_d3q19_isotropy();
IsotropyReport check_d3q27_isotropy();

}  // namespace lbm
