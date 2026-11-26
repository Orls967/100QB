#pragma once

#include "lbm/collision.hpp"
#include "lbm/equilibrium.hpp"
#include "lbm/lattice.hpp"
#include "lbm/moments.hpp"
#include "lbm/rng.hpp"
#include <cstdint>
#include <vector>

namespace lbm {

enum class Scheme { PUSH_2GRID, PULL_2GRID, AA_PATTERN };
enum class MemLayout { AOS, SOA, TILE };

inline constexpr const char* scheme_name(Scheme s) noexcept {
    switch (s) {
        case Scheme::PUSH_2GRID: return "two_grid_push";
        case Scheme::PULL_2GRID: return "two_grid_pull";
        case Scheme::AA_PATTERN: return "aa_pattern";
    }
    return "two_grid_push";
}

inline constexpr const char* layout_name(MemLayout l) noexcept {
    switch (l) {
        case MemLayout::AOS: return "AoS";
        case MemLayout::SOA: return "SoA";
        case MemLayout::TILE: return "Tile";
    }
    return "SoA";
}

// Unified 2D D2Q9 simulation state supporting AoS, SoA, Blocked Tile, Two-Grid Push/Pull, and In-Place AA-pattern.
struct Grid2D {
    int nx = 0;
    int ny = 0;
    int tile_w = 16;
    int tile_h = 16;
    Scheme scheme = Scheme::PULL_2GRID;
    MemLayout layout = MemLayout::SOA;
    std::size_t step_count = 0;

    std::vector<double> f0;
    std::vector<double> f1;
    std::vector<std::uint8_t> solid;  // 0 = fluid, 1 = solid

    Grid2D() = default;
    Grid2D(int nx_, int ny_, Scheme sc = Scheme::PULL_2GRID, MemLayout lay = MemLayout::SOA,
           int tw = 16, int th = 16);

    std::size_t num_cells() const noexcept {
        return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);
    }

    inline std::size_t idx(int x, int y, int q) const noexcept {
        const std::size_t cell = static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                                 static_cast<std::size_t>(x);
        if (layout == MemLayout::AOS) {
            return cell * 9ULL + static_cast<std::size_t>(q);
        }
        if (layout == MemLayout::SOA) {
            return static_cast<std::size_t>(q) * num_cells() + cell;
        }
        // Blocked tile layout: tile index * (tile_area * 9) + SoA-within-tile
        const int ntx = (nx + tile_w - 1) / tile_w;
        const int tx = x / tile_w;
        const int ty = y / tile_h;
        const int lx = x % tile_w;
        const int ly = y % tile_h;
        const std::size_t tile_area = static_cast<std::size_t>(tile_w) * static_cast<std::size_t>(tile_h);
        const std::size_t tid = static_cast<std::size_t>(ty) * static_cast<std::size_t>(ntx) +
                                static_cast<std::size_t>(tx);
        const std::size_t lcell = static_cast<std::size_t>(ly) * static_cast<std::size_t>(tile_w) +
                                  static_cast<std::size_t>(lx);
        return (tid * 9ULL + static_cast<std::size_t>(q)) * tile_area + lcell;
    }

    // Logical population read (accounts for AA-pattern odd-step slot permutation if step_count is odd)
    double get_pop(int x, int y, int q) const noexcept;
    void set_pop(int x, int y, int q, double val) noexcept;

    void init_uniform(double rho, double ux, double uy);
    void set_cell_array(int x, int y, const std::array<double, 9>& vals);
    std::array<double, 9> get_cell_array(int x, int y) const noexcept;

    Macro2D cell_macro(int x, int y, double fx = 0.0, double fy = 0.0) const noexcept;
    std::uint64_t field_hash_fnv1a() const noexcept;
    void totals(double& out_mass, double& out_px, double& out_py, double& out_ek,
                double& out_rho_min, double& out_rho_max, double& out_u_max,
                bool& out_has_nan, double fx = 0.0, double fy = 0.0) const noexcept;
};

// 3D D3Q19 state supporting Push/Pull/AA-pattern & deterministic std::thread execution
struct Grid3D {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    Scheme scheme = Scheme::PULL_2GRID;
    std::size_t step_count = 0;
    std::vector<double> f0;
    std::vector<double> f1;
    std::vector<std::uint8_t> solid;

    Grid3D() = default;
    Grid3D(int nx_, int ny_, int nz_, Scheme sc = Scheme::PULL_2GRID);

    std::size_t num_cells() const noexcept {
        return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
               static_cast<std::size_t>(nz);
    }

    inline std::size_t idx(int x, int y, int z, int q) const noexcept {
        const std::size_t cell =
            (static_cast<std::size_t>(z) * static_cast<std::size_t>(ny) +
             static_cast<std::size_t>(y)) *
                static_cast<std::size_t>(nx) +
            static_cast<std::size_t>(x);
        return static_cast<std::size_t>(q) * num_cells() + cell;
    }

    double get_pop(int x, int y, int z, int q) const noexcept;
    void set_pop(int x, int y, int z, int q, double val) noexcept;
    void init_uniform(double rho, double ux, double uy, double uz);
    std::uint64_t field_hash_fnv1a() const noexcept;
    void totals(double& out_mass, double& out_px, double& out_py, double& out_pz,
                double& out_ek, double& out_rho_min, bool& out_has_nan,
                double fx = 0.0, double fy = 0.0, double fz = 0.0) const noexcept;
};

// Step kernels (implemented in src/collision_stream.cpp)
void step_periodic_2d(Grid2D& g, const CollisionParams& cp, int num_threads = 1);
void step_periodic_3d(Grid3D& g, const CollisionParams& cp, int num_threads = 1);

// Channel / enclosure 2D steppers supporting bounce-back and Zou-He boundaries
void step_channel_bb_2d(Grid2D& g, const CollisionParams& cp, bool periodic_x, bool halfway_bb,
                        double top_ux = 0.0, double bot_ux = 0.0, int num_threads = 1);

}  // namespace lbm
