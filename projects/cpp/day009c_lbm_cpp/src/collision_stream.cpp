#include "lbm/grid.hpp"
#include "lbm/io.hpp"
#include "lbm/parallel.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace lbm {

Grid2D::Grid2D(int nx_, int ny_, Scheme sc, MemLayout lay, int tw, int th)
    : nx(nx_), ny(ny_), tile_w(tw), tile_h(th), scheme(sc), layout(lay), step_count(0) {
    std::size_t alloc_cells = num_cells();
    if (layout == MemLayout::TILE) {
        const int ntx = (nx + tile_w - 1) / tile_w;
        const int nty = (ny + tile_h - 1) / tile_h;
        alloc_cells = static_cast<std::size_t>(ntx) * static_cast<std::size_t>(nty) *
                      static_cast<std::size_t>(tile_w) * static_cast<std::size_t>(tile_h);
    }
    f0.assign(alloc_cells * 9ULL, 0.0);
    if (scheme != Scheme::AA_PATTERN) {
        f1.assign(alloc_cells * 9ULL, 0.0);
    }
    solid.assign(num_cells(), 0);
}

double Grid2D::get_pop(int x, int y, int q) const noexcept {
    if (scheme == Scheme::AA_PATTERN && (step_count & 1ULL)) {
        const int px = (x + D2Q9::cx[static_cast<std::size_t>(D2Q9::opp[static_cast<std::size_t>(q)])] + nx) % nx;
        const int py = (y + D2Q9::cy[static_cast<std::size_t>(D2Q9::opp[static_cast<std::size_t>(q)])] + ny) % ny;
        return f0[idx(px, py, D2Q9::opp[static_cast<std::size_t>(q)])];
    }
    return f0[idx(x, y, q)];
}

void Grid2D::set_pop(int x, int y, int q, double val) noexcept {
    if (scheme == Scheme::AA_PATTERN && (step_count & 1ULL)) {
        const int m = D2Q9::opp[static_cast<std::size_t>(q)];
        const int px = (x + D2Q9::cx[static_cast<std::size_t>(m)] + nx) % nx;
        const int py = (y + D2Q9::cy[static_cast<std::size_t>(m)] + ny) % ny;
        f0[idx(px, py, m)] = val;
        return;
    }
    f0[idx(x, y, q)] = val;
}

void Grid2D::init_uniform(double rho, double ux, double uy) {
    step_count = 0;
    const auto feq = compute_feq<D2Q9>(rho, ux, uy);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            set_cell_array(x, y, feq);
        }
    }
}

void Grid2D::set_cell_array(int x, int y, const std::array<double, 9>& vals) {
    for (int q = 0; q < 9; ++q) {
        set_pop(x, y, q, vals[static_cast<std::size_t>(q)]);
    }
}

std::array<double, 9> Grid2D::get_cell_array(int x, int y) const noexcept {
    std::array<double, 9> vals{};
    for (int q = 0; q < 9; ++q) {
        vals[static_cast<std::size_t>(q)] = get_pop(x, y, q);
    }
    return vals;
}

Macro2D Grid2D::cell_macro(int x, int y, double fx, double fy) const noexcept {
    const auto vals = get_cell_array(x, y);
    return compute_macro_2d(std::span<const double, 9>(vals), fx, fy);
}

std::uint64_t Grid2D::field_hash_fnv1a() const noexcept {
    std::uint64_t h = 14695981039346656037ULL;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            for (int q = 0; q < 9; ++q) {
                const double v = get_pop(x, y, q);
                h = fnv1a_bytes(&v, sizeof(double), h);
            }
        }
    }
    return h;
}

void Grid2D::totals(double& out_mass, double& out_px, double& out_py, double& out_ek,
                    double& out_rho_min, double& out_rho_max, double& out_u_max,
                    bool& out_has_nan, double fx, double fy) const noexcept {
    out_mass = 0.0;
    out_px = 0.0;
    out_py = 0.0;
    out_ek = 0.0;
    out_rho_min = 1e300;
    out_rho_max = -1e300;
    out_u_max = 0.0;
    out_has_nan = false;

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            if (solid[static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                      static_cast<std::size_t>(x)]) {
                continue;
            }
            const Macro2D mc = cell_macro(x, y, fx, fy);
            if (!std::isfinite(mc.rho) || !std::isfinite(mc.ux) || !std::isfinite(mc.uy)) {
                out_has_nan = true;
                return;
            }
            out_mass += mc.rho;
            out_px += mc.rho * mc.ux;
            out_py += mc.rho * mc.uy;
            const double u2 = mc.ux * mc.ux + mc.uy * mc.uy;
            out_ek += 0.5 * mc.rho * u2;
            out_rho_min = std::min(out_rho_min, mc.rho);
            out_rho_max = std::max(out_rho_max, mc.rho);
            out_u_max = std::max(out_u_max, std::sqrt(u2));
        }
    }
}

void step_periodic_2d(Grid2D& g, const CollisionParams& cp, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;

    if (g.scheme == Scheme::PUSH_2GRID) {
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            std::array<double, 9> loc{};
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) {
                        loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                    }
                    collide_d2q9_cell(loc.data(), cp);
                    for (int q = 0; q < 9; ++q) {
                        const int tx = (ix + D2Q9::cx[static_cast<std::size_t>(q)] + nx) % nx;
                        const int ty = (iy + D2Q9::cy[static_cast<std::size_t>(q)] + ny) % ny;
                        g.f1[g.idx(tx, ty, q)] = loc[static_cast<std::size_t>(q)];
                    }
                }
            }
        });
        g.f0.swap(g.f1);
        g.step_count += 1;
        return;
    }

    if (g.scheme == Scheme::PULL_2GRID) {
        // Two-phase Collide -> Pull-Stream so operator ordering (Stream o Collide) matches Push & AA bitwise
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            std::array<double, 9> loc{};
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                    collide_d2q9_cell(loc.data(), cp);
                    for (int q = 0; q < 9; ++q) g.f1[g.idx(ix, iy, q)] = loc[static_cast<std::size_t>(q)];
                }
            }
        });
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) {
                        const int sx = (ix - D2Q9::cx[static_cast<std::size_t>(q)] + nx) % nx;
                        const int sy = (iy - D2Q9::cy[static_cast<std::size_t>(q)] + ny) % ny;
                        g.f0[g.idx(ix, iy, q)] = g.f1[g.idx(sx, sy, q)];
                    }
                }
            }
        });
        g.step_count += 1;
        return;
    }

    // In-place single-grid AA-pattern (bit-identical to Push & Pull at every integer step)
    if ((g.step_count & 1ULL) == 0ULL) {
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            std::array<double, 9> loc{};
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) {
                        loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                    }
                    collide_d2q9_cell(loc.data(), cp);
                    for (int q = 0; q < 9; ++q) {
                        g.f0[g.idx(ix, iy, D2Q9::opp[static_cast<std::size_t>(q)])] =
                            loc[static_cast<std::size_t>(q)];
                    }
                }
            }
        });
    } else {
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            std::array<double, 9> loc{};
            std::array<std::size_t, 9> addrs{};
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) {
                        const int m = D2Q9::opp[static_cast<std::size_t>(q)];
                        const int mx = (ix + D2Q9::cx[static_cast<std::size_t>(m)] + nx) % nx;
                        const int my = (iy + D2Q9::cy[static_cast<std::size_t>(m)] + ny) % ny;
                        const std::size_t a = g.idx(mx, my, m);
                        addrs[static_cast<std::size_t>(m)] = a;
                        loc[static_cast<std::size_t>(q)] = g.f0[a];
                    }
                    collide_d2q9_cell(loc.data(), cp);
                    for (int m = 0; m < 9; ++m) {
                        g.f0[addrs[static_cast<std::size_t>(m)]] = loc[static_cast<std::size_t>(m)];
                    }
                }
            }
        });
    }
    g.step_count += 1;
}

void step_channel_bb_2d(Grid2D& g, const CollisionParams& cp, bool periodic_x, bool halfway_bb,
                        double top_ux, double bot_ux, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;

    if (halfway_bb) {
        // Collide all cells then pull-stream with virtual walls at y = -0.5 and y = ny - 0.5
        parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
            std::array<double, 9> loc{};
            for (std::size_t y = y0; y < y1; ++y) {
                const int iy = static_cast<int>(y);
                for (int ix = 0; ix < nx; ++ix) {
                    for (int q = 0; q < 9; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                    collide_d2q9_cell(loc.data(), cp);
                    for (int q = 0; q < 9; ++q) {
                        const int tx = ix + D2Q9::cx[static_cast<std::size_t>(q)];
                        const int ty = iy + D2Q9::cy[static_cast<std::size_t>(q)];
                        if (ty < 0) {
                            const double cw = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(q)]) * bot_ux;
                            const int op = D2Q9::opp[static_cast<std::size_t>(q)];
                            g.f1[g.idx(ix, iy, op)] =
                                loc[static_cast<std::size_t>(q)] -
                                6.0 * D2Q9::w[static_cast<std::size_t>(q)] * 1.0 * cw;
                        } else if (ty >= ny) {
                            const double cw = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(q)]) * top_ux;
                            const int op = D2Q9::opp[static_cast<std::size_t>(q)];
                            g.f1[g.idx(ix, iy, op)] =
                                loc[static_cast<std::size_t>(q)] -
                                6.0 * D2Q9::w[static_cast<std::size_t>(q)] * 1.0 * cw;
                        } else {
                            int wx = tx;
                            if (periodic_x) {
                                wx = (tx + nx) % nx;
                            } else if (tx < 0 || tx >= nx) {
                                continue;
                            }
                            g.f1[g.idx(wx, ty, q)] = loc[static_cast<std::size_t>(q)];
                        }
                    }
                }
            }
        });
        g.f0.swap(g.f1);
        g.step_count += 1;
        return;
    }

    // Full-way bounce-back on solid rows y=0 and y=ny-1
    parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
        std::array<double, 9> loc{};
        for (std::size_t y = y0; y < y1; ++y) {
            const int iy = static_cast<int>(y);
            for (int ix = 0; ix < nx; ++ix) {
                for (int q = 0; q < 9; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                if (iy == 0 || iy == ny - 1) {
                    std::array<double, 9> bb{};
                    for (int q = 0; q < 9; ++q) {
                        bb[static_cast<std::size_t>(D2Q9::opp[static_cast<std::size_t>(q)])] =
                            loc[static_cast<std::size_t>(q)];
                    }
                    loc = bb;
                } else {
                    collide_d2q9_cell(loc.data(), cp);
                }
                for (int q = 0; q < 9; ++q) {
                    const int tx = periodic_x ? ((ix + D2Q9::cx[static_cast<std::size_t>(q)] + nx) % nx)
                                              : (ix + D2Q9::cx[static_cast<std::size_t>(q)]);
                    const int ty = iy + D2Q9::cy[static_cast<std::size_t>(q)];
                    if (tx >= 0 && tx < nx && ty >= 0 && ty < ny) {
                        g.f1[g.idx(tx, ty, q)] = loc[static_cast<std::size_t>(q)];
                    }
                }
            }
        }
    });
    g.f0.swap(g.f1);
    g.step_count += 1;
}

// ---------- Grid3D D3Q19 ----------

Grid3D::Grid3D(int nx_, int ny_, int nz_, Scheme sc)
    : nx(nx_), ny(ny_), nz(nz_), scheme(sc), step_count(0) {
    f0.assign(num_cells() * 19ULL, 0.0);
    if (scheme != Scheme::AA_PATTERN) {
        f1.assign(num_cells() * 19ULL, 0.0);
    }
    solid.assign(num_cells(), 0);
}

double Grid3D::get_pop(int x, int y, int z, int q) const noexcept {
    if (scheme == Scheme::AA_PATTERN && (step_count & 1ULL)) {
        const auto m = static_cast<std::size_t>(D3Q19::opp[static_cast<std::size_t>(q)]);
        const int px = (x + D3Q19::cx[m] + nx) % nx;
        const int py = (y + D3Q19::cy[m] + ny) % ny;
        const int pz = (z + D3Q19::cz[m] + nz) % nz;
        return f0[idx(px, py, pz, static_cast<int>(m))];
    }
    return f0[idx(x, y, z, q)];
}

void Grid3D::set_pop(int x, int y, int z, int q, double val) noexcept {
    if (scheme == Scheme::AA_PATTERN && (step_count & 1ULL)) {
        const auto m = static_cast<std::size_t>(D3Q19::opp[static_cast<std::size_t>(q)]);
        const int px = (x + D3Q19::cx[m] + nx) % nx;
        const int py = (y + D3Q19::cy[m] + ny) % ny;
        const int pz = (z + D3Q19::cz[m] + nz) % nz;
        f0[idx(px, py, pz, static_cast<int>(m))] = val;
        return;
    }
    f0[idx(x, y, z, q)] = val;
}

void Grid3D::init_uniform(double rho, double ux, double uy, double uz) {
    step_count = 0;
    const auto feq = compute_feq<D3Q19>(rho, ux, uy, uz);
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                for (int q = 0; q < 19; ++q) {
                    f0[idx(x, y, z, q)] = feq[static_cast<std::size_t>(q)];
                }
            }
        }
    }
}

std::uint64_t Grid3D::field_hash_fnv1a() const noexcept {
    std::uint64_t h = 14695981039346656037ULL;
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                for (int q = 0; q < 19; ++q) {
                    const double v = get_pop(x, y, z, q);
                    h = fnv1a_bytes(&v, sizeof(double), h);
                }
            }
        }
    }
    return h;
}

void Grid3D::totals(double& out_mass, double& out_px, double& out_py, double& out_pz,
                    double& out_ek, double& out_rho_min, bool& out_has_nan,
                    double fx, double fy, double fz) const noexcept {
    out_mass = out_px = out_py = out_pz = out_ek = 0.0;
    out_rho_min = 1e300;
    out_has_nan = false;
    std::array<double, 19> loc{};
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = get_pop(x, y, z, q);
                const Macro3D mc = compute_macro_3d<D3Q19>(loc.data(), fx, fy, fz);
                if (!std::isfinite(mc.rho) || !std::isfinite(mc.ux)) {
                    out_has_nan = true;
                    return;
                }
                out_mass += mc.rho;
                out_px += mc.rho * mc.ux;
                out_py += mc.rho * mc.uy;
                out_pz += mc.rho * mc.uz;
                out_ek += 0.5 * mc.rho * (mc.ux * mc.ux + mc.uy * mc.uy + mc.uz * mc.uz);
                out_rho_min = std::min(out_rho_min, mc.rho);
            }
        }
    }
}

void step_periodic_3d(Grid3D& g, const CollisionParams& cp, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;
    const int nz = g.nz;

    if (g.scheme == Scheme::PUSH_2GRID) {
        parallel_for_range(static_cast<std::size_t>(nz), num_threads, [&](std::size_t z0, std::size_t z1, int) {
            std::array<double, 19> loc{};
            for (std::size_t z = z0; z < z1; ++z) {
                const int iz = static_cast<int>(z);
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, iz, q)];
                        collide_d3q19_cell(loc.data(), cp);
                        for (int q = 0; q < 19; ++q) {
                            const auto qq = static_cast<std::size_t>(q);
                            const int tx = (ix + D3Q19::cx[qq] + nx) % nx;
                            const int ty = (iy + D3Q19::cy[qq] + ny) % ny;
                            const int tz = (iz + D3Q19::cz[qq] + nz) % nz;
                            g.f1[g.idx(tx, ty, tz, q)] = loc[qq];
                        }
                    }
                }
            }
        });
        g.f0.swap(g.f1);
        g.step_count += 1;
        return;
    }

    if (g.scheme == Scheme::PULL_2GRID) {
        // Two-phase Collide -> Pull-Stream matching Push & AA bitwise
        parallel_for_range(static_cast<std::size_t>(nz), num_threads, [&](std::size_t z0, std::size_t z1, int) {
            std::array<double, 19> loc{};
            for (std::size_t z = z0; z < z1; ++z) {
                const int iz = static_cast<int>(z);
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, iz, q)];
                        collide_d3q19_cell(loc.data(), cp);
                        for (int q = 0; q < 19; ++q) g.f1[g.idx(ix, iy, iz, q)] = loc[static_cast<std::size_t>(q)];
                    }
                }
            }
        });
        parallel_for_range(static_cast<std::size_t>(nz), num_threads, [&](std::size_t z0, std::size_t z1, int) {
            for (std::size_t z = z0; z < z1; ++z) {
                const int iz = static_cast<int>(z);
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        for (int q = 0; q < 19; ++q) {
                            const auto qq = static_cast<std::size_t>(q);
                            const int sx = (ix - D3Q19::cx[qq] + nx) % nx;
                            const int sy = (iy - D3Q19::cy[qq] + ny) % ny;
                            const int sz = (iz - D3Q19::cz[qq] + nz) % nz;
                            g.f0[g.idx(ix, iy, iz, q)] = g.f1[g.idx(sx, sy, sz, q)];
                        }
                    }
                }
            }
        });
        g.step_count += 1;
        return;
    }

    // D3Q19 AA-pattern in-place
    if ((g.step_count & 1ULL) == 0ULL) {
        parallel_for_range(static_cast<std::size_t>(nz), num_threads, [&](std::size_t z0, std::size_t z1, int) {
            std::array<double, 19> loc{};
            for (std::size_t z = z0; z < z1; ++z) {
                const int iz = static_cast<int>(z);
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        for (int q = 0; q < 19; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, iz, q)];
                        collide_d3q19_cell(loc.data(), cp);
                        for (int q = 0; q < 19; ++q) {
                            g.f0[g.idx(ix, iy, iz, D3Q19::opp[static_cast<std::size_t>(q)])] =
                                loc[static_cast<std::size_t>(q)];
                        }
                    }
                }
            }
        });
    } else {
        parallel_for_range(static_cast<std::size_t>(nz), num_threads, [&](std::size_t z0, std::size_t z1, int) {
            std::array<double, 19> loc{};
            std::array<std::size_t, 19> addrs{};
            for (std::size_t z = z0; z < z1; ++z) {
                const int iz = static_cast<int>(z);
                for (int iy = 0; iy < ny; ++iy) {
                    for (int ix = 0; ix < nx; ++ix) {
                        for (int q = 0; q < 19; ++q) {
                            const int m = D3Q19::opp[static_cast<std::size_t>(q)];
                            const auto mm = static_cast<std::size_t>(m);
                            const int mx = (ix + D3Q19::cx[mm] + nx) % nx;
                            const int my = (iy + D3Q19::cy[mm] + ny) % ny;
                            const int mz = (iz + D3Q19::cz[mm] + nz) % nz;
                            const std::size_t a = g.idx(mx, my, mz, m);
                            addrs[mm] = a;
                            loc[static_cast<std::size_t>(q)] = g.f0[a];
                        }
                        collide_d3q19_cell(loc.data(), cp);
                        for (int m = 0; m < 19; ++m) {
                            g.f0[addrs[static_cast<std::size_t>(m)]] = loc[static_cast<std::size_t>(m)];
                        }
                    }
                }
            }
        });
    }
    g.step_count += 1;
}

}  // namespace lbm
