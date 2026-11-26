#include "lbm/boundary.hpp"
#include "lbm/parallel.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace lbm {

void apply_zou_he_velocity_west(Grid2D& g, int y_min, int y_max, double ux, double uy) {
    const int x = 0;
    for (int y = y_min; y <= y_max; ++y) {
        const double f0 = g.f0[g.idx(x, y, 0)];
        const double f2 = g.f0[g.idx(x, y, 2)];
        const double f3 = g.f0[g.idx(x, y, 3)];
        const double f4 = g.f0[g.idx(x, y, 4)];
        const double f6 = g.f0[g.idx(x, y, 6)];
        const double f7 = g.f0[g.idx(x, y, 7)];
        const double rho = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - ux);
        const double ru_x = rho * ux;
        const double ru_y = rho * uy;
        const double f1 = f3 + (2.0 / 3.0) * ru_x;
        const double f5 = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * ru_x + 0.5 * ru_y;
        const double f8 = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * ru_x - 0.5 * ru_y;
        g.f0[g.idx(x, y, 1)] = f1;
        g.f0[g.idx(x, y, 5)] = f5;
        g.f0[g.idx(x, y, 8)] = f8;
    }
}

void apply_zou_he_pressure_west(Grid2D& g, int y_min, int y_max, double rho_in, double uy) {
    const int x = 0;
    for (int y = y_min; y <= y_max; ++y) {
        const double f0 = g.f0[g.idx(x, y, 0)];
        const double f2 = g.f0[g.idx(x, y, 2)];
        const double f3 = g.f0[g.idx(x, y, 3)];
        const double f4 = g.f0[g.idx(x, y, 4)];
        const double f6 = g.f0[g.idx(x, y, 6)];
        const double f7 = g.f0[g.idx(x, y, 7)];
        const double ux = 1.0 - (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / rho_in;
        const double ru_x = rho_in * ux;
        const double ru_y = rho_in * uy;
        g.f0[g.idx(x, y, 1)] = f3 + (2.0 / 3.0) * ru_x;
        g.f0[g.idx(x, y, 5)] = f7 - 0.5 * (f2 - f4) + (1.0 / 6.0) * ru_x + 0.5 * ru_y;
        g.f0[g.idx(x, y, 8)] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * ru_x - 0.5 * ru_y;
    }
}

void apply_zou_he_pressure_east(Grid2D& g, int y_min, int y_max, double rho_out, double uy) {
    const int x = g.nx - 1;
    for (int y = y_min; y <= y_max; ++y) {
        const double f0 = g.f0[g.idx(x, y, 0)];
        const double f1 = g.f0[g.idx(x, y, 1)];
        const double f2 = g.f0[g.idx(x, y, 2)];
        const double f4 = g.f0[g.idx(x, y, 4)];
        const double f5 = g.f0[g.idx(x, y, 5)];
        const double f8 = g.f0[g.idx(x, y, 8)];
        const double ux = -1.0 + (f0 + f2 + f4 + 2.0 * (f1 + f5 + f8)) / rho_out;
        const double ru_x = rho_out * ux;
        const double ru_y = rho_out * uy;
        g.f0[g.idx(x, y, 3)] = f1 - (2.0 / 3.0) * ru_x;
        g.f0[g.idx(x, y, 7)] = f5 + 0.5 * (f2 - f4) - (1.0 / 6.0) * ru_x - 0.5 * ru_y;
        g.f0[g.idx(x, y, 6)] = f8 - 0.5 * (f2 - f4) - (1.0 / 6.0) * ru_x + 0.5 * ru_y;
    }
}

void apply_convective_outflow_east(Grid2D& g, double /*u_conv*/) {
    const int x = g.nx - 1;
    const int xm = g.nx - 2;
    for (int y = 0; y < g.ny; ++y) {
        for (int q : {3, 6, 7}) {
            g.f0[g.idx(x, y, q)] = g.f0[g.idx(xm, y, q)];
        }
    }
}

void step_cavity_2d(Grid2D& g, const CollisionParams& cp, double u_lid, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;
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
                    if (ty >= ny) {
                        // Moving top lid at y = ny - 0.5 with velocity (u_lid, 0)
                        const double cu = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(q)]) * u_lid;
                        const int op = D2Q9::opp[static_cast<std::size_t>(q)];
                        g.f1[g.idx(ix, iy, op)] =
                            loc[static_cast<std::size_t>(q)] - 6.0 * D2Q9::w[static_cast<std::size_t>(q)] * cu;
                    } else if (tx < 0 || tx >= nx || ty < 0) {
                        // Stationary half-way bounce-back on Left, Right, Bottom
                        const int op = D2Q9::opp[static_cast<std::size_t>(q)];
                        g.f1[g.idx(ix, iy, op)] = loc[static_cast<std::size_t>(q)];
                    } else {
                        g.f1[g.idx(tx, ty, q)] = loc[static_cast<std::size_t>(q)];
                    }
                }
            }
        }
    });
    g.f0.swap(g.f1);
    g.step_count += 1;
}

CylinderDomain build_cylinder_domain(Grid2D& g, double cx, double cy, double radius, bool use_bouzidi) {
    CylinderDomain dom;
    dom.cx = cx;
    dom.cy = cy;
    dom.radius = radius;
    dom.use_bouzidi = use_bouzidi;
    std::fill(g.solid.begin(), g.solid.end(), 0);

    for (int y = 0; y < g.ny; ++y) {
        for (int x = 0; x < g.nx; ++x) {
            const double dx = static_cast<double>(x) - cx;
            const double dy = static_cast<double>(y) - cy;
            if (dx * dx + dy * dy <= radius * radius) {
                g.solid[static_cast<std::size_t>(y) * static_cast<std::size_t>(g.nx) +
                        static_cast<std::size_t>(x)] = 1;
            }
        }
    }

    // Collect boundary links sorted deterministically by (fy, fx, dir)
    for (int y = 0; y < g.ny; ++y) {
        for (int x = 0; x < g.nx; ++x) {
            if (g.solid[static_cast<std::size_t>(y) * static_cast<std::size_t>(g.nx) +
                        static_cast<std::size_t>(x)]) {
                continue;
            }
            for (int q = 1; q < 9; ++q) {
                const int sx = x + D2Q9::cx[static_cast<std::size_t>(q)];
                const int sy = y + D2Q9::cy[static_cast<std::size_t>(q)];
                if (sx < 0 || sx >= g.nx || sy < 0 || sy >= g.ny) continue;
                if (!g.solid[static_cast<std::size_t>(sy) * static_cast<std::size_t>(g.nx) +
                             static_cast<std::size_t>(sx)]) {
                    continue;
                }
                double q_frac = 0.5;
                if (use_bouzidi) {
                    // Solve | (x, y) + q*(c_x, c_y) - (cx, cy) |^2 = radius^2 for q in (0, 1]
                    const double ox = static_cast<double>(x) - cx;
                    const double oy = static_cast<double>(y) - cy;
                    const double vx = static_cast<double>(D2Q9::cx[static_cast<std::size_t>(q)]);
                    const double vy = static_cast<double>(D2Q9::cy[static_cast<std::size_t>(q)]);
                    const double a = vx * vx + vy * vy;
                    const double b = 2.0 * (ox * vx + oy * vy);
                    const double c = ox * ox + oy * oy - radius * radius;
                    const double disc = std::max(0.0, b * b - 4.0 * a * c);
                    const double q1 = (-b - std::sqrt(disc)) / (2.0 * a);
                    const double q2 = (-b + std::sqrt(disc)) / (2.0 * a);
                    if (q1 > 0.0 && q1 <= 1.0001) {
                        q_frac = std::clamp(q1, 1e-4, 1.0);
                    } else {
                        q_frac = std::clamp(q2, 1e-4, 1.0);
                    }
                }
                dom.links.push_back({x, y, q, q_frac});
            }
        }
    }
    return dom;
}

ForceResult step_cylinder_2d(Grid2D& g, const CylinderDomain& cyl, const CollisionParams& cp,
                             double u_inlet, double rho_out, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;

    // 1. Collision at all fluid nodes (in-place into f0 before streaming)
    parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
        std::array<double, 9> loc{};
        for (std::size_t y = y0; y < y1; ++y) {
            const int iy = static_cast<int>(y);
            for (int ix = 0; ix < nx; ++ix) {
                if (g.solid[static_cast<std::size_t>(iy) * static_cast<std::size_t>(nx) +
                            static_cast<std::size_t>(ix)]) {
                    continue;
                }
                for (int q = 0; q < 9; ++q) loc[static_cast<std::size_t>(q)] = g.f0[g.idx(ix, iy, q)];
                collide_d2q9_cell(loc.data(), cp);
                for (int q = 0; q < 9; ++q) g.f0[g.idx(ix, iy, q)] = loc[static_cast<std::size_t>(q)];
            }
        }
    });

    // 2. Streaming on fluid nodes (periodic top/bottom, open west/east; links into solid skipped)
    parallel_for_range(static_cast<std::size_t>(ny), num_threads, [&](std::size_t y0, std::size_t y1, int) {
        for (std::size_t y = y0; y < y1; ++y) {
            const int iy = static_cast<int>(y);
            for (int ix = 0; ix < nx; ++ix) {
                if (g.solid[static_cast<std::size_t>(iy) * static_cast<std::size_t>(nx) +
                            static_cast<std::size_t>(ix)]) {
                    continue;
                }
                for (int q = 0; q < 9; ++q) {
                    const int sx = ix - D2Q9::cx[static_cast<std::size_t>(q)];
                    const int sy = (iy - D2Q9::cy[static_cast<std::size_t>(q)] + ny) % ny;
                    if (sx >= 0 && sx < nx &&
                        !g.solid[static_cast<std::size_t>(sy) * static_cast<std::size_t>(nx) +
                                 static_cast<std::size_t>(sx)]) {
                        g.f1[g.idx(ix, iy, q)] = g.f0[g.idx(sx, sy, q)];
                    }
                }
            }
        }
    });

    // 3. Curved Bouzidi / staircase bounce-back + deterministic force summation
    ForceResult res{};
    for (const auto& lnk : cyl.links) {
        const int dir = lnk.dir;
        const int op = D2Q9::opp[static_cast<std::size_t>(dir)];
        const double f_star_out = g.f0[g.idx(lnk.fx, lnk.fy, dir)];
        double f_bounced = f_star_out;
        if (cyl.use_bouzidi) {
            const int bx = lnk.fx - D2Q9::cx[static_cast<std::size_t>(dir)];
            const int by = (lnk.fy - D2Q9::cy[static_cast<std::size_t>(dir)] + ny) % ny;
            const bool bx_valid = (bx >= 0 && bx < nx) &&
                                  !g.solid[static_cast<std::size_t>(by) * static_cast<std::size_t>(nx) +
                                           static_cast<std::size_t>(bx)];
            if (lnk.q < 0.5 && bx_valid) {
                const double f_star_back = g.f0[g.idx(bx, by, dir)];
                f_bounced = 2.0 * lnk.q * f_star_out + (1.0 - 2.0 * lnk.q) * f_star_back;
            } else if (lnk.q >= 0.5) {
                const double f_star_opp = g.f0[g.idx(lnk.fx, lnk.fy, op)];
                f_bounced = (0.5 / lnk.q) * f_star_out + ((2.0 * lnk.q - 1.0) / (2.0 * lnk.q)) * f_star_opp;
            }
        }
        g.f1[g.idx(lnk.fx, lnk.fy, op)] = f_bounced;

        // Momentum exchange on solid body: c_dir * (f*_dir(xf) + f_opp(xf))
        const double mom = f_star_out + f_bounced;
        res.fx += static_cast<double>(D2Q9::cx[static_cast<std::size_t>(dir)]) * mom;
        res.fy += static_cast<double>(D2Q9::cy[static_cast<std::size_t>(dir)]) * mom;
    }

    g.f0.swap(g.f1);
    g.step_count += 1;

    // 4. Gentle downstream acoustic sponge layer in buffer [nx-14, nx-2] to damp standing waves
    const int sponge_start = std::max(1, nx - 14);
    for (int y = 0; y < ny; ++y) {
        for (int x = sponge_start; x < nx - 1; ++x) {
            if (g.solid[static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                        static_cast<std::size_t>(x)]) {
                continue;
            }
            const double frac = static_cast<double>(x - sponge_start) / static_cast<double>(nx - 1 - sponge_start);
            const double sig = 0.08 * frac * frac;
            const auto mc = g.cell_macro(x, y);
            const auto feq = compute_feq<D2Q9>(1.0, mc.ux, mc.uy);
            for (int q = 0; q < 9; ++q) {
                const double f_cur = g.f0[g.idx(x, y, q)];
                g.f0[g.idx(x, y, q)] = (1.0 - sig) * f_cur + sig * feq[static_cast<std::size_t>(q)];
            }
        }
    }

    // 5. Inlet Zou-He uniform velocity (west) & Outlet Zou-He pressure (east)
    apply_zou_he_velocity_west(g, 0, ny - 1, u_inlet, 0.0);
    apply_zou_he_pressure_east(g, 0, ny - 1, rho_out, 0.0);

    const double diameter = 2.0 * cyl.radius;
    const double q_dyn = 0.5 * 1.0 * u_inlet * u_inlet * diameter;
    res.cd = (q_dyn > 1e-30) ? (res.fx / q_dyn) : 0.0;
    res.cl = (q_dyn > 1e-30) ? (res.fy / q_dyn) : 0.0;
    return res;
}

void step_duct_3d(Grid3D& g, const CollisionParams& cp, int num_threads) {
    const int nx = g.nx;
    const int ny = g.ny;
    const int nz = g.nz;
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
                        const int ty = iy + D3Q19::cy[qq];
                        const int tz = iz + D3Q19::cz[qq];
                        if (ty < 0 || ty >= ny || tz < 0 || tz >= nz) {
                            const int op = D3Q19::opp[qq];
                            g.f1[g.idx(ix, iy, iz, op)] = loc[qq];
                        } else {
                            g.f1[g.idx(tx, ty, tz, q)] = loc[qq];
                        }
                    }
                }
            }
        }
    });
    g.f0.swap(g.f1);
    g.step_count += 1;
}

}  // namespace lbm
