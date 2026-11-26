#pragma once

#include "lbm/collision.hpp"
#include "lbm/grid.hpp"
#include <cstdint>
#include <vector>

namespace lbm {

// Zou-He boundary conditions on 2D D2Q9 boundary rows/columns (post-streaming)
void apply_zou_he_velocity_west(Grid2D& g, int y_min, int y_max, double ux, double uy);
void apply_zou_he_pressure_west(Grid2D& g, int y_min, int y_max, double rho_in, double uy = 0.0);
void apply_zou_he_pressure_east(Grid2D& g, int y_min, int y_max, double rho_out, double uy = 0.0);
void apply_convective_outflow_east(Grid2D& g, double u_conv);

// Lid-driven cavity stepper (half-way bounce-back on Left/Right/Bottom, moving top lid u_lid)
void step_cavity_2d(Grid2D& g, const CollisionParams& cp, double u_lid, int num_threads = 1);

// Curved cylinder geometry & Bouzidi / Staircase boundary links
struct CurvedLink {
    int fx = 0;
    int fy = 0;
    int dir = 0;     // Direction from fluid cell into solid cell (f_bar = opp[dir])
    double q = 0.5;  // Wall distance fraction along lattice link in (0, 1]
};

struct CylinderDomain {
    double cx = 0.0;
    double cy = 0.0;
    double radius = 1.0;
    bool use_bouzidi = true;
    std::vector<CurvedLink> links;
};

CylinderDomain build_cylinder_domain(Grid2D& g, double cx, double cy, double radius, bool use_bouzidi);

struct ForceResult {
    double fx = 0.0;
    double fy = 0.0;
    double cd = 0.0;
    double cl = 0.0;
};

// Advance cylinder cross-flow 1 step (uniform/parabolic inlet west, pressure outlet east, periodic or slip top/bottom)
// Returns deterministic momentum-exchange aerodynamic forces summed in fixed order.
ForceResult step_cylinder_2d(Grid2D& g, const CylinderDomain& cyl, const CollisionParams& cp,
                             double u_inlet, double rho_out = 1.0, int num_threads = 1);

// 3D D3Q19 square duct (periodic x, half-way bounce-back walls on y=0, y=ny-1, z=0, z=nz-1)
void step_duct_3d(Grid3D& g, const CollisionParams& cp, int num_threads = 1);

}  // namespace lbm
