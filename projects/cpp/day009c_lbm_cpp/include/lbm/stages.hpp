#pragma once

#include "lbm/io.hpp"
#include <cstdint>

namespace lbm {

bool run_selftest();

void run_group_a(RunContext& ctx);  // S01-S08: Foundation
void run_group_b(RunContext& ctx);  // S09-S16: Collision, shear wave, stability, Galilean
void run_group_c(RunContext& ctx);  // S17-S28: Taylor-Green, Ma scaling, Poiseuille, Couette, Zou-He, Acoustic
void run_group_d(RunContext& ctx);  // S29-S40: Cavity & Cylinder benchmarks
void run_group_e(RunContext& ctx);  // S41-S46: 3D D3Q19 shear wave, square duct, TG3D
void run_group_f(RunContext& ctx);  // S47-S54: Schemes (Push/Pull/AA), Layouts (AoS/SoA/Tile), Thread parity
void run_group_g(RunContext& ctx);  // S55-S60: Optimization ladder V0-V6 + Bandwidth/Roofline + Grid scaling
void run_group_h(RunContext& ctx);  // S61-S65: Determinism, Delta+RLE, Gallery SVGs, Manifest & results.json

// Helper returning bitwise hash of 2D/3D test lattice evolved deterministically with n_threads
std::uint64_t benchmark_thread_parity_hash(int n_threads);

}  // namespace lbm
