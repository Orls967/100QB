#pragma once

#include "lbm/grid.hpp"
#include <string>
#include <utility>
#include <vector>

namespace lbm {

struct OlsResult {
    double slope = 0.0;
    double slope_se = 0.0;
    double intercept = 0.0;
    double r2 = 0.0;
};

struct MedianMad {
    double median = 0.0;
    double mad = 0.0;
};

struct VortexCenter {
    double x = 0.0;
    double y = 0.0;
    double psi_ext = 0.0;
    int corner_eddies = 0;
    double u_center_min = 0.0;
    double v_center_min = 0.0;
    double v_center_max = 0.0;
};

struct StrouhalStats {
    double freq_zc = 0.0;
    double freq_dft = 0.0;
    double strouhal = 0.0;
    double cl_amp = 0.0;
    int zero_crossings = 0;
};

OlsResult ols_fit(const std::vector<double>& xs, const std::vector<double>& ys);
MedianMad median_mad(std::vector<double> xs);

// Cavity primary vortex via integrated streamfunction + sub-grid quadratic peak + corner vortices check
VortexCenter analyze_cavity_vortices(const Grid2D& g);

// Wake recirculation bubble length Lr / D behind cylinder (normalized by diameter D)
double measure_recirculation_length(const Grid2D& g, double cx, double cy, double diameter);

// Frequency & Strouhal analysis on scalar time series
StrouhalStats analyze_strouhal(const std::vector<double>& ts, double diameter, double u_inf);

// Exact 3D Poiseuille square duct analytical Fourier series velocity u_x(y_phys, z_phys) on [0, H] x [0, H]
double exact_square_duct_ux(double y_pos, double z_pos, double height, double fx, double rho0, double nu);

// Measure reference host memory copy bandwidth (GB/s) for Roofline normalization
double measure_copy_bandwidth_gbps(std::size_t num_doubles, int repeats);

// Lossless quantized Delta + RLE compression on 1D field vector
struct RleMetrics {
    std::size_t raw_bytes = 0;
    std::size_t compressed_bytes = 0;
    double ratio = 1.0;
    std::int64_t max_roundtrip_q_err = 0;
};

RleMetrics compress_delta_rle(const std::vector<double>& field, double quant_step);

}  // namespace lbm
