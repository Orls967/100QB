#pragma once
#include "common.hpp"
#include <cmath>
#include <complex>

struct GL16 {
    double x[16], w[16];
    GL16() {
        constexpr double pi = 3.14159265358979323846;
        for (int i = 0; i < 8; ++i) {
            double z = std::cos(pi * (i + 0.75) / 16.5), pp = 0.0;
            for (int it = 0; it < 30; ++it) {
                double p1 = 1.0, p2 = 0.0;
                for (int k = 1; k <= 16; ++k) {
                    double p3 = p2; p2 = p1;
                    p1 = ((2.0 * k - 1.0) * z * p2 - (k - 1.0) * p3) / k;
                }
                pp = 16.0 * (z * p1 - p2) / (z * z - 1.0);
                double dz = p1 / pp; z -= dz;
                if (std::abs(dz) < 1e-16) break;
            }
            x[i] = -z; x[15 - i] = z;
            w[i] = w[15 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
        }
    }
};

struct HestonPsi {
    Params p; double S0, T, v0;
    std::complex<double> operator()(std::complex<double> z) const {
        using cd = std::complex<double>;
        cd iz(-z.imag(), z.real()), b = p.kappa - p.rho * p.sigma * iz;
        cd d = std::sqrt(b * b + p.sigma * p.sigma * (iz + z * z)), g = (b - d) / (b + d), edt = std::exp(-d * T);
        double s2 = p.sigma * p.sigma;
        cd C = iz * (std::log(S0) + p.r * T) + (p.kappa * p.theta / s2) * ((b - d) * T - 2.0 * std::log((1.0 - g * edt) / (1.0 - g)));
        return std::exp(C + ((b - d) / s2) * ((1.0 - edt) / (1.0 - g * edt)) * v0);
    }
};

struct BsPsi {
    double S0, r, T, sig;
    std::complex<double> operator()(std::complex<double> z) const {
        std::complex<double> iz(-z.imag(), z.real());
        return std::exp(iz * (std::log(S0) + (r - 0.5 * sig * sig) * T) - 0.5 * sig * sig * z * z * T);
    }
};

template <class Psi>
double gp_price(const Psi& psi, double S0, double K, double r, double T, bool is_call, const GL16& gl = GL16()) {
    using cd = std::complex<double>;
    constexpr double pi = 3.14159265358979323846;
    double lk = std::log(K), df = std::exp(-r * T), sum = 0.0;
    cd phi_mi = psi(cd(0.0, -1.0));
    for (int panel = 0; panel < 200; ++panel) {
        double mid = panel + 0.5;
        for (int q = 0; q < 16; ++q) {
            double u = mid + 0.5 * gl.x[q];
            cd e = std::exp(cd(0.0, -u * lk));
            double p1 = (e * psi(cd(u, -1.0)) / (cd(0.0, u) * phi_mi)).real();
            double p2 = (e * psi(cd(u, 0.0)) / cd(0.0, u)).real();
            sum += 0.5 * gl.w[q] * (S0 * p1 - K * df * p2);
        }
    }
    return 0.5 * (is_call ? (S0 - K * df) : (K * df - S0)) + sum / pi;
}
template <class Psi> double gp_call(const Psi& psi, double S0, double K, double r, double T, const GL16& gl = GL16()) {
    return gp_price(psi, S0, K, r, T, true, gl);
}
template <class Psi> double gp_put(const Psi& psi, double S0, double K, double r, double T, const GL16& gl = GL16()) {
    return gp_price(psi, S0, K, r, T, false, gl);
}

template <class Psi>
double lewis_call(const Psi& psi, double S0, double K, double r, double T, const GL16& gl = GL16()) {
    using cd = std::complex<double>;
    constexpr double pi = 3.14159265358979323846;
    double lk = std::log(K), sum = 0.0;
    for (int panel = 0; panel < 200; ++panel) {
        double mid = panel + 0.5;
        for (int q = 0; q < 16; ++q) {
            double u = mid + 0.5 * gl.x[q];
            sum += 0.5 * gl.w[q] * (std::exp(cd(0.0, -u * lk)) * psi(cd(u, -0.5))).real() / (u * u + 0.25);
        }
    }
    return S0 - (std::sqrt(K) * std::exp(-r * T) / pi) * sum;
}

inline double bs_call(double S0, double K, double r, double T, double sig) {
    double srt = sig * std::sqrt(T), d1 = (std::log(S0 / K) + (r + 0.5 * sig * sig) * T) / srt;
    auto ncdf = [](double z) { return 0.5 * std::erfc(-z * 0.70710678118654752440); };
    return S0 * ncdf(d1) - K * std::exp(-r * T) * ncdf(d1 - srt);
}
