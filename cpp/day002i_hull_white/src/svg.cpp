#include "hw/svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hw::svg {

namespace {

std::string fmt(double v, int prec = 2) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(prec) << v;
    return oss.str();
}

void write_header(std::ofstream& out, int w, int h, const std::string& title) {
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
        << "\" viewBox=\"0 0 " << w << " " << h << "\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
        << "<style>text{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;fill:#1f2937;}"
        << ".title{font-size:16px;font-weight:700;}.axis{stroke:#4b5563;stroke-width:1.2;}"
        << ".grid{stroke:#e5e7eb;stroke-width:0.8;}.lbl{font-size:11px;fill:#4b5563;}"
        << ".leg{font-size:12px;}</style>\n"
        << "<text x=\"" << w / 2 << "\" y=\"28\" text-anchor=\"middle\" class=\"title\">"
        << title << "</text>\n";
}

}  // namespace

void write_tree_convergence_svg(const std::string& path,
                                const std::vector<int>& m_values,
                                const std::vector<double>& zbp_errors,
                                const std::vector<double>& zbc_errors,
                                const std::vector<double>& cap_errors,
                                double slope_zbp,
                                double slope_zbc) {
    std::ofstream out(path);
    const int w = 780;
    const int h = 440;
    write_header(out, w, h, "Hull-White Trinomial Tree Convergence (Log-Log Error vs Steps/Year)");

    const double px0 = 85.0;
    const double px1 = 730.0;
    const double py0 = 60.0;
    const double py1 = 375.0;

    double xmin = std::log10(static_cast<double>(m_values.front()));
    double xmax = std::log10(static_cast<double>(m_values.back()));
    double ymin = -6.5;
    double ymax = -2.5;

    for (std::size_t i = 0; i < m_values.size(); ++i) {
        for (double err : {zbp_errors[i], zbc_errors[i], cap_errors[i]}) {
            double ly = std::log10(std::max(1e-12, err));
            ymin = std::min(ymin, std::floor(ly * 2.0) / 2.0 - 0.3);
            ymax = std::max(ymax, std::ceil(ly * 2.0) / 2.0 + 0.3);
        }
    }

    auto map_x = [&](double lx) { return px0 + (lx - xmin) / (xmax - xmin) * (px1 - px0); };
    auto map_y = [&](double ly) { return py1 - (ly - ymin) / (ymax - ymin) * (py1 - py0); };

    out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
    out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

    for (int m : m_values) {
        double x = map_x(std::log10(static_cast<double>(m)));
        out << "<line x1=\"" << x << "\" y1=\"" << py0 << "\" x2=\"" << x << "\" y2=\"" << py1 << "\" class=\"grid\"/>\n";
        out << "<text x=\"" << x << "\" y=\"" << py1 + 18 << "\" text-anchor=\"middle\" class=\"lbl\">m=" << m << "</text>\n";
    }
    for (int p = static_cast<int>(std::ceil(ymin)); p <= static_cast<int>(std::floor(ymax)); ++p) {
        double y = map_y(static_cast<double>(p));
        out << "<line x1=\"" << px0 << "\" y1=\"" << y << "\" x2=\"" << px1 << "\" y2=\"" << y << "\" class=\"grid\"/>\n";
        out << "<text x=\"" << px0 - 10 << "\" y=\"" << y + 4 << "\" text-anchor=\"end\" class=\"lbl\">1e" << p << "</text>\n";
    }

    auto draw_series = [&](const std::vector<double>& errs, const std::string& color) {
        out << "<polyline fill=\"none\" stroke=\"" << color << "\" stroke-width=\"2.2\" points=\"";
        for (std::size_t i = 0; i < m_values.size(); ++i) {
            double x = map_x(std::log10(static_cast<double>(m_values[i])));
            double y = map_y(std::log10(std::max(1e-12, errs[i])));
            out << fmt(x) << "," << fmt(y) << " ";
        }
        out << "\"/>\n";
        for (std::size_t i = 0; i < m_values.size(); ++i) {
            double x = map_x(std::log10(static_cast<double>(m_values[i])));
            double y = map_y(std::log10(std::max(1e-12, errs[i])));
            out << "<circle cx=\"" << fmt(x) << "\" cy=\"" << fmt(y) << "\" r=\"4\" fill=\"" << color << "\"/>\n";
        }
    };

    draw_series(zbp_errors, "#2563eb");
    draw_series(zbc_errors, "#dc2626");
    draw_series(cap_errors, "#059669");

    out << "<text x=\"" << (px0 + px1) / 2 << "\" y=\"418\" text-anchor=\"middle\" class=\"lbl\">Steps per year (m, log scale)</text>\n";
    out << "<rect x=\"480\" y=\"65\" width=\"240\" height=\"76\" fill=\"#ffffff\" stroke=\"#d1d5db\" rx=\"4\"/>\n";
    out << "<line x1=\"492\" y1=\"84\" x2=\"518\" y2=\"84\" stroke=\"#2563eb\" stroke-width=\"2.5\"/>"
        << "<text x=\"526\" y=\"88\" class=\"leg\">ZBP (slope " << fmt(slope_zbp, 2) << ")</text>\n";
    out << "<line x1=\"492\" y1=\"106\" x2=\"518\" y2=\"106\" stroke=\"#dc2626\" stroke-width=\"2.5\"/>"
        << "<text x=\"526\" y=\"110\" class=\"leg\">ZBC (slope " << fmt(slope_zbc, 2) << ")</text>\n";
    out << "<line x1=\"492\" y1=\"128\" x2=\"518\" y2=\"128\" stroke=\"#059669\" stroke-width=\"2.5\"/>"
        << "<text x=\"526\" y=\"132\" class=\"leg\">5y Cap</text>\n";
    out << "</svg>\n";
}

void write_curve_and_swaptions_svg(const std::string& path,
                                   const std::vector<double>& maturities,
                                   const std::vector<double>& zero_rates,
                                   const std::vector<double>& fwd_rates,
                                   const std::vector<double>& expiries,
                                   const std::vector<double>& jamshidian_prices,
                                   const std::vector<double>& tree_prices,
                                   const std::vector<double>& mc_prices) {
    std::ofstream out(path);
    const int w = 880;
    const int h = 420;
    write_header(out, w, h, "Nelson-Siegel Term Structure &amp; 10y Coterminal Payer Swaptions");

    // Left panel: Curve
    {
        const double px0 = 65.0, px1 = 410.0, py0 = 65.0, py1 = 360.0;
        double ymin = 0.015, ymax = 0.065;
        auto mx = [&](double t) { return px0 + (t / 10.0) * (px1 - px0); };
        auto my = [&](double r) { return py1 - ((r - ymin) / (ymax - ymin)) * (py1 - py0); };

        out << "<text x=\"237\" y=\"52\" text-anchor=\"middle\" class=\"leg\">Initial Curve (Zero &amp; Inst. Forward)</text>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

        for (double r = 0.02; r <= 0.061; r += 0.01) {
            double y = my(r);
            out << "<line x1=\"" << px0 << "\" y1=\"" << y << "\" x2=\"" << px1 << "\" y2=\"" << y << "\" class=\"grid\"/>\n";
            out << "<text x=\"" << px0 - 8 << "\" y=\"" << y + 4 << "\" text-anchor=\"end\" class=\"lbl\">" << fmt(r * 100.0, 1) << "%</text>\n";
        }

        auto draw_c = [&](const std::vector<double>& ys, const std::string& col) {
            out << "<polyline fill=\"none\" stroke=\"" << col << "\" stroke-width=\"2.2\" points=\"";
            for (std::size_t i = 0; i < maturities.size(); ++i) {
                out << fmt(mx(maturities[i])) << "," << fmt(my(ys[i])) << " ";
            }
            out << "\"/>\n";
        };
        draw_c(zero_rates, "#2563eb");
        draw_c(fwd_rates, "#d97706");
        out << "<text x=\"237\" y=\"392\" text-anchor=\"middle\" class=\"lbl\">Maturity (years)</text>\n";
        out << "<text x=\"280\" y=\"92\" class=\"lbl\" fill=\"#2563eb\">Blue: z(t) | Orange: f(0,t)</text>\n";
    }

    // Right panel: Swaptions
    {
        const double px0 = 490.0, px1 = 845.0, py0 = 65.0, py1 = 360.0;
        double ymax = 0.04;
        for (double v : jamshidian_prices) {
            ymax = std::max(ymax, v * 1.15);
        }
        auto mx = [&](double t) { return px0 + ((t - 0.5) / 9.0) * (px1 - px0); };
        auto my = [&](double p) { return py1 - (p / ymax) * (py1 - py0); };

        out << "<text x=\"667\" y=\"52\" text-anchor=\"middle\" class=\"leg\">Coterminal Swaptions (Jamshidian / Tree / MC)</text>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

        for (std::size_t i = 0; i < expiries.size(); ++i) {
            double x = mx(expiries[i]);
            out << "<text x=\"" << x << "\" y=\"" << py1 + 16 << "\" text-anchor=\"middle\" class=\"lbl\">" << fmt(expiries[i], 0) << "y</text>\n";
            out << "<circle cx=\"" << fmt(x) << "\" cy=\"" << fmt(my(jamshidian_prices[i])) << "\" r=\"5\" fill=\"#111827\"/>\n";
            out << "<circle cx=\"" << fmt(x) << "\" cy=\"" << fmt(my(tree_prices[i])) << "\" r=\"3.5\" fill=\"#2563eb\"/>\n";
            out << "<rect x=\"" << fmt(x - 2.5) << "\" y=\"" << fmt(my(mc_prices[i]) - 2.5) << "\" width=\"5\" height=\"5\" fill=\"#10b981\"/>\n";
        }
        out << "<text x=\"667\" y=\"392\" text-anchor=\"middle\" class=\"lbl\">Option Expiry T0 (10y coterminal, K=4.5%)</text>\n";
    }
    out << "</svg>\n";
}

void write_calib_scatter_svg(const std::string& path,
                             double true_a,
                             double true_sigma,
                             const std::vector<double>& a_samples,
                             const std::vector<double>& sigma_samples,
                             double corr_a_sigma,
                             double cond_number) {
    std::ofstream out(path);
    const int w = 740;
    const int h = 440;
    write_header(out, w, h, "Noisy Calibration Scatter (R=100, 0.5% Rel Price Noise)");

    double amin = true_a * 0.7, amax = true_a * 1.3;
    double smin = true_sigma * 0.7, smax = true_sigma * 1.3;
    for (double a : a_samples) {
        amin = std::min(amin, a);
        amax = std::max(amax, a);
    }
    for (double s : sigma_samples) {
        smin = std::min(smin, s);
        smax = std::max(smax, s);
    }
    double apad = 0.08 * (amax - amin);
    double spad = 0.08 * (smax - smin);
    amin -= apad;
    amax += apad;
    smin -= spad;
    smax += spad;

    const double px0 = 90.0, px1 = 690.0, py0 = 65.0, py1 = 375.0;
    auto mx = [&](double a) { return px0 + (a - amin) / (amax - amin) * (px1 - px0); };
    auto my = [&](double s) { return py1 - (s - smin) / (smax - smin) * (py1 - py0); };

    out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
    out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

    for (std::size_t i = 0; i < a_samples.size(); ++i) {
        out << "<circle cx=\"" << fmt(mx(a_samples[i])) << "\" cy=\"" << fmt(my(sigma_samples[i]))
            << "\" r=\"3.2\" fill=\"#2563eb\" fill-opacity=\"0.65\"/>\n";
    }
    double tx = mx(true_a);
    double ty = my(true_sigma);
    out << "<line x1=\"" << fmt(tx - 8) << "\" y1=\"" << fmt(ty) << "\" x2=\"" << fmt(tx + 8) << "\" y2=\"" << fmt(ty) << "\" stroke=\"#dc2626\" stroke-width=\"2.5\"/>\n";
    out << "<line x1=\"" << fmt(tx) << "\" y1=\"" << fmt(ty - 8) << "\" x2=\"" << fmt(tx) << "\" y2=\"" << fmt(ty + 8) << "\" stroke=\"#dc2626\" stroke-width=\"2.5\"/>\n";

    out << "<text x=\"390\" y=\"415\" text-anchor=\"middle\" class=\"lbl\">Mean reversion a (true = " << fmt(true_a, 3)
        << ") | corr(a, sigma) = " << fmt(corr_a_sigma, 3) << " | cond(J) = " << fmt(cond_number, 1) << "</text>\n";
    out << "<text x=\"25\" y=\"220\" text-anchor=\"middle\" transform=\"rotate(-90 25 220)\" class=\"lbl\">Volatility sigma</text>\n";
    out << "</svg>\n";
}

void write_bermudan_bounds_dv01_svg(const std::string& path,
                                    const std::vector<int>& m_values,
                                    const std::vector<double>& bermudan_tree,
                                    double bermudan_rich,
                                    double lsm_price,
                                    double lsm_se,
                                    double max_euro,
                                    double sum_euro,
                                    double dv01_euro,
                                    double dv01_berm) {
    std::ofstream out(path);
    const int w = 860;
    const int h = 420;
    write_header(out, w, h, "10y Bermudan Swaption: Tree vs LSM &amp; Parallel +1bp DV01");

    // Left: Bermudan convergence + LSM
    {
        const double px0 = 75.0, px1 = 470.0, py0 = 65.0, py1 = 360.0;
        double ymin = std::min(max_euro, lsm_price - 3.0 * lsm_se) * 0.96;
        double ymax = std::max(bermudan_rich, lsm_price + 3.0 * lsm_se) * 1.04;
        auto mx = [&](std::size_t idx) {
            return px0 + (static_cast<double>(idx) / static_cast<double>(m_values.size() - 1)) * (px1 - px0);
        };
        auto my = [&](double v) { return py1 - ((v - ymin) / (ymax - ymin)) * (py1 - py0); };

        out << "<text x=\"272\" y=\"52\" text-anchor=\"middle\" class=\"leg\">Bermudan Price (Tree vs LSM, max Euro = "
            << fmt(max_euro, 4) << ")</text>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

        out << "<line x1=\"" << px0 << "\" y1=\"" << fmt(my(bermudan_rich)) << "\" x2=\"" << px1 << "\" y2=\"" << fmt(my(bermudan_rich))
            << "\" stroke=\"#059669\" stroke-dasharray=\"5,4\" stroke-width=\"1.8\"/>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << fmt(my(lsm_price)) << "\" x2=\"" << px1 << "\" y2=\"" << fmt(my(lsm_price))
            << "\" stroke=\"#dc2626\" stroke-dasharray=\"3,3\" stroke-width=\"1.8\"/>\n";

        out << "<polyline fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.2\" points=\"";
        for (std::size_t i = 0; i < m_values.size(); ++i) {
            out << fmt(mx(i)) << "," << fmt(my(bermudan_tree[i])) << " ";
            out << "";
        }
        out << "\"/>\n";
        for (std::size_t i = 0; i < m_values.size(); ++i) {
            out << "<circle cx=\"" << fmt(mx(i)) << "\" cy=\"" << fmt(my(bermudan_tree[i])) << "\" r=\"4\" fill=\"#2563eb\"/>\n";
            out << "<text x=\"" << fmt(mx(i)) << "\" y=\"" << py1 + 16 << "\" text-anchor=\"middle\" class=\"lbl\">m=" << m_values[i] << "</text>\n";
        }
        out << "<text x=\"272\" y=\"392\" text-anchor=\"middle\" class=\"lbl\">Blue: Tree | Green: Richardson | Red: LSM (" << fmt(lsm_price, 4) << ") | SumEuro=" << fmt(sum_euro, 3) << "</text>\n";
    }

    // Right: DV01 comparison bar chart
    {
        const double px0 = 540.0, px1 = 810.0, py0 = 65.0, py1 = 360.0;
        double ymax = std::max(std::abs(dv01_euro), std::abs(dv01_berm)) * 1.25;
        auto my = [&](double v) { return py1 - (std::abs(v) / ymax) * (py1 - py0); };

        out << "<text x=\"675\" y=\"52\" text-anchor=\"middle\" class=\"leg\">Parallel DV01 (+-1 bp)</text>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py1 << "\" x2=\"" << px1 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";
        out << "<line x1=\"" << px0 << "\" y1=\"" << py0 << "\" x2=\"" << px0 << "\" y2=\"" << py1 << "\" class=\"axis\"/>\n";

        double y_e = my(dv01_euro);
        double y_b = my(dv01_berm);
        out << "<rect x=\"590\" y=\"" << fmt(y_e) << "\" width=\"65\" height=\"" << fmt(py1 - y_e) << "\" fill=\"#2563eb\"/>\n";
        out << "<text x=\"622\" y=\"" << fmt(y_e - 8) << "\" text-anchor=\"middle\" class=\"lbl\">" << fmt(dv01_euro * 1e4, 2) << " bp</text>\n";
        out << "<text x=\"622\" y=\"" << py1 + 16 << "\" text-anchor=\"middle\" class=\"lbl\">5x5 Euro</text>\n";

        out << "<rect x=\"705\" y=\"" << fmt(y_b) << "\" width=\"65\" height=\"" << fmt(py1 - y_b) << "\" fill=\"#059669\"/>\n";
        out << "<text x=\"737\" y=\"" << fmt(y_b - 8) << "\" text-anchor=\"middle\" class=\"lbl\">" << fmt(dv01_berm * 1e4, 2) << " bp</text>\n";
        out << "<text x=\"737\" y=\"" << py1 + 16 << "\" text-anchor=\"middle\" class=\"lbl\">10y Berm</text>\n";
    }
    out << "</svg>\n";
}

}  // namespace hw::svg
