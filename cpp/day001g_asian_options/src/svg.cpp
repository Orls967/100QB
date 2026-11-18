#include "asian/svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace asian {

namespace {

std::string fmt_fixed(double v, int prec = 1) {
  std::ostringstream os;
  os << std::fixed << std::setprecision(prec) << v;
  return os.str();
}

std::string fmt_sci(double v) {
  std::ostringstream os;
  os << std::scientific << std::setprecision(0) << v;
  return os.str();
}

} // namespace

std::string render_convergence_svg(std::span<const ConvergencePoint> points, const Slopes &slopes) {
  constexpr int w = 900;
  constexpr int h = 520;
  constexpr int pad_l = 85;
  constexpr int pad_r = 30;
  constexpr int pad_t = 55;
  constexpr int pad_b = 65;
  constexpr int plot_w = w - pad_l - pad_r;
  constexpr int plot_h = h - pad_t - pad_b;

  const double min_log_x = std::log10(900.0);
  const double max_log_x = std::log10(230000.0);
  double min_se = std::numeric_limits<double>::infinity();
  double max_se = 0.0;
  for (const auto &pt : points) {
    min_se = std::min(min_se, pt.se_control_variate);
    max_se = std::max(max_se, pt.se_naive);
  }
  const double min_log_y = std::log10(min_se * 0.7);
  const double max_log_y = std::log10(max_se * 1.35);

  auto map_x = [&](std::size_t n) {
    return pad_l + ((std::log10(static_cast<double>(n)) - min_log_x) / (max_log_x - min_log_x)) *
                       plot_w;
  };
  auto map_y = [&](double se) {
    return pad_t + (1.0 - (std::log10(se) - min_log_y) / (max_log_y - min_log_y)) * plot_h;
  };

  std::ostringstream os;
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
     << "\" viewBox=\"0 0 " << w << " " << h
     << "\" font-family=\"system-ui, -apple-system, sans-serif\">\n";
  os << "<rect width=\"" << w << "\" height=\"" << h << "\" fill=\"#ffffff\"/>\n";
  os << "<text x=\"" << (w / 2)
     << "\" y=\"30\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\" fill=\"#111827\">"
     << "Konvergensi Standard Error vs Jumlah Path (Log-Log, C++20)</text>\n";

  for (const std::size_t nt : {1000ULL, 10000ULL, 100000ULL}) {
    const double x = map_x(nt);
    os << "<line x1=\"" << fmt_fixed(x) << "\" y1=\"" << pad_t << "\" x2=\"" << fmt_fixed(x)
       << "\" y2=\"" << (pad_t + plot_h)
       << "\" stroke=\"#e5e7eb\" stroke-dasharray=\"4,4\"/>\n";
    os << "<text x=\"" << fmt_fixed(x) << "\" y=\"" << (pad_t + plot_h + 22)
       << "\" text-anchor=\"middle\" font-size=\"12\" fill=\"#374151\">" << nt << "</text>\n";
  }

  for (int exp = static_cast<int>(std::ceil(min_log_y));
       exp <= static_cast<int>(std::floor(max_log_y)); ++exp) {
    const double val = std::pow(10.0, exp);
    const double y = map_y(val);
    os << "<line x1=\"" << pad_l << "\" y1=\"" << fmt_fixed(y) << "\" x2=\"" << (pad_l + plot_w)
       << "\" y2=\"" << fmt_fixed(y) << "\" stroke=\"#e5e7eb\" stroke-dasharray=\"4,4\"/>\n";
    os << "<text x=\"" << (pad_l - 10) << "\" y=\"" << fmt_fixed(y + 4)
       << "\" text-anchor=\"end\" font-size=\"12\" fill=\"#374151\">" << fmt_sci(val)
       << "</text>\n";
  }

  os << "<rect x=\"" << pad_l << "\" y=\"" << pad_t << "\" width=\"" << plot_w << "\" height=\""
     << plot_h << "\" fill=\"none\" stroke=\"#9ca3af\"/>\n";

  struct SeriesInfo {
    std::string label;
    const char *color;
    double (*getter)(const ConvergencePoint &);
  };
  const std::array<SeriesInfo, 3> series{{
      {"Naive MC (slope=" + fmt_fixed(slopes.naive, 3) + ")", "#1f6feb",
       [](const ConvergencePoint &p) { return p.se_naive; }},
      {"Antithetic (slope=" + fmt_fixed(slopes.antithetic, 3) + ")", "#d1495b",
       [](const ConvergencePoint &p) { return p.se_antithetic; }},
      {"Control Variate (slope=" + fmt_fixed(slopes.control_variate, 3) + ")", "#2a9d8f",
       [](const ConvergencePoint &p) { return p.se_control_variate; }},
  }};

  for (const auto &s : series) {
    os << "<polyline fill=\"none\" stroke=\"" << s.color << "\" stroke-width=\"2.5\" points=\"";
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (i > 0) os << ' ';
      os << fmt_fixed(map_x(points[i].n)) << ',' << fmt_fixed(map_y(s.getter(points[i])));
    }
    os << "\"/>\n";
    for (const auto &pt : points) {
      os << "<circle cx=\"" << fmt_fixed(map_x(pt.n)) << "\" cy=\""
         << fmt_fixed(map_y(s.getter(pt))) << "\" r=\"4\" fill=\"" << s.color << "\"/>\n";
    }
  }

  os << "<rect x=\"" << (pad_l + plot_w - 252) << "\" y=\"" << (pad_t + 8)
     << "\" width=\"244\" height=\"72\" rx=\"6\" fill=\"#ffffff\" fill-opacity=\"0.9\" stroke=\"#d1d5db\"/>\n";
  for (std::size_t i = 0; i < series.size(); ++i) {
    const int ly = pad_t + 22 + static_cast<int>(i) * 22;
    const int lx = pad_l + plot_w - 240;
    os << "<line x1=\"" << lx << "\" y1=\"" << ly << "\" x2=\"" << (lx + 24) << "\" y2=\"" << ly
       << "\" stroke=\"" << series[i].color << "\" stroke-width=\"2.5\"/>\n";
    os << "<text x=\"" << (lx + 32) << "\" y=\"" << (ly + 4)
       << "\" font-size=\"12\" fill=\"#111827\">" << series[i].label << "</text>\n";
  }

  os << "<text x=\"" << (pad_l + plot_w / 2) << "\" y=\"" << (h - 16)
     << "\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\">Jumlah Path N (skala log)</text>\n";
  os << "<text x=\"22\" y=\"" << (pad_t + plot_h / 2)
     << "\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\" transform=\"rotate(-90 22 "
     << (pad_t + plot_h / 2) << ")\">Standard Error (SE)</text>\n";
  os << "</svg>\n";
  return os.str();
}

std::string render_estimators_svg(std::span<const EstimatorSummary> estimators) {
  constexpr int w = 900;
  constexpr int h = 440;
  constexpr int pad_l = 210;
  constexpr int pad_r = 60;
  constexpr int pad_t = 60;
  constexpr int pad_b = 60;
  constexpr int plot_w = w - pad_l - pad_r;
  constexpr int plot_h = h - pad_t - pad_b;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  for (const auto &e : estimators) {
    min_x = std::min(min_x, e.ci_low);
    max_x = std::max(max_x, e.ci_high);
  }
  const double margin = (max_x - min_x) * 0.15;
  min_x -= margin;
  max_x += margin;

  auto map_x = [&](double val) { return pad_l + ((val - min_x) / (max_x - min_x)) * plot_w; };
  constexpr std::array<const char *, 3> colors{"#1f6feb", "#d1495b", "#2a9d8f"};

  std::ostringstream os;
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
     << "\" viewBox=\"0 0 " << w << " " << h
     << "\" font-family=\"system-ui, -apple-system, sans-serif\">\n";
  os << "<rect width=\"" << w << "\" height=\"" << h << "\" fill=\"#ffffff\"/>\n";
  os << "<text x=\"" << (w / 2)
     << "\" y=\"32\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"bold\" fill=\"#111827\">"
     << "Perbandingan Estimasi Harga + CI 95% &amp; VRF (200k Path, C++20)</text>\n";
  os << "<rect x=\"" << pad_l << "\" y=\"" << pad_t << "\" width=\"" << plot_w << "\" height=\""
     << plot_h << "\" fill=\"none\" stroke=\"#9ca3af\"/>\n";

  for (std::size_t idx = 0; idx < estimators.size(); ++idx) {
    const auto &e = estimators[idx];
    const double cy = pad_t + ((static_cast<double>(idx) + 0.5) /
                               static_cast<double>(estimators.size())) *
                                  plot_h;
    const double x_lo = map_x(e.ci_low);
    const double x_hi = map_x(e.ci_high);
    const double x_mid = map_x(e.price);
    const char *col = colors[idx % colors.size()];

    os << "<text x=\"" << (pad_l - 12) << "\" y=\"" << fmt_fixed(cy + 4)
       << "\" text-anchor=\"end\" font-size=\"12\" fill=\"#111827\">" << e.name << "</text>\n";
    os << "<line x1=\"" << fmt_fixed(x_lo) << "\" y1=\"" << fmt_fixed(cy) << "\" x2=\""
       << fmt_fixed(x_hi) << "\" y2=\"" << fmt_fixed(cy) << "\" stroke=\"" << col
       << "\" stroke-width=\"3\"/>\n";
    os << "<line x1=\"" << fmt_fixed(x_lo) << "\" y1=\"" << fmt_fixed(cy - 6) << "\" x2=\""
       << fmt_fixed(x_lo) << "\" y2=\"" << fmt_fixed(cy + 6) << "\" stroke=\"" << col
       << "\" stroke-width=\"2\"/>\n";
    os << "<line x1=\"" << fmt_fixed(x_hi) << "\" y1=\"" << fmt_fixed(cy - 6) << "\" x2=\""
       << fmt_fixed(x_hi) << "\" y2=\"" << fmt_fixed(cy + 6) << "\" stroke=\"" << col
       << "\" stroke-width=\"2\"/>\n";
    os << "<circle cx=\"" << fmt_fixed(x_mid) << "\" cy=\"" << fmt_fixed(cy) << "\" r=\"5\" fill=\""
       << col << "\"/>\n";
    os << "<text x=\"" << fmt_fixed(x_hi + 8) << "\" y=\"" << fmt_fixed(cy + 4)
       << "\" font-size=\"11\" fill=\"#374151\">" << fmt_fixed(e.price, 4) << " (VRF "
       << fmt_fixed(e.vrf, 1) << "x)</text>\n";
  }

  os << "<text x=\"" << (pad_l + plot_w / 2) << "\" y=\"" << (h - 18)
     << "\" text-anchor=\"middle\" font-size=\"13\" fill=\"#111827\">Harga Opsi Asia</text>\n";
  os << "</svg>\n";
  return os.str();
}

} // namespace asian
