#include "aad/svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace aad::svg {

namespace {

void write_header(std::ostringstream& os, int w, int h, const std::string& title) {
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
     << "\" viewBox=\"0 0 " << w << " " << h << "\">\n"
     << "  <rect width=\"" << w << "\" height=\"" << h << "\" fill=\"#ffffff\"/>\n"
     << "  <text x=\"" << (w / 2) << "\" y=\"28\" text-anchor=\"middle\" font-family=\"sans-serif\" "
     << "font-size=\"15\" font-weight=\"bold\" fill=\"#0f172a\">" << title << "</text>\n";
}

} // namespace

void write_aad_vs_analytic_svg(const std::string& path, const Gbm1StepResult& gbm,
                               double an_delta, double an_vega, double an_rho, double an_theta) {
  std::ostringstream os;
  write_header(os, 820, 420, "AAD Pathwise (200k) vs Analitik (z-score selisih)");
  const char* names[4] = {"Delta", "Vega", "Rho", "Theta"};
  const double z[4] = {
      std::abs(gbm.delta.mean - an_delta) / gbm.delta.se,
      std::abs(gbm.vega.mean - an_vega) / gbm.vega.se,
      std::abs(gbm.rho.mean - an_rho) / gbm.rho.se,
      std::abs(gbm.theta.mean - an_theta) / gbm.theta.se,
  };
  os << "  <rect x=\"80\" y=\"50\" width=\"700\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  const double y_3se = 360.0 - (3.0 / 3.5) * 310.0;
  os << "  <line x1=\"80\" y1=\"" << y_3se << "\" x2=\"780\" y2=\"" << y_3se
     << "\" stroke=\"#dc2626\" stroke-dasharray=\"6,4\"/>\n"
     << "  <text x=\"775\" y=\"" << (y_3se - 6)
     << "\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#dc2626\">Batas 3 SE</text>\n";
  for (int i = 0; i < 4; ++i) {
    const double bh = std::min(z[i], 3.5) / 3.5 * 310.0;
    const double x = 140.0 + i * 160.0;
    os << "  <rect x=\"" << x << "\" y=\"" << (360.0 - bh) << "\" width=\"70\" height=\"" << bh
       << "\" fill=\"#1d4ed8\"/>\n"
       << "  <text x=\"" << (x + 35) << "\" y=\"382\" text-anchor=\"middle\" font-family=\"sans-serif\" "
       << "font-size=\"12\" fill=\"#334155\">" << names[i] << "</text>\n"
       << "  <text x=\"" << (x + 35) << "\" y=\"" << (354.0 - bh)
       << "\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#0f172a\">"
       << std::fixed << std::setprecision(2) << z[i] << "</text>\n";
  }
  os << "</svg>\n";
  std::ofstream(path) << os.str();
}

void write_digital_bias_var_svg(const std::string& path, const DigitalResult& dig) {
  std::ostringstream os;
  write_header(os, 820, 420, "Digital Call Delta: Tradeoff Bias^2, Varians & MSE vs Epsilon");
  os << "  <rect x=\"80\" y=\"50\" width=\"700\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  double max_mse = 0.0;
  for (const auto& row : dig.smoothed) {
    max_mse = std::max(max_mse, row.mse);
  }
  max_mse = std::max(max_mse, 1e-12);
  std::ostringstream pts;
  for (std::size_t i = 0; i < dig.smoothed.size(); ++i) {
    const double x = 120.0 + (static_cast<double>(i) / 4.0) * 620.0;
    const double y = 360.0 - (dig.smoothed[i].mse / max_mse) * 290.0;
    pts << x << "," << y << " ";
    os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"4\" fill=\"#059669\"/>\n"
       << "  <text x=\"" << x << "\" y=\"382\" text-anchor=\"middle\" font-family=\"monospace\" "
       << "font-size=\"11\" fill=\"#334155\">eps=" << dig.smoothed[i].eps << "</text>\n";
  }
  os << "  <polyline fill=\"none\" stroke=\"#059669\" stroke-width=\"2.2\" points=\"" << pts.str()
     << "\"/>\n</svg>\n";
  std::ofstream(path) << os.str();
}

void write_basket_speedup_svg(const std::string& path, const BasketComparison& arith) {
  std::ostringstream os;
  write_header(os, 820, 420, "Keranjang 10 Aset (67 Sensitivitas): Waktu AAD vs Beda Hingga CRN");
  os << "  <rect x=\"80\" y=\"50\" width=\"700\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  const double max_t = std::max(arith.fd_time_sec, arith.aad_time_sec);
  const double h_aad = (arith.aad_time_sec / max_t) * 290.0;
  const double h_fd = (arith.fd_time_sec / max_t) * 290.0;
  os << "  <rect x=\"220\" y=\"" << (360.0 - h_aad) << "\" width=\"120\" height=\"" << h_aad
     << "\" fill=\"#059669\"/>\n"
     << "  <rect x=\"480\" y=\"" << (360.0 - h_fd) << "\" width=\"120\" height=\"" << h_fd
     << "\" fill=\"#dc2626\"/>\n"
     << "  <text x=\"280\" y=\"382\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">AAD ("
     << std::fixed << std::setprecision(2) << arith.aad_time_sec << " s)</text>\n"
     << "  <text x=\"540\" y=\"382\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\">FD CRN ("
     << std::fixed << std::setprecision(2) << arith.fd_time_sec << " s)</text>\n"
     << "  <text x=\"410\" y=\"85\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\" "
     << "font-weight=\"bold\" fill=\"#0f172a\">Speedup = " << std::setprecision(1) << arith.speedup_vs_fd
     << "x</text>\n</svg>\n";
  std::ofstream(path) << os.str();
}

void write_tape_growth_svg(const std::string& path, std::size_t gbm_nodes,
                           std::size_t basket_nodes, std::size_t asian_nodes) {
  std::ostringstream os;
  write_header(os, 820, 420, "Ukuran Tape Per Path (Node) per Model");
  os << "  <rect x=\"80\" y=\"50\" width=\"700\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  const double max_n = static_cast<double>(std::max({gbm_nodes, basket_nodes, asian_nodes}));
  const std::size_t vals[3] = {gbm_nodes, basket_nodes, asian_nodes};
  const char* labels[3] = {"GBM 1-step", "Basket 10-asset", "Asian 252-step"};
  for (int i = 0; i < 3; ++i) {
    const double bh = (static_cast<double>(vals[i]) / max_n) * 290.0;
    const double x = 160.0 + i * 200.0;
    os << "  <rect x=\"" << x << "\" y=\"" << (360.0 - bh) << "\" width=\"90\" height=\"" << bh
       << "\" fill=\"#6366f1\"/>\n"
       << "  <text x=\"" << (x + 45) << "\" y=\"382\" text-anchor=\"middle\" font-family=\"sans-serif\" "
       << "font-size=\"12\">" << labels[i] << "</text>\n"
       << "  <text x=\"" << (x + 45) << "\" y=\"" << (353.0 - bh)
       << "\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\">" << vals[i] << "</text>\n";
  }
  os << "</svg>\n";
  std::ofstream(path) << os.str();
}

} // namespace aad::svg
