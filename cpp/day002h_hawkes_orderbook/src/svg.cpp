#include "hawkes_ob/svg.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hob::svg {

namespace {

void write_header(std::ostringstream& os, int w, int h, const std::string& title) {
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
     << "\" viewBox=\"0 0 " << w << " " << h << "\">\n"
     << "  <rect width=\"" << w << "\" height=\"" << h << "\" fill=\"#ffffff\"/>\n"
     << "  <text x=\"" << (w / 2) << "\" y=\"28\" text-anchor=\"middle\" font-family=\"sans-serif\" "
     << "font-size=\"15\" font-weight=\"bold\" fill=\"#0f172a\">" << title << "</text>\n";
}

} // namespace

void write_fano_and_qq_svg(const std::string& path, const std::vector<double>& win_sizes,
                           const std::vector<double>& fanos, const std::vector<double>& res_exp1) {
  std::ostringstream os;
  write_header(os, 900, 420, "Hawkes: Fano Factor vs Ukuran Jendela & QQ Residu Exp(1)");
  os << "  <rect x=\"60\" y=\"50\" width=\"360\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n"
     << "  <rect x=\"490\" y=\"50\" width=\"360\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";

  const double y4 = 360.0 - (4.0 / 5.0) * 310.0;
  os << "  <line x1=\"60\" y1=\"" << y4 << "\" x2=\"420\" y2=\"" << y4
     << "\" stroke=\"#dc2626\" stroke-dasharray=\"5,4\"/>\n";
  std::ostringstream pts1;
  for (std::size_t i = 0; i < win_sizes.size(); ++i) {
    const double x = 80.0 + (static_cast<double>(i) / std::max<std::size_t>(win_sizes.size() - 1, 1)) * 320.0;
    const double y = 360.0 - std::min(fanos[i], 5.0) / 5.0 * 310.0;
    pts1 << x << "," << y << " ";
    os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"3.5\" fill=\"#1d4ed8\"/>\n";
  }
  os << "  <polyline fill=\"none\" stroke=\"#1d4ed8\" stroke-width=\"2\" points=\"" << pts1.str() << "\"/>\n";

  std::vector<double> s = res_exp1;
  std::sort(s.begin(), s.end());
  const std::size_t step = std::max<std::size_t>(s.size() / 80, 1);
  os << "  <line x1=\"490\" y1=\"360\" x2=\"850\" y2=\"50\" stroke=\"#0f172a\" stroke-dasharray=\"5,4\"/>\n";
  for (std::size_t i = 0; i < s.size(); i += step) {
    const double p = (static_cast<double>(i) + 0.5) / static_cast<double>(s.size());
    const double q_theo = -std::log(1.0 - p);
    const double x = 490.0 + std::min(q_theo, 6.0) / 6.0 * 360.0;
    const double y = 360.0 - std::min(s[i], 6.0) / 6.0 * 310.0;
    os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"2.2\" fill=\"#059669\"/>\n";
  }
  os << "</svg>\n";
  std::ofstream(path) << os.str();
}

void write_mle_scatter_svg(const std::string& path, const HawkesMleMcSummary& mle) {
  std::ostringstream os;
  write_header(os, 820, 440, "Hawkes MLE (R=200, T=5000): Sebaran Estimasi alpha vs beta");
  os << "  <rect x=\"80\" y=\"50\" width=\"680\" height=\"330\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  for (const auto& e : mle.estimates) {
    const double x = 80.0 + std::clamp((e.alpha - 0.45) / 0.30, 0.0, 1.0) * 680.0;
    const double y = 380.0 - std::clamp((e.beta - 0.90) / 0.60, 0.0, 1.0) * 330.0;
    os << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"3\" fill=\"#1d4ed8\" fill-opacity=\"0.6\"/>\n";
  }
  const double tx = 80.0 + ((0.60 - 0.45) / 0.30) * 680.0;
  const double ty = 380.0 - ((1.20 - 0.90) / 0.60) * 330.0;
  os << "  <circle cx=\"" << tx << "\" cy=\"" << ty << "\" r=\"6\" fill=\"#dc2626\"/>\n</svg>\n";
  std::ofstream(path) << os.str();
}

void write_latency_svg(const std::string& path, const OrderFlowBenchmark& ob) {
  std::ostringstream os;
  write_header(os, 820, 420, "Order Book Latensi (ns) pada 1.000.000 Order");
  os << "  <rect x=\"80\" y=\"50\" width=\"700\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";
  const double vals[4] = {ob.p50_ns, ob.p90_ns, ob.p99_ns, ob.p999_ns};
  const char* lbls[4] = {"p50", "p90", "p99", "p99.9"};
  const double max_v = std::max(ob.p999_ns, 1.0);
  for (int i = 0; i < 4; ++i) {
    const double bh = (vals[i] / max_v) * 280.0;
    const double x = 140.0 + i * 160.0;
    os << "  <rect x=\"" << x << "\" y=\"" << (360.0 - bh) << "\" width=\"75\" height=\"" << bh
       << "\" fill=\"#059669\"/>\n"
       << "  <text x=\"" << (x + 37.5) << "\" y=\"382\" text-anchor=\"middle\" font-family=\"sans-serif\" "
       << "font-size=\"12\">" << lbls[i] << "</text>\n"
       << "  <text x=\"" << (x + 37.5) << "\" y=\"" << (353.0 - bh)
       << "\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\">"
       << static_cast<long long>(vals[i]) << " ns</text>\n";
  }
  os << "</svg>\n";
  std::ofstream(path) << os.str();
}

void write_market_stylized_svg(const std::string& path, const MarketSimResult& mkt) {
  std::ostringstream os;
  write_header(os, 900, 420, "Fakta Stilistik: Profil Kedalaman Buku & ACF |Return|");
  os << "  <rect x=\"60\" y=\"50\" width=\"360\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n"
     << "  <rect x=\"490\" y=\"50\" width=\"360\" height=\"310\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";

  double max_d = 1.0;
  for (double d : mkt.depth_profile) {
    max_d = std::max(max_d, d);
  }
  for (std::size_t i = 0; i < mkt.depth_profile.size(); ++i) {
    const double bh = (mkt.depth_profile[i] / max_d) * 280.0;
    const double x = 75.0 + i * 33.0;
    os << "  <rect x=\"" << x << "\" y=\"" << (360.0 - bh) << "\" width=\"22\" height=\"" << bh
       << "\" fill=\"#6366f1\"/>\n";
  }

  double max_acf = 0.05;
  for (std::size_t i = 1; i < mkt.abs_ret_acf.size(); ++i) {
    max_acf = std::max(max_acf, std::abs(mkt.abs_ret_acf[i]));
  }
  for (std::size_t lag = 1; lag < mkt.abs_ret_acf.size(); ++lag) {
    const double bh = (std::max(mkt.abs_ret_acf[lag], 0.0) / max_acf) * 280.0;
    const double x = 495.0 + (lag - 1) * 7.0;
    os << "  <rect x=\"" << x << "\" y=\"" << (360.0 - bh) << "\" width=\"5\" height=\"" << bh
       << "\" fill=\"#dc2626\"/>\n";
  }
  os << "</svg>\n";
  std::ofstream(path) << os.str();
}

} // namespace hob::svg
