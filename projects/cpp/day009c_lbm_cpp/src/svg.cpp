#include "lbm/svg.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace lbm {

namespace {

std::string fmt_num(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3g", v);
    return buf;
}

std::string fmt_coord(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.1f", v);
    return buf;
}

double safe_tx(double v, bool is_log) {
    if (!is_log) return v;
    return std::log10(v > 1e-30 ? v : 1e-30);
}

void write_header(std::ostringstream& os, int w, int h, const std::string& title,
                  const std::string& subtitle) {
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
       << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << w << "\" height=\"" << h
       << "\" viewBox=\"0 0 " << w << " " << h << "\">\n"
       << "  <rect x=\"0\" y=\"0\" width=\"" << w << "\" height=\"" << h << "\" fill=\"#fafbfc\"/>\n"
       << "  <rect x=\"1\" y=\"1\" width=\"" << (w - 2) << "\" height=\"" << (h - 2)
       << "\" fill=\"none\" stroke=\"#d8dee9\" stroke-width=\"1\"/>\n"
       << "  <text x=\"" << (w / 2) << "\" y=\"28\" text-anchor=\"middle\" "
       << "font-family=\"Helvetica,Arial,sans-serif\" font-size=\"15\" font-weight=\"bold\" fill=\"#1f2937\">"
       << SvgWriter::escape_xml(title) << "</text>\n";
    if (!subtitle.empty()) {
        os << "  <text x=\"" << (w / 2) << "\" y=\"46\" text-anchor=\"middle\" "
           << "font-family=\"Helvetica,Arial,sans-serif\" font-size=\"11\" fill=\"#4b5563\">"
           << SvgWriter::escape_xml(subtitle) << "</text>\n";
    }
}

bool flush_file(const std::string& path, const std::string& content) {
    const std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) return false;
    ofs.write(content.data(), static_cast<std::streamsize>(content.size()));
    return ofs.good();
}

}  // namespace

std::string SvgWriter::escape_xml(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out.push_back(c); break;
        }
    }
    return out;
}

bool SvgWriter::write_xy(const std::string& path, const std::string& title,
                         const std::string& subtitle, const std::string& xlabel,
                         const std::string& ylabel, const std::vector<SvgSeries>& series,
                         bool log_x, bool log_y, bool has_ref_y, double ref_y,
                         const std::string& ref_label) {
    constexpr int W = 760, H = 460;
    constexpr int ML = 78, MR = 32, MT = 62, MB = 60;
    constexpr int PW = W - ML - MR;
    constexpr int PH = H - MT - MB;

    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    bool any = false;
    for (const auto& sr : series) {
        const std::size_t n = std::min(sr.x.size(), sr.y.size());
        for (std::size_t i = 0; i < n; ++i) {
            const double xv = sr.x[i];
            const double yv = sr.y[i];
            if (!std::isfinite(xv) || !std::isfinite(yv)) continue;
            if (log_x && xv <= 0.0) continue;
            if (log_y && yv <= 0.0) continue;
            const double tx = safe_tx(xv, log_x);
            const double ty = safe_tx(yv, log_y);
            xmin = std::min(xmin, tx);
            xmax = std::max(xmax, tx);
            ymin = std::min(ymin, ty);
            ymax = std::max(ymax, ty);
            if (!log_y && i < sr.yerr.size() && std::isfinite(sr.yerr[i])) {
                ymin = std::min(ymin, ty - std::abs(sr.yerr[i]));
                ymax = std::max(ymax, ty + std::abs(sr.yerr[i]));
            }
            any = true;
        }
    }
    if (has_ref_y && (!log_y || ref_y > 0.0)) {
        const double ry = safe_tx(ref_y, log_y);
        ymin = std::min(ymin, ry);
        ymax = std::max(ymax, ry);
        any = true;
    }
    if (!any) {
        xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0;
    }
    if (xmin >= xmax) { xmin -= 1.0; xmax += 1.0; }
    if (ymin >= ymax) { ymin -= 1.0; ymax += 1.0; }
    const double xpad = 0.04 * (xmax - xmin);
    const double ypad = 0.08 * (ymax - ymin);
    xmin -= xpad; xmax += xpad;
    ymin -= ypad; ymax += ypad;

    auto map_x = [&](double tx) { return ML + (tx - xmin) / (xmax - xmin) * PW; };
    auto map_y = [&](double ty) { return MT + PH - (ty - ymin) / (ymax - ymin) * PH; };

    std::ostringstream os;
    write_header(os, W, H, title, subtitle);
    os << "  <rect x=\"" << ML << "\" y=\"" << MT << "\" width=\"" << PW << "\" height=\"" << PH
       << "\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n";

    for (int t = 0; t <= 5; ++t) {
        const double frac = t / 5.0;
        const double gx = ML + frac * PW;
        const double vx = xmin + frac * (xmax - xmin);
        const double dx = log_x ? std::pow(10.0, vx) : vx;
        os << "  <line x1=\"" << fmt_coord(gx) << "\" y1=\"" << MT << "\" x2=\"" << fmt_coord(gx)
           << "\" y2=\"" << (MT + PH) << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n"
           << "  <text x=\"" << fmt_coord(gx) << "\" y=\"" << (MT + PH + 16)
           << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
           << fmt_num(dx) << "</text>\n";

        const double gy = MT + PH - frac * PH;
        const double vy = ymin + frac * (ymax - ymin);
        const double dy = log_y ? std::pow(10.0, vy) : vy;
        os << "  <line x1=\"" << ML << "\" y1=\"" << fmt_coord(gy) << "\" x2=\"" << (ML + PW)
           << "\" y2=\"" << fmt_coord(gy) << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n"
           << "  <text x=\"" << (ML - 8) << "\" y=\"" << fmt_coord(gy + 3.5)
           << "\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
           << fmt_num(dy) << "</text>\n";
    }

    if (has_ref_y && (!log_y || ref_y > 0.0)) {
        const double py = map_y(safe_tx(ref_y, log_y));
        os << "  <line x1=\"" << ML << "\" y1=\"" << fmt_coord(py) << "\" x2=\"" << (ML + PW)
           << "\" y2=\"" << fmt_coord(py)
           << "\" stroke=\"#6b7280\" stroke-width=\"1.4\" stroke-dasharray=\"5,4\"/>\n";
        if (!ref_label.empty()) {
            os << "  <text x=\"" << (ML + PW - 6) << "\" y=\"" << fmt_coord(py - 5.0)
               << "\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
               << escape_xml(ref_label) << "</text>\n";
        }
    }

    for (const auto& sr : series) {
        const std::size_t n = std::min(sr.x.size(), sr.y.size());
        std::vector<std::pair<double, double>> pts;
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(sr.x[i]) || !std::isfinite(sr.y[i])) continue;
            if (log_x && sr.x[i] <= 0.0) continue;
            if (log_y && sr.y[i] <= 0.0) continue;
            pts.emplace_back(map_x(safe_tx(sr.x[i], log_x)), map_y(safe_tx(sr.y[i], log_y)));
        }
        if (sr.line && pts.size() >= 2) {
            os << "  <polyline fill=\"none\" stroke=\"" << sr.color << "\" stroke-width=\"2\"";
            if (sr.dashed) os << " stroke-dasharray=\"6,4\"";
            os << " points=\"";
            for (std::size_t i = 0; i < pts.size(); ++i) {
                if (i) os << ' ';
                os << fmt_coord(pts[i].first) << ',' << fmt_coord(pts[i].second);
            }
            os << "\"/>\n";
        }
        if (!log_y && !sr.yerr.empty()) {
            for (std::size_t i = 0; i < n && i < sr.yerr.size(); ++i) {
                if (!std::isfinite(sr.x[i]) || !std::isfinite(sr.y[i]) || !std::isfinite(sr.yerr[i])) continue;
                const double px = map_x(safe_tx(sr.x[i], log_x));
                const double py1 = map_y(sr.y[i] - sr.yerr[i]);
                const double py2 = map_y(sr.y[i] + sr.yerr[i]);
                os << "  <line x1=\"" << fmt_coord(px) << "\" y1=\"" << fmt_coord(py1)
                   << "\" x2=\"" << fmt_coord(px) << "\" y2=\"" << fmt_coord(py2)
                   << "\" stroke=\"" << sr.color << "\" stroke-width=\"1.2\"/>\n";
            }
        }
        if (sr.points) {
            for (const auto& pt : pts) {
                os << "  <circle cx=\"" << fmt_coord(pt.first) << "\" cy=\"" << fmt_coord(pt.second)
                   << "\" r=\"3.2\" fill=\"" << sr.color << "\"/>\n";
            }
        }
    }

    int ly = MT + 14;
    for (const auto& sr : series) {
        if (sr.label.empty()) continue;
        const int lx = ML + PW - 190;
        os << "  <line x1=\"" << lx << "\" y1=\"" << ly << "\" x2=\"" << (lx + 18) << "\" y2=\"" << ly
           << "\" stroke=\"" << sr.color << "\" stroke-width=\"2.2\"/>\n"
           << "  <text x=\"" << (lx + 24) << "\" y=\"" << (ly + 4)
           << "\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10.5\" fill=\"#1f2937\">"
           << escape_xml(sr.label) << "</text>\n";
        ly += 16;
    }

    os << "  <text x=\"" << (ML + PW / 2) << "\" y=\"" << (H - 16)
       << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\">"
       << escape_xml(xlabel) << "</text>\n"
       << "  <text x=\"20\" y=\"" << (MT + PH / 2)
       << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" "
       << "transform=\"rotate(-90 20 " << (MT + PH / 2) << ")\">" << escape_xml(ylabel) << "</text>\n"
       << "</svg>\n";
    return flush_file(path, os.str());
}

bool SvgWriter::write_bar(const std::string& path, const std::string& title,
                          const std::string& subtitle, const std::string& ylabel,
                          const std::vector<SvgBar>& bars) {
    constexpr int W = 760, H = 440;
    constexpr int ML = 78, MR = 28, MT = 62, MB = 78;
    constexpr int PW = W - ML - MR;
    constexpr int PH = H - MT - MB;

    double ymax = 1e-12;
    for (const auto& b : bars) ymax = std::max(ymax, b.value + std::abs(b.err));
    ymax *= 1.15;

    std::ostringstream os;
    write_header(os, W, H, title, subtitle);
    os << "  <rect x=\"" << ML << "\" y=\"" << MT << "\" width=\"" << PW << "\" height=\"" << PH
       << "\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n";

    for (int t = 0; t <= 5; ++t) {
        const double frac = t / 5.0;
        const double gy = MT + PH - frac * PH;
        os << "  <line x1=\"" << ML << "\" y1=\"" << fmt_coord(gy) << "\" x2=\"" << (ML + PW)
           << "\" y2=\"" << fmt_coord(gy) << "\" stroke=\"#e5e7eb\" stroke-width=\"1\"/>\n"
           << "  <text x=\"" << (ML - 8) << "\" y=\"" << fmt_coord(gy + 3.5)
           << "\" text-anchor=\"end\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
           << fmt_num(frac * ymax) << "</text>\n";
    }

    const int n = static_cast<int>(bars.size());
    const double step = static_cast<double>(PW) / std::max(1, n);
    const double bw = step * 0.62;
    for (int i = 0; i < n; ++i) {
        const auto& b = bars[static_cast<std::size_t>(i)];
        const double x = ML + i * step + 0.5 * (step - bw);
        const double bh = (std::max(0.0, b.value) / ymax) * PH;
        const double y = MT + PH - bh;
        os << "  <rect x=\"" << fmt_coord(x) << "\" y=\"" << fmt_coord(y)
           << "\" width=\"" << fmt_coord(bw) << "\" height=\"" << fmt_coord(bh)
           << "\" fill=\"" << b.color << "\" rx=\"2\"/>\n"
           << "  <text x=\"" << fmt_coord(x + 0.5 * bw) << "\" y=\"" << (MT + PH + 16)
           << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
           << escape_xml(b.label) << "</text>\n";
    }
    os << "  <text x=\"20\" y=\"" << (MT + PH / 2)
       << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" "
       << "transform=\"rotate(-90 20 " << (MT + PH / 2) << ")\">" << escape_xml(ylabel) << "</text>\n"
       << "</svg>\n";
    return flush_file(path, os.str());
}

bool SvgWriter::write_heatmap(const std::string& path, const std::string& title,
                              const std::string& subtitle, const std::string& xlabel,
                              const std::string& ylabel, int nx, int ny,
                              const std::vector<double>& values, const std::string& cbar_label) {
    constexpr int W = 720, H = 560;
    constexpr int ML = 78, MR = 80, MT = 62, MB = 60;
    constexpr int PW = W - ML - MR;
    constexpr int PH = H - MT - MB;

    double vmin = 1e300, vmax = -1e300;
    for (double v : values) {
        if (!std::isfinite(v)) continue;
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }
    if (vmin >= vmax) { vmin -= 1e-6; vmax += 1e-6; }

    std::ostringstream os;
    write_header(os, W, H, title, subtitle);
    os << "  <rect x=\"" << ML << "\" y=\"" << MT << "\" width=\"" << PW << "\" height=\"" << PH
       << "\" fill=\"#ffffff\" stroke=\"#9ca3af\" stroke-width=\"1\"/>\n";

    // Downsample if > 64x64 to keep SVG lightweight and crisp
    const int sx = std::max(1, nx / 64);
    const int sy = std::max(1, ny / 64);
    const int gx_n = (nx + sx - 1) / sx;
    const int gy_n = (ny + sy - 1) / sy;
    const double cw = static_cast<double>(PW) / std::max(1, gx_n);
    const double ch = static_cast<double>(PH) / std::max(1, gy_n);

    for (int gy = 0; gy < gy_n; ++gy) {
        for (int gx = 0; gx < gx_n; ++gx) {
            const int x = std::min(nx - 1, gx * sx);
            const int y = std::min(ny - 1, gy * sy);
            const double v = values[static_cast<std::size_t>(y) * static_cast<std::size_t>(nx) +
                                    static_cast<std::size_t>(x)];
            const double t = std::clamp((v - vmin) / (vmax - vmin), 0.0, 1.0);
            // Blue -> White/Cyan -> Amber/Red colormap
            const int r = static_cast<int>(std::clamp(255.0 * (0.15 + 0.85 * t), 0.0, 255.0));
            const int g = static_cast<int>(std::clamp(255.0 * (0.25 + 0.55 * (1.0 - std::abs(2.0 * t - 1.0))), 0.0, 255.0));
            const int b = static_cast<int>(std::clamp(255.0 * (0.90 - 0.78 * t), 0.0, 255.0));
            char col[16];
            std::snprintf(col, sizeof(col), "#%02x%02x%02x", r, g, b);
            const double px = ML + gx * cw;
            const double py = MT + PH - (gy + 1) * ch;
            os << "  <rect x=\"" << fmt_coord(px) << "\" y=\"" << fmt_coord(py)
               << "\" width=\"" << fmt_coord(cw + 0.6) << "\" height=\"" << fmt_coord(ch + 0.6)
               << "\" fill=\"" << col << "\"/>\n";
        }
    }

    // Colorbar legend
    os << "  <text x=\"" << (ML + PW + 12) << "\" y=\"" << (MT + 10)
       << "\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#111827\">"
       << fmt_num(vmax) << "</text>\n"
       << "  <text x=\"" << (ML + PW + 12) << "\" y=\"" << (MT + PH)
       << "\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#111827\">"
       << fmt_num(vmin) << "</text>\n";
    if (!cbar_label.empty()) {
        os << "  <text x=\"" << (ML + PW + 12) << "\" y=\"" << (MT + PH / 2)
           << "\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10\" fill=\"#374151\">"
           << escape_xml(cbar_label) << "</text>\n";
    }

    os << "  <text x=\"" << (ML + PW / 2) << "\" y=\"" << (H - 16)
       << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\">"
       << escape_xml(xlabel) << "</text>\n"
       << "  <text x=\"20\" y=\"" << (MT + PH / 2)
       << "\" text-anchor=\"middle\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"12\" fill=\"#111827\" "
       << "transform=\"rotate(-90 20 " << (MT + PH / 2) << ")\">" << escape_xml(ylabel) << "</text>\n"
       << "</svg>\n";
    return flush_file(path, os.str());
}

}  // namespace lbm
