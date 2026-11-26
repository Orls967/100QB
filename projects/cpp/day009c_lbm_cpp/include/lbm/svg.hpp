#pragma once

#include <string>
#include <vector>

namespace lbm {

struct SvgSeries {
    std::string label;
    std::string color;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> yerr;
    bool line = true;
    bool points = true;
    bool dashed = false;
};

struct SvgBar {
    std::string label;
    double value = 0.0;
    double err = 0.0;
    std::string color = "#2563eb";
};

class SvgWriter {
public:
    static std::string escape_xml(const std::string& s);

    static bool write_xy(const std::string& path,
                         const std::string& title,
                         const std::string& subtitle,
                         const std::string& xlabel,
                         const std::string& ylabel,
                         const std::vector<SvgSeries>& series,
                         bool log_x = false,
                         bool log_y = false,
                         bool has_ref_y = false,
                         double ref_y = 0.0,
                         const std::string& ref_label = "");

    static bool write_bar(const std::string& path,
                          const std::string& title,
                          const std::string& subtitle,
                          const std::string& ylabel,
                          const std::vector<SvgBar>& bars);

    static bool write_heatmap(const std::string& path,
                              const std::string& title,
                              const std::string& subtitle,
                              const std::string& xlabel,
                              const std::string& ylabel,
                              int nx, int ny,
                              const std::vector<double>& values,
                              const std::string& cbar_label = "");
};

}  // namespace lbm
