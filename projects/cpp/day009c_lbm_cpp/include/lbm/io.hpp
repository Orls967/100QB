#pragma once

#include "lbm/svg.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace lbm {

std::uint64_t fnv1a_bytes(const void* data, std::size_t n, std::uint64_t h = 14695981039346656037ULL) noexcept;
std::string hex64(std::uint64_t val);

struct StageRecord {
    std::string id;
    std::string name;
    double elapsed_sec = 0.0;
    int files_written = 0;
};

struct CriterionRecord {
    std::string id;
    std::string status;
    std::string measured;
    std::string threshold;
};

struct VacuousFlag {
    std::string stage;
    std::string file;
    std::string reason;
};

struct FileEntry {
    std::string path;
    std::string kind;
    std::size_t bytes = 0;
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::uint64_t fnv1a = 0;
    bool deterministic = true;
};

struct RunContext {
    std::string scale = "full";
    int threads = 4;
    std::size_t cpp_lines = 0;
    std::string compiler_version = "c++20";
    std::string arch = "arm64";
    int cores = 10;
    double total_elapsed_sec = 0.0;
    bool selftest_ok = false;

    std::uint64_t det_csv_hash = 14695981039346656037ULL;
    bool threads_1_2_4_8_identical = false;
    bool two_runs_identical = false;
    std::uint64_t thread_hash_1 = 0;
    std::uint64_t thread_hash_2 = 0;
    std::uint64_t thread_hash_4 = 0;
    std::uint64_t thread_hash_8 = 0;

    std::vector<StageRecord> stages;
    std::vector<CriterionRecord> criteria;
    std::vector<VacuousFlag> vacuous_flags;
    std::vector<FileEntry> files;

    bool is_small() const noexcept { return scale == "small"; }

    void write_csv(const std::string& stage_id,
                   const std::string& rel_path,
                   const std::vector<std::string>& headers,
                   const std::vector<std::vector<std::string>>& rows,
                   const std::vector<std::vector<double>>& numeric_cols,
                   bool deterministic = true);

    void write_svg_xy(const std::string& stage_id,
                      const std::string& rel_path,
                      const std::string& title,
                      const std::string& subtitle,
                      const std::string& xlabel,
                      const std::string& ylabel,
                      const std::vector<SvgSeries>& series,
                      bool log_x = false,
                      bool log_y = false,
                      bool has_ref_y = false,
                      double ref_y = 0.0,
                      const std::string& ref_label = "",
                      bool deterministic = true);

    void write_svg_bar(const std::string& stage_id,
                       const std::string& rel_path,
                       const std::string& title,
                       const std::string& subtitle,
                       const std::string& ylabel,
                       const std::vector<SvgBar>& bars,
                       bool deterministic = true);

    void write_svg_heatmap(const std::string& stage_id,
                           const std::string& rel_path,
                           const std::string& title,
                           const std::string& subtitle,
                           const std::string& xlabel,
                           const std::string& ylabel,
                           int nx, int ny,
                           const std::vector<double>& values,
                           const std::string& cbar_label = "",
                           bool deterministic = true);

    void record_stage(const std::string& id, const std::string& name, double elapsed_sec, int files_written);
    void set_criterion(const std::string& id, const std::string& status,
                       const std::string& measured, const std::string& threshold);
    void add_flag(const std::string& stage_id, const std::string& file, const std::string& reason);

    void write_results_json(const std::string& path) const;
};

std::string fmt_f8(double v);
std::string fmt_e6(double v);

}  // namespace lbm
