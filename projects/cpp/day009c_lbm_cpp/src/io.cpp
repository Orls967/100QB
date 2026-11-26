#include "lbm/io.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>

namespace lbm {

namespace {

bool is_exempt_col(const std::string& h) {
    static const std::set<std::string> exempt = {
        "stage", "step", "steps", "t", "n", "nx", "ny", "nz", "diameter", "res",
        "op", "operator", "tau", "lambda", "re", "ma", "u0", "u_wind", "u_lid",
        "u_in", "fx", "rho0", "scheme", "layout", "tile_w", "threads", "repeat",
        "version", "label", "dir", "q_frac", "bouzidi", "quant_step", "n_cells",
        "bytes", "rows", "cols", "files", "pass", "ok", "violations", "eddies",
        "survived", "bin", "y", "x", "z", "k", "df", "exact_roundtrip", "max_q_err",
        "raw_bytes", "comp_bytes", "file", "kind", "fnv1a", "id", "name", "status", "reason"
    };
    return exempt.count(h) > 0;
}

std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out.push_back(c); break;
        }
    }
    return out;
}

}  // namespace

std::uint64_t fnv1a_bytes(const void* data, std::size_t n, std::uint64_t h) noexcept {
    const auto* p = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

std::string hex64(std::uint64_t val) {
    char buf[20];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(val));
    return buf;
}

std::string fmt_f8(double v) {
    if (!std::isfinite(v)) return std::isnan(v) ? "NaN" : (v > 0 ? "Inf" : "-Inf");
    if (v != 0.0 && (std::abs(v) < 1e-4 || std::abs(v) >= 1e7)) return fmt_e6(v);
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.8f", v);
    return buf;
}

std::string fmt_e6(double v) {
    if (!std::isfinite(v)) return std::isnan(v) ? "NaN" : (v > 0 ? "Inf" : "-Inf");
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.6e", v);
    return buf;
}

void RunContext::write_csv(const std::string& stage_id, const std::string& rel_path,
                           const std::vector<std::string>& headers,
                           const std::vector<std::vector<std::string>>& rows,
                           const std::vector<std::vector<double>>& numeric_cols,
                           bool deterministic) {
    const std::filesystem::path p(rel_path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

    std::ostringstream os;
    for (std::size_t c = 0; c < headers.size(); ++c) {
        if (c) os << ',';
        os << headers[c];
    }
    os << '\n';
    for (const auto& r : rows) {
        for (std::size_t c = 0; c < r.size(); ++c) {
            if (c) os << ',';
            os << r[c];
        }
        os << '\n';
    }
    const std::string content = os.str();
    std::ofstream ofs(rel_path, std::ios::binary);
    ofs.write(content.data(), static_cast<std::streamsize>(content.size()));

    const std::uint64_t h = fnv1a_bytes(content.data(), content.size());
    if (deterministic) {
        det_csv_hash = fnv1a_bytes(&h, sizeof(h), det_csv_hash);
    }

    // Anti-vacuous audit across numeric columns
    for (std::size_t c = 0; c < numeric_cols.size() && c < headers.size(); ++c) {
        const auto& col = numeric_cols[c];
        if (col.empty()) continue;
        bool has_bad = false;
        double vmin = col[0], vmax = col[0];
        for (double v : col) {
            if (!std::isfinite(v)) {
                has_bad = true;
                break;
            }
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
        if (has_bad) {
            add_flag(stage_id, rel_path, "Non-finite value in column " + headers[c]);
        }
        if (headers[c] == "rho" || headers[c] == "rho_min" || headers[c] == "rho_mid" ||
            headers[c] == "rho_in" || headers[c] == "rho_out") {
            if (vmin <= 0.0) {
                add_flag(stage_id, rel_path, "Non-positive density in column " + headers[c]);
            }
        }
        if (rows.size() >= 4 && !is_exempt_col(headers[c]) &&
            headers[c].find("exact") == std::string::npos &&
            headers[c].find("target") == std::string::npos &&
            headers[c].find("analyt") == std::string::npos &&
            headers[c].find("nu_th") == std::string::npos &&
            headers[c].find("err") == std::string::npos &&
            headers[c].find("delta_") == std::string::npos &&
            headers[c].find("drift") == std::string::npos &&
            headers[c].find("pass") == std::string::npos &&
            headers[c].find("violations") == std::string::npos &&
            std::abs(vmax - vmin) <= 1e-12 * std::max(1e-30, std::abs(vmax))) {
            add_flag(stage_id, rel_path, "Constant non-parameter column " + headers[c]);
        }
    }

    files.push_back({rel_path, "csv", content.size(), rows.size(), headers.size(), h, deterministic});
}

void RunContext::write_svg_xy(const std::string& stage_id, const std::string& rel_path,
                              const std::string& title, const std::string& subtitle,
                              const std::string& xlabel, const std::string& ylabel,
                              const std::vector<SvgSeries>& series, bool log_x, bool log_y,
                              bool has_ref_y, double ref_y, const std::string& ref_label,
                              bool deterministic) {
    (void)stage_id;
    SvgWriter::write_xy(rel_path, title, subtitle, xlabel, ylabel, series, log_x, log_y, has_ref_y,
                        ref_y, ref_label);
    std::ifstream ifs(rel_path, std::ios::binary);
    std::string data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    const std::uint64_t h = fnv1a_bytes(data.data(), data.size());
    files.push_back({rel_path, "svg", data.size(), 0, 0, h, deterministic});
}

void RunContext::write_svg_bar(const std::string& stage_id, const std::string& rel_path,
                               const std::string& title, const std::string& subtitle,
                               const std::string& ylabel, const std::vector<SvgBar>& bars,
                               bool deterministic) {
    (void)stage_id;
    SvgWriter::write_bar(rel_path, title, subtitle, ylabel, bars);
    std::ifstream ifs(rel_path, std::ios::binary);
    std::string data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    const std::uint64_t h = fnv1a_bytes(data.data(), data.size());
    files.push_back({rel_path, "svg", data.size(), 0, 0, h, deterministic});
}

void RunContext::write_svg_heatmap(const std::string& stage_id, const std::string& rel_path,
                                   const std::string& title, const std::string& subtitle,
                                   const std::string& xlabel, const std::string& ylabel,
                                   int nx, int ny, const std::vector<double>& values,
                                   const std::string& cbar_label, bool deterministic) {
    (void)stage_id;
    SvgWriter::write_heatmap(rel_path, title, subtitle, xlabel, ylabel, nx, ny, values, cbar_label);
    std::ifstream ifs(rel_path, std::ios::binary);
    std::string data((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    const std::uint64_t h = fnv1a_bytes(data.data(), data.size());
    files.push_back({rel_path, "svg", data.size(), 0, 0, h, deterministic});
}

void RunContext::record_stage(const std::string& id, const std::string& name, double elapsed_sec,
                              int files_written) {
    stages.push_back({id, name, elapsed_sec, files_written});
}

void RunContext::set_criterion(const std::string& id, const std::string& status,
                               const std::string& measured, const std::string& threshold) {
    for (auto& c : criteria) {
        if (c.id == id) {
            c.status = status;
            c.measured = measured;
            c.threshold = threshold;
            return;
        }
    }
    criteria.push_back({id, status, measured, threshold});
}

void RunContext::add_flag(const std::string& stage_id, const std::string& file,
                          const std::string& reason) {
    vacuous_flags.push_back({stage_id, file, reason});
}

void RunContext::write_results_json(const std::string& path) const {
    const std::filesystem::path p(path);
    if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

    std::size_t n_csv = 0, n_svg = 0, csv_b = 0, svg_b = 0;
    for (const auto& f : files) {
        if (f.kind == "csv") { ++n_csv; csv_b += f.bytes; }
        else if (f.kind == "svg") { ++n_svg; svg_b += f.bytes; }
    }
    int n_ter = 0, n_seb = 0, n_oth = 0;
    for (const auto& c : criteria) {
        if (c.status == "tercapai") ++n_ter;
        else if (c.status == "sebagian") ++n_seb;
        else ++n_oth;
    }

    std::ostringstream os;
    os << "{\n"
       << "  \"session\": \"day009c\",\n"
       << "  \"project_path\": \"projects/cpp/day009c_lbm_cpp\",\n"
       << "  \"scale\": \"" << json_escape(scale) << "\",\n"
       << "  \"master_seed\": 20261127,\n"
       << "  \"threads\": " << threads << ",\n"
       << "  \"run_completed_without_errors\": " << (selftest_ok ? "true" : "false") << ",\n"
       << "  \"machine\": {\n"
       << "    \"arch\": \"" << json_escape(arch) << "\",\n"
       << "    \"cores\": " << cores << ",\n"
       << "    \"compiler\": \"" << json_escape(compiler_version) << "\"\n"
       << "  },\n"
       << "  \"code_metrics\": {\n"
       << "    \"cpp_lines\": " << cpp_lines << ",\n"
       << "    \"stages_completed\": " << stages.size() << ",\n"
       << "    \"total_output_files\": " << (files.size() + 1) << ",\n"
       << "    \"csv_files\": " << n_csv << ",\n"
       << "    \"svg_files\": " << n_svg << ",\n"
       << "    \"csv_total_bytes\": " << csv_b << ",\n"
       << "    \"svg_total_bytes\": " << svg_b << ",\n"
       << "    \"total_elapsed_sec\": " << fmt_f8(total_elapsed_sec) << "\n"
       << "  },\n"
       << "  \"determinism\": {\n"
       << "    \"combined_deterministic_csv_fnv1a\": \"" << hex64(det_csv_hash) << "\",\n"
       << "    \"threads_1_2_4_8_identical\": " << (threads_1_2_4_8_identical ? "true" : "false") << ",\n"
       << "    \"two_runs_identical\": " << (two_runs_identical ? "true" : "false") << ",\n"
       << "    \"thread_hash_1\": \"" << hex64(thread_hash_1) << "\",\n"
       << "    \"thread_hash_2\": \"" << hex64(thread_hash_2) << "\",\n"
       << "    \"thread_hash_4\": \"" << hex64(thread_hash_4) << "\",\n"
       << "    \"thread_hash_8\": \"" << hex64(thread_hash_8) << "\"\n"
       << "  },\n"
       << "  \"criteria_summary\": {\n"
       << "    \"total\": " << criteria.size() << ",\n"
       << "    \"tercapai\": " << n_ter << ",\n"
       << "    \"sebagian\": " << n_seb << ",\n"
       << "    \"lainnya\": " << n_oth << "\n"
       << "  },\n"
       << "  \"criteria\": {\n";
    for (std::size_t i = 0; i < criteria.size(); ++i) {
        const auto& c = criteria[i];
        os << "    \"" << json_escape(c.id) << "\": {\n"
           << "      \"status\": \"" << json_escape(c.status) << "\",\n"
           << "      \"measured\": \"" << json_escape(c.measured) << "\",\n"
           << "      \"threshold\": \"" << json_escape(c.threshold) << "\"\n"
           << "    }" << (i + 1 < criteria.size() ? "," : "") << "\n";
    }
    os << "  },\n"
       << "  \"vacuous_flags\": [";
    for (std::size_t i = 0; i < vacuous_flags.size(); ++i) {
        const auto& fl = vacuous_flags[i];
        os << (i ? ",\n" : "\n")
           << "    {\"stage\": \"" << json_escape(fl.stage) << "\", \"file\": \""
           << json_escape(fl.file) << "\", \"reason\": \"" << json_escape(fl.reason) << "\"}";
    }
    if (!vacuous_flags.empty()) os << "\n  ";
    os << "],\n"
       << "  \"stages\": [\n";
    for (std::size_t i = 0; i < stages.size(); ++i) {
        const auto& s = stages[i];
        os << "    {\"id\": \"" << json_escape(s.id) << "\", \"name\": \""
           << json_escape(s.name) << "\", \"elapsed_sec\": " << fmt_f8(s.elapsed_sec)
           << ", \"files\": " << s.files_written << "}"
           << (i + 1 < stages.size() ? "," : "") << "\n";
    }
    os << "  ]\n"
       << "}\n";

    std::ofstream ofs(path, std::ios::binary);
    const std::string out = os.str();
    ofs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

}  // namespace lbm
