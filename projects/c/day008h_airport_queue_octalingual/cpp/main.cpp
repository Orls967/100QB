#include "common.hpp"
#include <chrono>
#include <filesystem>
namespace fs = std::filesystem;

static void manifest_scan(const fs::path &base, const fs::path &cur, FILE *out) {
    std::vector<fs::directory_entry> ents;
    for (const auto &e : fs::directory_iterator(cur)) ents.push_back(e);
    std::sort(ents.begin(), ents.end(), [](const auto &a, const auto &b) { return a.path().filename() < b.path().filename(); });
    for (const auto &e : ents) {
        auto nm = e.path().filename().string();
        if (nm[0] == '.') continue;
        if (e.is_directory()) manifest_scan(base, e.path(), out);
        else if (nm != "files_manifest.csv") {
            FILE *f = std::fopen(e.path().c_str(), "rb"); if (!f) continue;
            std::fseek(f, 0, SEEK_END); size_t sz = static_cast<size_t>(std::ftell(f)); std::fseek(f, 0, SEEK_SET);
            std::vector<uint8_t> buf(sz); size_t rd = std::fread(buf.data(), 1, sz, f); std::fclose(f);
            int64_t lc = 0; for (size_t i = 0; i < rd; ++i) if (buf[i] == '\n') lc++;
            auto rel = fs::relative(e.path(), base).generic_string();
            std::fprintf(out, "%s,%lld,%zu,%016llx\n", rel.c_str(), static_cast<long long>(lc), rd, static_cast<unsigned long long>(fnv1a_bytes(buf.data(), rd)));
        }
    }
}

int main(int argc, char **argv) {
    const char *cfg_path = "config.kv", *out_dir = "out/cpp";
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--config") && i + 1 < argc) cfg_path = argv[++i];
        else if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out_dir = argv[++i];
    }
    fs::create_directories(std::string(out_dir) + "/figures");
    Config c = load_config(cfg_path);
    if (argc > 1 && !std::strcmp(argv[1], "s01")) { write_s01_to_s13(c, out_dir); return 0; }

    auto clk = []() { return std::chrono::steady_clock::now(); };
    double tot = 0.0;
    for (int st = 1; st <= 13; ++st) std::printf("[cpp] tahap %d/48 (s%02d)\n", st, st);
    auto t0 = clk(); write_s01_to_s13(c, out_dir); tot += std::chrono::duration<double>(clk() - t0).count();

    int64_t main_days = 365;
    std::printf("[cpp] tahap 14/48 (s14)\n");
    t0 = clk(); RunParams rpm = default_rp(c, main_days); rpm.audit_on = 1;
    auto s_main = sim_create(c, rpm); sim_run(*s_main); tot += std::chrono::duration<double>(clk() - t0).count();

    for (int st = 15; st <= 41; ++st) std::printf("[cpp] tahap %d/48 (s%02d)\n", st, st);
    t0 = clk(); write_s14_to_s45(c, out_dir, *s_main); tot += std::chrono::duration<double>(clk() - t0).count();

    std::printf("[cpp] tahap 42/48 (s42)\n");
    std::vector<double> vals(5), dev(5);
    for (int r = 0; r < 5; ++r) {
        auto sb = sim_create(c, default_rp(c, 3)); auto tb = clk(); sim_run(*sb);
        double dt = std::max(std::chrono::duration<double>(clk() - tb).count(), 1e-9);
        vals[r] = static_cast<double>(sb->sys_sum_L) / dt;
    }
    std::sort(vals.begin(), vals.end()); double med = vals[2];
    for (int r = 0; r < 5; ++r) dev[r] = std::fabs(vals[r] - med);
    std::sort(dev.begin(), dev.end());
    FILE *fb = std::fopen((std::string(out_dir) + "/s42_bench.csv").c_str(), "w");
    std::fprintf(fb, "slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,%.2f,%.2f\n", med, dev[2]);
    std::fclose(fb);

    std::printf("[cpp] tahap 43/48 (s43)\n[cpp] tahap 44/48 (s44)\n[cpp] tahap 45/48 (s45)\n[cpp] tahap 46/48 (s46)\n");
    t0 = clk(); write_s46_svgs(c, out_dir, *s_main); tot += std::chrono::duration<double>(clk() - t0).count();

    std::printf("[cpp] tahap 47/48 (s47)\n[cpp] tahap 48/48 (s48)\n");
    FILE *fs = std::fopen((std::string(out_dir) + "/s48_summary.json").c_str(), "w");
    std::fprintf(fs, "{\"lang\":\"cpp\",\"days\":%lld,\"sw_days\":3,\"toolchain\":\"Apple clang 16.0.0 (-std=c++20)\",\"total_sec\":%.3f}\n", static_cast<long long>(main_days), tot);
    std::fclose(fs);

    FILE *fm = std::fopen((std::string(out_dir) + "/files_manifest.csv").c_str(), "w");
    std::fprintf(fm, "file,lines,bytes,fnv1a\n");
    manifest_scan(out_dir, out_dir, fm);
    std::fclose(fm);
    return 0;
}
