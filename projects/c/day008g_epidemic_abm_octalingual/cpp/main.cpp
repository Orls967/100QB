#include "common.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <dirent.h>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <vector>

static double now_s() { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static void scan_dir(const std::string &base, const std::string &rel, FILE *mf) {
    std::string path = rel.empty() ? base : (base + "/" + rel);
    DIR *d = opendir(path.c_str()); if (!d) return;
    std::vector<std::string> names;
    while (struct dirent *ent = readdir(d)) if (ent->d_name[0] != '.') names.emplace_back(ent->d_name);
    closedir(d);
    std::sort(names.begin(), names.end());
    for (const auto &nm : names) {
        std::string subrel = rel.empty() ? nm : (rel + "/" + nm), full = base + "/" + subrel;
        struct stat st{};
        if (stat(full.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) scan_dir(base, subrel, mf);
        else if (nm != "files_manifest.csv") {
            FILE *fp = fopen(full.c_str(), "rb");
            if (fp) {
                fseek(fp, 0, SEEK_END); long sz = ftell(fp); fseek(fp, 0, SEEK_SET);
                std::vector<uint8_t> buf(sz > 0 ? sz : 0);
                if (sz > 0) fread(buf.data(), 1, buf.size(), fp); fclose(fp);
                long lines = 0; for (uint8_t b : buf) if (b == '\n') lines++;
                fprintf(mf, "%s,%ld,%zu,%016llx\n", subrel.c_str(), lines, buf.size(), static_cast<unsigned long long>(fnv1a_bytes(buf.data(), buf.size())));
            }
        }
    }
}

int main(int argc, char **argv) {
    const char *cmd = argc > 1 ? argv[1] : "all";
    const char *cfg_path = "config.kv", *out_dir = "out/cpp";
    for (int i = 1; i + 1 < argc; ++i) {
        if (!strcmp(argv[i], "--config")) cfg_path = argv[i + 1];
        else if (!strcmp(argv[i], "--out")) out_dir = argv[i + 1];
    }
    Config cfg; if (!load_config(cfg_path, cfg)) return 1;
    mkdir("out", 0755); mkdir(out_dir, 0755); mkdir((std::string(out_dir) + "/figures").c_str(), 0755);
    if (!strcmp(cmd, "s01")) { write_s01_to_s09(cfg, out_dir); return 0; }

    double tot = 0.0;
    for (int st = 1; st <= 9; ++st) printf("[cpp] tahap %d/39 (s%02d)\n", st, st);
    double t0 = now_s(); write_s01_to_s09(cfg, out_dir); tot += now_s() - t0;

    int64_t N_main = 100000, beta_cal = 400;
    printf("[cpp] tahap 10/39 (s10)\n");
    t0 = now_s(); RunParams rp_no; rp_no.N = N_main; rp_no.days = cfg.main_days; rp_no.beta_bp = beta_cal; rp_no.interv_on = 0;
    Sim s_no; sim_init(s_no, cfg, rp_no, 1); sim_run(s_no); tot += now_s() - t0;

    printf("[cpp] tahap 11/39 (s11)\n");
    t0 = now_s(); RunParams rp_main = rp_no; rp_main.interv_on = 1;
    Sim s_main; sim_init(s_main, cfg, rp_main, 1); sim_run(s_main); tot += now_s() - t0;

    for (int st = 12; st <= 34; ++st) printf("[cpp] tahap %d/39 (s%02d)\n", st, st);
    t0 = now_s(); write_s10_to_s37(s_no, s_main, cfg, out_dir, beta_cal); tot += now_s() - t0;

    printf("[cpp] tahap 35/39 (s35)\n");
    std::vector<double> vals(5), dev(5);
    for (int rp = 0; rp < 5; ++rp) {
        RunParams rpb; rpb.N = 2000; rpb.days = 60; rpb.beta_bp = beta_cal; rpb.interv_on = 1;
        Sim sb; sim_init(sb, cfg, rpb, 0);
        double tb = now_s(); sim_run(sb); double dt = now_s() - tb;
        vals[rp] = 120000.0 / (dt > 1e-9 ? dt : 1e-9);
    }
    std::sort(vals.begin(), vals.end()); double med = vals[2];
    for (int rp = 0; rp < 5; ++rp) dev[rp] = std::fabs(vals[rp] - med);
    std::sort(dev.begin(), dev.end());
    FILE *fb = fopen((std::string(out_dir) + "/s35_bench.csv").c_str(), "w");
    fprintf(fb, "person_days,pd_per_sec_med,pd_per_sec_mad\n120000,%.2f,%.2f\n", med, dev[2]); fclose(fb);

    printf("[cpp] tahap 36/39 (s36)\n[cpp] tahap 37/39 (s37)\n[cpp] tahap 38/39 (s38)\n");
    t0 = now_s(); write_s38_svgs(s_no, s_main, out_dir); tot += now_s() - t0;

    printf("[cpp] tahap 39/39 (s39)\n");
    FILE *fj = fopen((std::string(out_dir) + "/s39_summary.json").c_str(), "w");
    fprintf(fj, "{\"lang\":\"cpp\",\"N\":%lld,\"days\":%lld,\"N_sw\":1000,\"toolchain\":\"Apple clang 16.0.0 (-std=c++20)\",\"total_sec\":%.3f}\n",
            static_cast<long long>(N_main), static_cast<long long>(cfg.main_days), tot); fclose(fj);
    FILE *mf = fopen((std::string(out_dir) + "/files_manifest.csv").c_str(), "w"); fprintf(mf, "file,lines,bytes,fnv1a\n");
    scan_dir(out_dir, "", mf); fclose(mf);
    return 0;
}
