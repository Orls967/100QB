#include "common.hpp"
#include <chrono>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>

static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

static void manifest_scan(const char *dir, const char *prefix, FILE *mf) {
    DIR *d = opendir(dir); if (!d) return;
    struct dirent *de;
    while ((de = readdir(d)) != nullptr) {
        if (de->d_name[0] == '.') continue;
        char p[512]; snprintf(p, sizeof(p), "%s/%s", dir, de->d_name);
        struct stat st; if (stat(p, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) {
            char np[256]; snprintf(np, sizeof(np), "%s%s/", prefix, de->d_name); manifest_scan(p, np, mf);
        } else if (strcmp(de->d_name, "files_manifest.csv") != 0) {
            FILE *fp = fopen(p, "rb"); if (!fp) continue;
            uint64_t h = 0xCBF29CE484222325ULL; long lines = 0, bytes = 0; int ch;
            while ((ch = fgetc(fp)) != EOF) { bytes++; if (ch == '\n') lines++; h ^= static_cast<uint8_t>(ch); h *= 0x100000001B3ULL; }
            fclose(fp);
            fprintf(mf, "%s%s,%ld,%ld,%016llx\n", prefix, de->d_name, lines, bytes, static_cast<unsigned long long>(h));
        }
    }
    closedir(d);
}

void write_s40_manifest(const char *out_dir, const char *toolchain, const double *stage_s, int64_t ev_total) {
    double tot = 0.0; for (int i = 0; i < 40; ++i) tot += stage_s[i];
    char p[512]; snprintf(p, sizeof(p), "%s/s40_summary.json", out_dir);
    FILE *f = fopen(p, "w");
    fprintf(f, "{\"lang\":\"cpp\",\"days\":5,\"ev_per_day\":400000,\"ev_total\":%lld,\"toolchain\":\"%s\",\"total_sec\":%.3f}\n", static_cast<long long>(ev_total), toolchain, tot);
    fclose(f);
    snprintf(p, sizeof(p), "%s/files_manifest.csv", out_dir);
    FILE *mf = fopen(p, "w"); fprintf(mf, "file,lines,bytes,fnv1a\n"); manifest_scan(out_dir, "", mf); fclose(mf);
}

int main(int argc, char **argv) {
    const char *out_dir = "out/cpp";
    mkdir("out", 0755); mkdir(out_dir, 0755); mkdir("out/cpp/figures", 0755);
    Config cfg{}; if (!load_config("config.kv", cfg)) return 1;
    const char *cmd = argc > 1 ? argv[1] : "all";
    if (!strcmp(cmd, "s01")) { write_s01_to_s09(cfg, out_dir); return 0; }
    double stage_s[40]{}; double t0 = now_s();
    for (int st = 1; st <= 9; ++st) printf("[cpp] tahap %d/40 (s%02d)\n", st, st);
    write_s01_to_s09(cfg, out_dir); stage_s[0] = now_s() - t0;
    printf("[cpp] tahap 10/40 (s10)\n"); t0 = now_s();
    Sim sim; sim_init(sim, cfg, 1, out_dir); int64_t ev_pd = 400000;
    sim_run_days(sim, cfg.num_days, ev_pd, 0, -1, 1); stage_s[9] = now_s() - t0;
    for (int st = 11; st <= 38; ++st) printf("[cpp] tahap %d/40 (s%02d)\n", st, st);
    t0 = now_s(); write_s11_to_s38(sim, cfg, out_dir); stage_s[10] = now_s() - t0;
    printf("[cpp] tahap 39/40 (s39)\n"); t0 = now_s(); write_s39_svgs(sim, out_dir); stage_s[38] = now_s() - t0;
    printf("[cpp] tahap 40/40 (s40)\n"); write_s40_manifest(out_dir, "Apple clang 16.0.0 (-std=c++20)", stage_s, cfg.num_days * ev_pd);
    sim_free(sim); return 0;
}
