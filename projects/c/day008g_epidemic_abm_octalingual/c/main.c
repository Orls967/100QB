#include "common.h"
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }
static int cmp_d(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }

static void ensure_dir(const char *p) { mkdir(p, 0755); }

static void scan_dir(const char *base, const char *rel, FILE *mf) {
    char path[512];
    if (rel[0]) snprintf(path, sizeof(path), "%s/%s", base, rel);
    else snprintf(path, sizeof(path), "%s", base);
    DIR *d = opendir(path); if (!d) return;
    struct dirent *ent;
    char *names[256]; int n = 0;
    while ((ent = readdir(d)) != NULL && n < 256) {
        if (ent->d_name[0] == '.') continue;
        names[n++] = strdup(ent->d_name);
    }
    closedir(d);
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j)
        if (strcmp(names[i], names[j]) > 0) { char *t = names[i]; names[i] = names[j]; names[j] = t; }
    for (int i = 0; i < n; ++i) {
        char subrel[512], full[512];
        if (rel[0]) snprintf(subrel, sizeof(subrel), "%s/%s", rel, names[i]);
        else snprintf(subrel, sizeof(subrel), "%s", names[i]);
        snprintf(full, sizeof(full), "%s/%s", base, subrel);
        struct stat st;
        if (stat(full, &st) == 0 && S_ISDIR(st.st_mode)) {
            scan_dir(base, subrel, mf);
        } else if (strcmp(names[i], "files_manifest.csv") != 0) {
            FILE *fp = fopen(full, "rb");
            if (fp) {
                fseek(fp, 0, SEEK_END); long sz = ftell(fp); fseek(fp, 0, SEEK_SET);
                uint8_t *buf = (uint8_t *)malloc((size_t)(sz > 0 ? sz : 1));
                size_t rd = sz > 0 ? fread(buf, 1, (size_t)sz, fp) : 0; fclose(fp);
                long lines = 0; for (size_t k = 0; k < rd; ++k) if (buf[k] == '\n') lines++;
                fprintf(mf, "%s,%ld,%zu,%016llx\n", subrel, lines, rd, (unsigned long long)fnv1a_bytes(buf, rd));
                free(buf);
            }
        }
        free(names[i]);
    }
}

int main(int argc, char **argv) {
    const char *cmd = argc > 1 ? argv[1] : "all";
    const char *cfg_path = "config.kv";
    const char *out_dir = "out/c";
    for (int i = 1; i + 1 < argc; ++i) {
        if (!strcmp(argv[i], "--config")) cfg_path = argv[i + 1];
        else if (!strcmp(argv[i], "--out")) out_dir = argv[i + 1];
    }
    Config cfg;
    if (!load_config(cfg_path, &cfg)) { fprintf(stderr, "Cannot open config %s\n", cfg_path); return 1; }
    ensure_dir("out"); ensure_dir(out_dir);
    char fig_dir[512]; snprintf(fig_dir, sizeof(fig_dir), "%s/figures", out_dir); ensure_dir(fig_dir);

    if (!strcmp(cmd, "s01")) { write_s01_to_s09(&cfg, out_dir); return 0; }

    double tot = 0.0;
    for (int st = 1; st <= 9; ++st) printf("[c] tahap %d/39 (s%02d)\n", st, st);
    double t0 = now_s(); write_s01_to_s09(&cfg, out_dir); tot += now_s() - t0;

    int64_t N_main = 100000, beta_cal = 400;
    printf("[c] tahap 10/39 (s10)\n");
    t0 = now_s();
    RunParams rp_no = {N_main, cfg.main_days, beta_cal, -1, -1, -1, -1, -1, -1, -1, -1, 10000, 0, 0};
    Sim s_no; sim_init(&s_no, &cfg, &rp_no, 1); sim_run(&s_no);
    tot += now_s() - t0;

    printf("[c] tahap 11/39 (s11)\n");
    t0 = now_s();
    RunParams rp_main = {N_main, cfg.main_days, beta_cal, -1, -1, -1, -1, -1, -1, -1, -1, 10000, 0, 1};
    Sim s_main; sim_init(&s_main, &cfg, &rp_main, 1); sim_run(&s_main);
    tot += now_s() - t0;

    for (int st = 12; st <= 34; ++st) printf("[c] tahap %d/39 (s%02d)\n", st, st);
    t0 = now_s(); write_s10_to_s37(&s_no, &s_main, &cfg, out_dir, beta_cal); tot += now_s() - t0;

    printf("[c] tahap 35/39 (s35)\n");
    double vals[5], dev[5];
    for (int rp = 0; rp < 5; ++rp) {
        RunParams rpb = {2000, 60, beta_cal, -1, -1, -1, -1, -1, -1, -1, -1, 10000, 0, 1};
        Sim sb; sim_init(&sb, &cfg, &rpb, 0);
        double tb = now_s(); sim_run(&sb);
        double dt = now_s() - tb;
        vals[rp] = (2000.0 * 60.0) / (dt > 1e-9 ? dt : 1e-9);
        sim_free(&sb);
    }
    qsort(vals, 5, sizeof(double), cmp_d);
    double med = vals[2];
    for (int rp = 0; rp < 5; ++rp) dev[rp] = fabs(vals[rp] - med);
    qsort(dev, 5, sizeof(double), cmp_d);
    char p[512]; snprintf(p, sizeof(p), "%s/s35_bench.csv", out_dir);
    FILE *fb = fopen(p, "w");
    fprintf(fb, "person_days,pd_per_sec_med,pd_per_sec_mad\n120000,%.2f,%.2f\n", med, dev[2]);
    fclose(fb);

    printf("[c] tahap 36/39 (s36)\n");
    printf("[c] tahap 37/39 (s37)\n");
    printf("[c] tahap 38/39 (s38)\n");
    t0 = now_s(); write_s38_svgs(&s_no, &s_main, out_dir); tot += now_s() - t0;

    printf("[c] tahap 39/39 (s39)\n");
    snprintf(p, sizeof(p), "%s/s39_summary.json", out_dir);
    FILE *fj = fopen(p, "w");
    fprintf(fj, "{\"lang\":\"c\",\"N\":%lld,\"days\":%lld,\"N_sw\":1000,\"toolchain\":\"Apple clang 16.0.0 (-std=c11)\",\"total_sec\":%.3f}\n",
            (long long)N_main, (long long)cfg.main_days, tot);
    fclose(fj);

    snprintf(p, sizeof(p), "%s/files_manifest.csv", out_dir);
    FILE *mf = fopen(p, "w"); fprintf(mf, "file,lines,bytes,fnv1a\n");
    scan_dir(out_dir, "", mf); fclose(mf);

    sim_free(&s_no); sim_free(&s_main);
    return 0;
}
