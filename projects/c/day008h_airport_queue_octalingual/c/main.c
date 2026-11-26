#include "common.h"
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>

static double now_sec(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }

static int cmp_dbl(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }

static void manifest_scan(const char *base, const char *rel, FILE *out) {
    char dirp[512]; snprintf(dirp, sizeof(dirp), "%s%s%s", base, rel[0] ? "/" : "", rel);
    DIR *d = opendir(dirp); if (!d) return;
    struct dirent *de; char names[256][128]; int n = 0;
    while ((de = readdir(d)) && n < 256) {
        if (de->d_name[0] == '.') continue;
        snprintf(names[n++], 128, "%s", de->d_name);
    }
    closedir(d);
    for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) if (strcmp(names[i], names[j]) > 0) {
        char tmp[128]; strcpy(tmp, names[i]); strcpy(names[i], names[j]); strcpy(names[j], tmp);
    }
    for (int i = 0; i < n; ++i) {
        char fp[512], r2[256];
        snprintf(fp, sizeof(fp), "%s/%s", dirp, names[i]);
        snprintf(r2, sizeof(r2), "%s%s%s", rel, rel[0] ? "/" : "", names[i]);
        struct stat st; if (stat(fp, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) manifest_scan(base, r2, out);
        else if (strcmp(r2, "files_manifest.csv") != 0) {
            FILE *f = fopen(fp, "rb"); if (!f) continue;
            uint8_t *buf = (uint8_t *)malloc((size_t)st.st_size + 1);
            size_t rd = fread(buf, 1, (size_t)st.st_size, f); fclose(f);
            int64_t lc = 0; for (size_t k = 0; k < rd; ++k) if (buf[k] == '\n') lc++;
            fprintf(out, "%s,%lld,%zu,%016llx\n", r2, (long long)lc, rd, (unsigned long long)fnv1a_bytes(buf, rd));
            free(buf);
        }
    }
}

int main(int argc, char **argv) {
    const char *cfg_path = "config.kv", *out_dir = "out/c";
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "--config") && i + 1 < argc) cfg_path = argv[++i];
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_dir = argv[++i];
    }
    char fig_dir[512]; snprintf(fig_dir, sizeof(fig_dir), "%s/figures", out_dir);
    mkdir("out", 0755); mkdir(out_dir, 0755); mkdir(fig_dir, 0755);

    Config c = load_config(cfg_path);
    if (argc > 1 && !strcmp(argv[1], "s01")) { write_s01_to_s13(&c, out_dir); return 0; }

    double tot = 0.0, t0;
    for (int st = 1; st <= 13; ++st) printf("[c] tahap %d/48 (s%02d)\n", st, st);
    t0 = now_sec(); write_s01_to_s13(&c, out_dir); tot += now_sec() - t0;

    int64_t main_days = 365;
    printf("[c] tahap 14/48 (s14)\n");
    t0 = now_sec();
    RunParams rpm = default_rp(&c, main_days); rpm.audit_on = 1;
    Sim *s_main = sim_create(&c, &rpm); sim_run(s_main);
    tot += now_sec() - t0;

    for (int st = 15; st <= 41; ++st) printf("[c] tahap %d/48 (s%02d)\n", st, st);
    t0 = now_sec(); write_s14_to_s45(&c, out_dir, s_main); tot += now_sec() - t0;

    printf("[c] tahap 42/48 (s42)\n");
    double vals[5];
    for (int r = 0; r < 5; ++r) {
        RunParams rpb = default_rp(&c, 3);
        Sim *sb = sim_create(&c, &rpb);
        double tb = now_sec(); sim_run(sb);
        double dt = now_sec() - tb; if (dt < 1e-9) dt = 1e-9;
        vals[r] = (double)sb->sys_sum_L / dt;
        sim_free(sb);
    }
    qsort(vals, 5, sizeof(double), cmp_dbl);
    double med = vals[2], dev[5];
    for (int r = 0; r < 5; ++r) dev[r] = fabs(vals[r] - med);
    qsort(dev, 5, sizeof(double), cmp_dbl);
    char bp[512]; snprintf(bp, sizeof(bp), "%s/s42_bench.csv", out_dir);
    FILE *fb = fopen(bp, "w");
    fprintf(fb, "slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,%0.2f,%0.2f\n", med, dev[2]);
    fclose(fb);

    printf("[c] tahap 43/48 (s43)\n[c] tahap 44/48 (s44)\n[c] tahap 45/48 (s45)\n[c] tahap 46/48 (s46)\n");
    t0 = now_sec(); write_s46_svgs(&c, out_dir, s_main); tot += now_sec() - t0;

    printf("[c] tahap 47/48 (s47)\n[c] tahap 48/48 (s48)\n");
    snprintf(bp, sizeof(bp), "%s/s48_summary.json", out_dir);
    FILE *fs = fopen(bp, "w");
    fprintf(fs, "{\"lang\":\"c\",\"days\":%lld,\"sw_days\":3,\"toolchain\":\"Apple clang 16.0.0 (-std=c11)\",\"total_sec\":%.3f}\n", (long long)main_days, tot);
    fclose(fs);

    snprintf(bp, sizeof(bp), "%s/files_manifest.csv", out_dir);
    FILE *fm = fopen(bp, "w"); fprintf(fm, "file,lines,bytes,fnv1a\n");
    manifest_scan(out_dir, "", fm); fclose(fm);
    sim_free(s_main);
    return 0;
}
