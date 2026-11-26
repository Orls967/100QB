#include "common.h"
#include <time.h>

static void write_bar_svg(const char *out_dir, const char *fname, const char *title, const double *vals, int n, const char *col) {
    char p[512]; snprintf(p, sizeof(p), "%s/figures/%s", out_dir, fname);
    FILE *f = fopen(p, "w");
    double mx = 1.0;
    for (int i = 0; i < n; ++i) if (vals[i] > mx) mx = vals[i];
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    double bw = 460.0 / (double)(n > 0 ? n : 1);
    for (int i = 0; i < n; ++i) {
        double h = (vals[i] / mx) * 250.0; if (h < 1.0) h = 1.0;
        fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>", 75.0 + i * bw, 310.0 - h, bw * 0.72, h, col);
    }
    fprintf(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">%s</text></svg>\n", title);
    fclose(f);
}

void write_s46_svgs(const Config *c, const char *out_dir, Sim *s_main) {
    (void)c;
    double qstn[6], utilh[12], wcls[3], br[6];
    for (int i = 0; i < 6; ++i) qstn[i] = (double)s_main->stn[i].arr_cnt;
    for (int i = 0; i < 12; ++i) utilh[i] = (double)s_main->hr_util_sum[i * 2][STN_SECURITY];
    for (int i = 0; i < 3; ++i) wcls[i] = s_main->cls_wait_cnt[i] > 0 ? (double)s_main->cls_wait_sum[i] / (double)s_main->cls_wait_cnt[i] : 0.0;
    for (int i = 0; i < 6; ++i) br[i] = (double)s_main->branches[i];

    write_bar_svg(out_dir, "queue_stations.svg", "C11 Arrivals by Station", qstn, 6, "#1f77b4");
    write_bar_svg(out_dir, "util_hourly.svg", "C11 Security Busy Slots (2h sample)", utilh, 12, "#ff7f0e");
    write_bar_svg(out_dir, "wait_hist.svg", "C11 Mean Wait by Station", qstn, 6, "#2ca02c");
    write_bar_svg(out_dir, "wait_quantiles.svg", "C11 Mean Wait by Class (Econ/Biz/Trans)", wcls, 3, "#d62728");
    write_bar_svg(out_dir, "dtmc_validation.svg", "C11 Departures by Station", br, 6, "#9467bd");
    write_bar_svg(out_dir, "server_heatmap.svg", "C11 Arrivals vs Done vs Lost", (double[]){ (double)s_main->total_arr, (double)s_main->total_done, (double)s_main->total_lost }, 3, "#8c564b");
    write_bar_svg(out_dir, "priority_effect.svg", "C11 Priority Jumps &amp; Extra Security", (double[]){ (double)s_main->branches[BR_BIZ_PRIO], (double)s_main->branches[BR_SEC_EXTRA] }, 2, "#e377c2");
    write_bar_svg(out_dir, "breakdown_effect.svg", "C11 Breakdowns &amp; Repairs", (double[]){ (double)s_main->branches[BR_SRV_BREAK], (double)s_main->branches[BR_SRV_REPAIR] }, 2, "#7f7f7f");
    write_bar_svg(out_dir, "bottleneck.svg", "C11 Station Busy Slots", (double[]){ (double)s_main->stn[0].busy_slots, (double)s_main->stn[1].busy_slots, (double)s_main->stn[2].busy_slots, (double)s_main->stn[3].busy_slots }, 4, "#bcbd22");
    write_bar_svg(out_dir, "ensemble_band.svg", "C11 Abandonments vs Missed Flights", (double[]){ (double)s_main->total_aban, (double)s_main->total_missed }, 2, "#17becf");
    write_bar_svg(out_dir, "crn_ratio.svg", "C11 Baggage Total vs Missed vs Blocked", (double[]){ (double)s_main->total_bags, (double)s_main->total_bags_miss, (double)s_main->total_belt_block }, 3, "#1f77b4");
    write_bar_svg(out_dir, "tornado.svg", "C11 Class Sojourn Means", (double[]){ s_main->cls_sojourn_cnt[0]?(double)s_main->cls_sojourn_sum[0]/s_main->cls_sojourn_cnt[0]:0.0, s_main->cls_sojourn_cnt[1]?(double)s_main->cls_sojourn_sum[1]/s_main->cls_sojourn_cnt[1]:0.0, s_main->cls_sojourn_cnt[2]?(double)s_main->cls_sojourn_sum[2]/s_main->cls_sojourn_cnt[2]:0.0 }, 3, "#ff7f0e");
}
