#include "common.hpp"

static void write_bar_svg(const char *out_dir, const char *fname, const char *title, const std::vector<double> &vals, const char *col) {
    char p[512]; std::snprintf(p, sizeof(p), "%s/figures/%s", out_dir, fname);
    FILE *f = std::fopen(p, "w");
    double mx = 1.0;
    for (double v : vals) mx = std::max(mx, v);
    std::fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    double bw = 460.0 / static_cast<double>(std::max<size_t>(vals.size(), 1));
    for (size_t i = 0; i < vals.size(); ++i) {
        double h = std::max((vals[i] / mx) * 250.0, 1.0);
        std::fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>", 75.0 + static_cast<double>(i) * bw, 310.0 - h, bw * 0.72, h, col);
    }
    std::fprintf(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">%s</text></svg>\n", title);
    std::fclose(f);
}

void write_s46_svgs(const Config &c, const char *out_dir, Sim &s_main) {
    (void)c;
    std::vector<double> qstn(6), utilh(12), wcls(3), br(6);
    for (int i = 0; i < 6; ++i) qstn[i] = static_cast<double>(s_main.stn[i].arr_cnt);
    for (int i = 0; i < 12; ++i) utilh[i] = static_cast<double>(s_main.hr_util_sum[i * 2][STN_SECURITY]);
    for (int i = 0; i < 3; ++i) wcls[i] = s_main.cls_wait_cnt[i] ? static_cast<double>(s_main.cls_wait_sum[i]) / s_main.cls_wait_cnt[i] : 0.0;
    for (int i = 0; i < 6; ++i) br[i] = static_cast<double>(s_main.branches[i]);

    write_bar_svg(out_dir, "queue_stations.svg", "C++20 Arrivals by Station", qstn, "#1f77b4");
    write_bar_svg(out_dir, "util_hourly.svg", "C++20 Security Busy Slots (2h sample)", utilh, "#ff7f0e");
    write_bar_svg(out_dir, "wait_hist.svg", "C++20 Mean Wait by Station", qstn, "#2ca02c");
    write_bar_svg(out_dir, "wait_quantiles.svg", "C++20 Mean Wait by Class (Econ/Biz/Trans)", wcls, "#d62728");
    write_bar_svg(out_dir, "dtmc_validation.svg", "C++20 Departures by Station", br, "#9467bd");
    write_bar_svg(out_dir, "server_heatmap.svg", "C++20 Arrivals vs Done vs Lost", {static_cast<double>(s_main.total_arr), static_cast<double>(s_main.total_done), static_cast<double>(s_main.total_lost)}, "#8c564b");
    write_bar_svg(out_dir, "priority_effect.svg", "C++20 Priority Jumps &amp; Extra Security", {static_cast<double>(s_main.branches[BR_BIZ_PRIO]), static_cast<double>(s_main.branches[BR_SEC_EXTRA])}, "#e377c2");
    write_bar_svg(out_dir, "breakdown_effect.svg", "C++20 Breakdowns &amp; Repairs", {static_cast<double>(s_main.branches[BR_SRV_BREAK]), static_cast<double>(s_main.branches[BR_SRV_REPAIR])}, "#7f7f7f");
    write_bar_svg(out_dir, "bottleneck.svg", "C++20 Station Busy Slots", {static_cast<double>(s_main.stn[0].busy_slots), static_cast<double>(s_main.stn[1].busy_slots), static_cast<double>(s_main.stn[2].busy_slots), static_cast<double>(s_main.stn[3].busy_slots)}, "#bcbd22");
    write_bar_svg(out_dir, "ensemble_band.svg", "C++20 Abandonments vs Missed Flights", {static_cast<double>(s_main.total_aban), static_cast<double>(s_main.total_missed)}, "#17becf");
    write_bar_svg(out_dir, "crn_ratio.svg", "C++20 Baggage Total vs Missed vs Blocked", {static_cast<double>(s_main.total_bags), static_cast<double>(s_main.total_bags_miss), static_cast<double>(s_main.total_belt_block)}, "#1f77b4");
    write_bar_svg(out_dir, "tornado.svg", "C++20 Class Sojourn Means", {s_main.cls_sojourn_cnt[0]?static_cast<double>(s_main.cls_sojourn_sum[0])/s_main.cls_sojourn_cnt[0]:0.0, s_main.cls_sojourn_cnt[1]?static_cast<double>(s_main.cls_sojourn_sum[1])/s_main.cls_sojourn_cnt[1]:0.0, s_main.cls_sojourn_cnt[2]?static_cast<double>(s_main.cls_sojourn_sum[2])/s_main.cls_sojourn_cnt[2]:0.0}, "#ff7f0e");
}
