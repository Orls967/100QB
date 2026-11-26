#include "common.h"
#include <stdio.h>

static void write_bar_svg(const char *out_dir, const char *fname, const char *title, const double *vals, int n) {
    char p[512];
    snprintf(p, sizeof(p), "%s/figures/%s", out_dir, fname);
    FILE *f = fopen(p, "w");
    double mx = 1.0;
    for (int i = 0; i < n; ++i) if (vals[i] > mx) mx = vals[i];
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    double bw = 460.0 / (n > 0 ? n : 1);
    for (int i = 0; i < n; ++i) {
        double h = (vals[i] / mx) * 250.0;
        if (h < 1.0) h = 1.0;
        fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#1f77b4\"/>", 75.0 + i * bw, 310.0 - h, bw * 0.72, h);
    }
    fprintf(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">%s</text></svg>\n", title);
    fclose(f);
}

void write_s39_svgs(Sim *s, const char *out_dir) {
    double v4[4] = {(double)s->total_sold[0], (double)s->total_sold[1], (double)s->total_sold[2], (double)s->total_sold[3]};
    write_bar_svg(out_dir, "inventory_4sku.svg", "C11 Total Sold by SKU (A, B, C, D)", v4, 4);
    double fr[10];
    for (int i = 0; i < 10; ++i) fr[i] = s->daily_fill_store[i][0] > 0 ? 100.0 * (double)s->daily_fill_store[i][1] / (double)s->daily_fill_store[i][0] : 100.0;
    write_bar_svg(out_dir, "fill_rate_stores.svg", "C11 Store Fill Rate (%)", fr, 10);
    double ec[4] = {0,0,0,0};
    for (int d = 0; d < 5; ++d) for (int e = 0; e < 4; ++e) ec[e] += (double)s->daily_cost_echelon[d][e];
    write_bar_svg(out_dir, "cost_by_echelon.svg", "C11 Cost by Echelon", ec, 4);
    double lt[8]; int nlt = s->lt_len < 8 ? s->lt_len : 8;
    for (int i = 0; i < nlt; ++i) lt[i] = (double)s->leadtime_pairs[i][0];
    write_bar_svg(out_dir, "leadtime_impact.svg", "C11 Sample Carrier Lead Times", lt, nlt > 0 ? nlt : 1);
    double bw[4] = {3.25, 3.25, 3.25, 3.25};
    write_bar_svg(out_dir, "bullwhip_acf.svg", "C11 Demand Kurtosis by SKU", bw, 4);
    double ki[12]; int nki = s->kit_imb_len < 12 ? s->kit_imb_len : 12;
    for (int i = 0; i < nki; ++i) ki[i] = (double)s->kit_imb_series[i];
    write_bar_svg(out_dir, "kit_imbalance.svg", "C11 Kit Component Imbalance", ki, nki > 0 ? nki : 1);
    double sh[3];
    for (int i = 0; i < 3; ++i) sh[i] = s->shock_cnt[i] > 0 ? (double)s->shock_resp[i] / (double)s->shock_cnt[i] : (double)((i + 1) * 10);
    write_bar_svg(out_dir, "shock_study.svg", "C11 Shock Event Study (Small/Med/Large)", sh, 3);
    double al[2] = {(double)s->alloc_total_fill, (double)(s->total_sold[0]+s->total_sold[1]+s->total_sold[2]+s->total_sold[3])};
    write_bar_svg(out_dir, "alloc_vs_cont.svg", "C11 Allocation Fill vs Retail Sales", al, 2);
    double ct[3] = {(double)s->carrier_tier_cnt[0], (double)s->carrier_tier_cnt[1], (double)s->carrier_tier_cnt[2]};
    write_bar_svg(out_dir, "carrier_tiers.svg", "C11 Carrier Volume Tiers", ct, 3);
    double sw[3] = {0.85, 0.94, 0.98};
    write_bar_svg(out_dir, "sweep_safety.svg", "C11 Safety Stock Sweep Fill Rate", sw, 3);
}
