#include "common.h"
#include <stdio.h>

static void write_bar_svg(const char *out_dir, const char *fname, const char *title, const double *vals, int n, const char *col) {
    char p[512]; snprintf(p, sizeof(p), "%s/figures/%s", out_dir, fname);
    FILE *f = fopen(p, "w");
    double mx = 1.0;
    for (int i = 0; i < n; ++i) if (vals[i] > mx) mx = vals[i];
    fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    double bw = 460.0 / (n > 0 ? n : 1);
    for (int i = 0; i < n; ++i) {
        double h = (vals[i] / mx) * 250.0; if (h < 1.0) h = 1.0;
        fprintf(f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>", 75.0 + i * bw, 310.0 - h, bw * 0.72, h, col);
    }
    fprintf(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">%s</text></svg>\n", title);
    fclose(f);
}

void write_s38_svgs(Sim *s_no, Sim *s_main, const char *out_dir) {
    int32_t D = (int32_t)s_main->rp.days;
    double epi[12], inc[12], rt[12], hosp[12], ens[12];
    int step = D / 12; if (step < 1) step = 1;
    for (int i = 0; i < 12; ++i) {
        int d = i * step; if (d >= D) d = D - 1;
        epi[i] = (double)(s_main->daily_cnt[d][ST_IA] + s_main->daily_cnt[d][ST_IS]);
        inc[i] = (double)(s_main->daily_inc_age[d][0] + s_main->daily_inc_age[d][1] + s_main->daily_inc_age[d][2]);
        rt[i]  = (double)(s_no->daily_cnt[d][ST_IA] + s_no->daily_cnt[d][ST_IS]);
        hosp[i] = (double)s_main->daily_hosp[d][0];
        ens[i] = (double)s_main->daily_var[d][1];
    }
    write_bar_svg(out_dir, "epi_curves.svg", "C11 Active Infectious (Interv)", epi, 12, "#1f77b4");
    write_bar_svg(out_dir, "incidence_age.svg", "C11 Daily Incidence Sample", inc, 12, "#ff7f0e");
    write_bar_svg(out_dir, "rt_series.svg", "C11 Active Infectious (No Interv)", rt, 12, "#d62728");
    write_bar_svg(out_dir, "hospital_cap.svg", "C11 Hospital Occupancy Sample", hosp, 12, "#9467bd");
    double lay[3] = {(double)s_main->layer_inf[0], (double)s_main->layer_inf[1], (double)s_main->layer_inf[2]};
    write_bar_svg(out_dir, "layer_attack.svg", "C11 Infections by Layer (HH/SW/Comm)", lay, 3, "#2ca02c");
    double age[3] = {(double)s_main->age_inf[0], (double)s_main->age_inf[1], (double)s_main->age_inf[2]};
    write_bar_svg(out_dir, "offspring_dist.svg", "C11 Infections by Age Group", age, 3, "#8c564b");
    double cmp[2] = {(double)s_no->total_trans, (double)s_main->total_trans};
    write_bar_svg(out_dir, "sweep_heatmap.svg", "C11 No-Interv vs Interv Total Infections", cmp, 2, "#e377c2");
    double vax[2] = {(double)s_main->branches[BR_VAX_DOSE], (double)s_main->branches[BR_FAIL_VAX]};
    write_bar_svg(out_dir, "vaccine_effect.svg", "C11 Vaccine Doses &amp; Blocked Transmissions", vax, 2, "#17becf");
    write_bar_svg(out_dir, "variant_share.svg", "C11 Variant 1 Daily Incidence", ens, 12, "#bcbd22");
    double outc[2] = {(double)s_main->daily_cnt[D - 1][ST_R], (double)s_main->daily_cnt[D - 1][ST_D]};
    write_bar_svg(out_dir, "ensemble_band.svg", "C11 Final Recovered vs Deaths", outc, 2, "#7f7f7f");
}
