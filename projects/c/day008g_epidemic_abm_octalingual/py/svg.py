from population import ST_IA, ST_IS, ST_R, ST_D


def write_bar_svg(out_dir: str, fname: str, title: str, vals, col: str):
    mx = max([1.0] + [float(v) for v in vals])
    bw = 460.0 / max(len(vals), 1)
    with open(f"{out_dir}/figures/{fname}", "w", encoding="utf-8") as f:
        f.write('<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>')
        for i, v in enumerate(vals):
            h = max((float(v) / mx) * 250.0, 1.0)
            f.write(f'<rect x="{75.0 + i * bw:.1f}" y="{310.0 - h:.1f}" width="{bw * 0.72:.1f}" height="{h:.1f}" fill="{col}"/>')
        f.write(f'<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">{title}</text></svg>\n')


def write_s38_svgs(s_no, s_main, out_dir: str):
    d_cnt = s_main.rp.days
    step = max(d_cnt // 12, 1)
    epi, inc, rt, hosp, ens = [], [], [], [], []
    for i in range(12):
        d = min(i * step, d_cnt - 1)
        epi.append(float(s_main.daily_cnt[d][ST_IA] + s_main.daily_cnt[d][ST_IS]))
        inc.append(float(sum(s_main.daily_inc_age[d])))
        rt.append(float(s_no.daily_cnt[d][ST_IA] + s_no.daily_cnt[d][ST_IS]))
        hosp.append(float(s_main.daily_hosp[d][0]))
        ens.append(float(s_main.daily_var[d][1]))
    write_bar_svg(out_dir, "epi_curves.svg", "Python Active Infectious (Interv)", epi, "#1f77b4")
    write_bar_svg(out_dir, "incidence_age.svg", "Python Daily Incidence Sample", inc, "#ff7f0e")
    write_bar_svg(out_dir, "rt_series.svg", "Python Active Infectious (No Interv)", rt, "#d62728")
    write_bar_svg(out_dir, "hospital_cap.svg", "Python Hospital Occupancy Sample", hosp, "#9467bd")
    write_bar_svg(out_dir, "layer_attack.svg", "Python Infections by Layer (HH/SW/Comm)", [float(x) for x in s_main.layer_inf], "#2ca02c")
    write_bar_svg(out_dir, "offspring_dist.svg", "Python Infections by Age Group", [float(x) for x in s_main.age_inf], "#8c564b")
    write_bar_svg(out_dir, "sweep_heatmap.svg", "Python No-Interv vs Interv Total Infections", [float(s_no.total_trans), float(s_main.total_trans)], "#e377c2")
    write_bar_svg(out_dir, "vaccine_effect.svg", "Python Vaccine Doses &amp; Blocked Transmissions", [float(s_main.branches[18]), float(s_main.branches[3])], "#17becf")
    write_bar_svg(out_dir, "variant_share.svg", "Python Variant 1 Daily Incidence", ens, "#bcbd22")
    write_bar_svg(out_dir, "ensemble_band.svg", "Python Final Recovered vs Deaths", [float(s_main.daily_cnt[d_cnt - 1][ST_R]), float(s_main.daily_cnt[d_cnt - 1][ST_D])], "#7f7f7f")
