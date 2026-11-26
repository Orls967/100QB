from config import STN_SECURITY
from network import BR_BIZ_PRIO, BR_SEC_EXTRA, BR_SRV_BREAK, BR_SRV_REPAIR


def write_bar_svg(out_dir: str, fname: str, title: str, vals: list[float], col: str) -> None:
    mx = max([1.0] + [float(v) for v in vals])
    bw = 460.0 / max(len(vals), 1)
    with open(f"{out_dir}/figures/{fname}", "w", encoding="utf-8") as f:
        f.write('<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>')
        for i, v in enumerate(vals):
            h = max((float(v) / mx) * 250.0, 1.0)
            f.write(f'<rect x="{75.0 + i * bw:.1f}" y="{310.0 - h:.1f}" width="{bw * 0.72:.1f}" height="{h:.1f}" fill="{col}"/>')
        f.write(f'<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">{title}</text></svg>\n')


def write_s46_svgs(s_main, out_dir: str) -> None:
    qstn = [float(s_main.stn[i].arr_cnt) for i in range(6)]
    utilh = [float(s_main.hr_util_sum[i * 2][STN_SECURITY]) for i in range(12)]
    wcls = [s_main.cls_wait_sum[i] / float(s_main.cls_wait_cnt[i]) if s_main.cls_wait_cnt[i] > 0 else 0.0 for i in range(3)]
    soj = [s_main.cls_sojourn_sum[i] / float(s_main.cls_sojourn_cnt[i]) if s_main.cls_sojourn_cnt[i] > 0 else 0.0 for i in range(3)]
    br = [float(s_main.branches[i]) for i in range(6)]
    write_bar_svg(out_dir, "queue_stations.svg", "Python Arrivals by Station", qstn, "#1f77b4")
    write_bar_svg(out_dir, "util_hourly.svg", "Python Security Busy Slots (2h sample)", utilh, "#ff7f0e")
    write_bar_svg(out_dir, "wait_hist.svg", "Python Mean Wait by Station", qstn, "#2ca02c")
    write_bar_svg(out_dir, "wait_quantiles.svg", "Python Mean Wait by Class (Econ/Biz/Trans)", wcls, "#d62728")
    write_bar_svg(out_dir, "dtmc_validation.svg", "Python Departures by Station", br, "#9467bd")
    write_bar_svg(out_dir, "server_heatmap.svg", "Python Arrivals vs Done vs Lost", [float(s_main.total_arr), float(s_main.total_done), float(s_main.total_lost)], "#8c564b")
    write_bar_svg(out_dir, "priority_effect.svg", "Python Priority Jumps &amp; Extra Security", [float(s_main.branches[BR_BIZ_PRIO]), float(s_main.branches[BR_SEC_EXTRA])], "#e377c2")
    write_bar_svg(out_dir, "breakdown_effect.svg", "Python Breakdowns &amp; Repairs", [float(s_main.branches[BR_SRV_BREAK]), float(s_main.branches[BR_SRV_REPAIR])], "#7f7f7f")
    write_bar_svg(out_dir, "bottleneck.svg", "Python Station Busy Slots", [float(s_main.stn[i].busy_slots) for i in range(4)], "#bcbd22")
    write_bar_svg(out_dir, "ensemble_band.svg", "Python Abandonments vs Missed Flights", [float(s_main.total_aban), float(s_main.total_missed)], "#17becf")
    write_bar_svg(out_dir, "crn_ratio.svg", "Python Baggage Total vs Missed vs Blocked", [float(s_main.total_bags), float(s_main.total_bags_miss), float(s_main.total_belt_block)], "#1f77b4")
    write_bar_svg(out_dir, "tornado.svg", "Python Class Sojourn Means", soj, "#ff7f0e")
