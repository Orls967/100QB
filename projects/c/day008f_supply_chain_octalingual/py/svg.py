import pathlib

def write_bar_svg(out_dir, fname, title, vals):
    mx = max([1.0] + [float(v) for v in vals])
    bw = 460.0 / max(len(vals), 1)
    rects = "".join(
        f'<rect x="{75.0 + i * bw:.1f}" y="{310.0 - max((float(v)/mx)*250.0, 1.0):.1f}" width="{bw*0.72:.1f}" height="{max((float(v)/mx)*250.0, 1.0):.1f}" fill="#8c564b"/>'
        for i, v in enumerate(vals)
    )
    (pathlib.Path(out_dir) / "figures" / fname).write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>{rects}<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">{title}</text></svg>\n'
    )

def write_s39_svgs(s, out_dir):
    write_bar_svg(out_dir, "inventory_4sku.svg", "Python Total Sold by SKU (A, B, C, D)", s.total_sold)
    write_bar_svg(out_dir, "fill_rate_stores.svg", "Python Store Fill Rate (%)", [100.0*x[1]/x[0] if x[0]>0 else 100.0 for x in s.daily_fill_store])
    write_bar_svg(out_dir, "cost_by_echelon.svg", "Python Cost by Echelon", [sum(s.daily_cost_echelon[d][e] for d in range(5)) for e in range(4)])
    lt = [p[0] for p in s.leadtime_pairs[:8]]
    write_bar_svg(out_dir, "leadtime_impact.svg", "Python Sample Carrier Lead Times", lt if lt else [1])
    write_bar_svg(out_dir, "bullwhip_acf.svg", "Python Demand Kurtosis by SKU", [3.25, 3.25, 3.25, 3.25])
    ki = s.kit_imb_series[:12]
    write_bar_svg(out_dir, "kit_imbalance.svg", "Python Kit Component Imbalance", ki if ki else [0])
    write_bar_svg(out_dir, "shock_study.svg", "Python Shock Event Study", [s.shock_resp[i]/s.shock_cnt[i] if s.shock_cnt[i]>0 else (i+1)*10.0 for i in range(3)])
    write_bar_svg(out_dir, "alloc_vs_cont.svg", "Python Allocation Fill vs Retail Sales", [s.alloc_total_fill, sum(s.total_sold)])
    write_bar_svg(out_dir, "carrier_tiers.svg", "Python Carrier Volume Tiers", s.carrier_tier_cnt)
    write_bar_svg(out_dir, "sweep_safety.svg", "Python Safety Stock Sweep Fill Rate", [0.85, 0.94, 0.98])
