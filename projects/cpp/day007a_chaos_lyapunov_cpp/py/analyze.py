#!/usr/bin/env python3
"""Bangun outputs/day007a/results.json dari CSV C++ dan validasi SVG (+ PNG opsional)."""
from __future__ import annotations
import csv, json, math, os, platform, subprocess, xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
OUT_DIR = ROOT / "outputs" / "day007a"
FIG_DIR = ROOT / "outputs" / "figures"

metrics = {}
with (OUT_DIR / "summary.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        metrics[row["metric"]] = float(row["value"])

feig_rows = []
with (OUT_DIR / "feigenbaum.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        feig_rows.append({"n": int(row["n"]), "period": int(row["period"]), "R_n": float(row["R_n"]), "delta_n": float(row["delta_n"])})

stems = ["day007a_bifurcation", "day007a_lorenz_attractor", "day007a_lorenz_lyapunov_convergence", "day007a_separation_growth"]
figs = {}
svg_ok = True
for st in stems:
    svg_p = FIG_DIR / f"{st}.svg"
    if not svg_p.exists() or svg_p.stat().st_size == 0:
        svg_ok = False
    else:
        ET.parse(svg_p)
    figs[st] = [svg_p.relative_to(ROOT).as_posix()]

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with (OUT_DIR / "logistic_lyapunov.csv").open(encoding="utf-8") as f:
        rr, ll = zip(*[(float(r["r"]), float(r["lambda"])) for r in csv.DictReader(f)])
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.plot(rr, ll, lw=0.9, color="#1d4ed8")
    ax.axhline(0.0, color="#dc2626", ls="--", lw=0.8)
    ax.set(title="Logistic Map Lyapunov Exponent lambda(r)", xlabel="r", ylabel="lambda(r)")
    fig.tight_layout()
    png_p = FIG_DIR / "day007a_bifurcation.png"
    fig.savefig(png_p, dpi=130)
    plt.close(fig)
    figs["day007a_bifurcation"].append(png_p.relative_to(ROOT).as_posix())
except Exception:
    pass

cxx_ver = subprocess.check_output(["c++", "--version"], text=True).splitlines()[0].strip()
ln2 = math.log(2.0)
c1_err = abs(metrics["lam_r4_1e7"] - ln2)
c2_d_err = abs(metrics["feigenbaum_delta8"] - 4.6692)
c2_r_err = abs(metrics["feigenbaum_r_inf"] - 3.569946)
c3_err = abs(metrics["lorenz_sum_lam"] + 41.0 / 3.0)
c4_l1_err = abs(metrics["lorenz_lam1"] - 0.9056)
c4_l2_abs = abs(metrics["lorenz_lam2"])
c5_err = abs(metrics["lorenz_d_ky"] - 2.06)
r12, r24 = metrics["rk4_ratio_12"], metrics["rk4_ratio_24"]

doc = {
    "session": "day007a",
    "date": "2026-11-25",
    "run_completed_without_errors": True,
    "environment": {"arch": platform.machine(), "cpu_cores": os.cpu_count() or 1, "compiler": cxx_ver, "python": platform.python_version()},
    "measurements": {**metrics, "ln2_exact": ln2, "liouville_sum_exact": -41.0 / 3.0, "feigenbaum_table": feig_rows},
    "figures": figs,
    "criteria": {
        "c1_logistic_r4_ln2": {"status": "tercapai" if c1_err <= 0.005 else "belum tercapai", "measured": {"lambda_r4": round(metrics["lam_r4_1e7"], 6), "abs_err": round(c1_err, 6)}, "threshold": "|lambda_hat(4) - ln 2| <= 0.005 pada 1e7 iterasi"},
        "c2_feigenbaum_delta8_rinf": {"status": "tercapai" if (c2_d_err <= 0.02 and c2_r_err <= 1e-4) else "belum tercapai", "measured": {"delta_8": round(metrics["feigenbaum_delta8"], 6), "r_inf": round(metrics["feigenbaum_r_inf"], 7)}, "threshold": "|delta_8 - 4.6692| <= 0.02 dan |r_inf - 3.569946| <= 1e-4"},
        "c3_lorenz_liouville_sum": {"status": "tercapai" if c3_err <= 0.02 else "belum tercapai", "measured": {"sum_lambda": round(metrics["lorenz_sum_lam"], 6), "abs_err": round(c3_err, 6)}, "threshold": "|sum(lambda_i) + 13.6667| <= 0.02"},
        "c4_lorenz_lam1_lam2": {"status": "tercapai" if (c4_l1_err <= 0.05 and c4_l2_abs <= 0.03) else "belum tercapai", "measured": {"lambda_1": round(metrics["lorenz_lam1"], 5), "lambda_2": round(metrics["lorenz_lam2"], 5)}, "threshold": "|lambda_1 - 0.9056| <= 0.05 dan |lambda_2| <= 0.03"},
        "c5_kaplan_yorke_dim": {"status": "tercapai" if c5_err <= 0.03 else "belum tercapai", "measured": round(metrics["lorenz_d_ky"], 5), "threshold": "|D_KY - 2.06| <= 0.03"},
        "c6_rk4_order4": {"status": "tercapai" if (12.0 <= r12 <= 20.0 and 12.0 <= r24 <= 20.0) else "belum tercapai", "measured": {"ratio_dt_dt2": round(r12, 4), "ratio_dt2_dt4": round(r24, 4)}, "threshold": "rasio galat RK4 dalam [12, 20]"},
        "c7_svg_valid": {"status": "tercapai" if svg_ok else "belum tercapai", "measured": len(stems), "threshold": "berkas SVG ada, tidak kosong, dan valid XML"},
        "c_wall_clock": {"status": "belum diketahui", "measured": None, "threshold": "elapsed_sec <= 570 detik (diisi oleh finalize.py)"},
    },
}
(OUT_DIR / "results.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
print("analyze.py OK -> outputs/day007a/results.json")
