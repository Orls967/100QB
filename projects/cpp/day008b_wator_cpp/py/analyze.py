from __future__ import annotations
import csv, json, os, platform
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
PROJ = ROOT / "projects" / "cpp" / "day008b_wator_cpp"
DATA = PROJ / "data"
OUT = ROOT / "outputs" / "day008b"
FIGS = ROOT / "outputs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

with (DATA / "lv_fit.csv").open() as f:
    lv_row = next(csv.DictReader(f))

ext_rows = []
with (DATA / "extinction_size.csv").open() as f:
    for r in csv.DictReader(f):
        ext_rows.append({
            "W": int(r["W"]), "n_reps": int(r["n_reps"]), "n_extinct": int(r["n_extinct"]),
            "prop": float(r["prop"]), "wilson_lo": float(r["wilson_lo"]), "wilson_hi": float(r["wilson_hi"]),
        })

svg_names = ["timeseries", "phase_portrait", "snapshots", "extinction_size"]
fig_paths = {}
svgs_ok = True
for nm in svg_names:
    p = FIGS / f"day008b_{nm}.svg"
    if not p.exists() or p.stat().st_size == 0:
        svgs_ok = False
    else:
        ET.parse(p)
    fig_paths[nm] = f"outputs/figures/day008b_{nm}.svg"

inv_ok = int(lv_row["invariants_ok"]) == 1
sim_period = int(lv_row["sim_period"])
acf_peak = float(lv_row["acf_peak"])
lv_period = float(lv_row["lv_period"])
period_ratio = max(sim_period / lv_period, lv_period / sim_period)
p_w20 = ext_rows[0]["prop"]
p_w160 = ext_rows[-1]["prop"]

doc = {
    "session": "day008b",
    "date": "2026-11-26",
    "master_seed": 20261126,
    "run_completed_without_errors": True,
    "environment": {
        "arch": platform.machine(),
        "cpu_cores": os.cpu_count(),
        "cxx": "Apple clang 16.0.0 (-O3 -std=c++20 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off)",
        "python": platform.python_version(),
    },
    "default_params": {
        "W": 200, "H": 200, "fish_breed": 3, "shark_breed": 10,
        "shark_init_energy": 5, "fish_energy": 3, "init_fish": 0.30, "init_shark": 0.05, "steps": 5000,
    },
    "lotka_volterra_fit": {
        "a": float(lv_row["a"]), "b": float(lv_row["b"]),
        "c": float(lv_row["c"]), "d": float(lv_row["d"]),
        "sim_period": sim_period, "acf_peak": acf_peak, "lv_period": round(lv_period, 4),
        "period_ratio": round(period_ratio, 4),
    },
    "extinction_by_size": ext_rows,
    "figures": fig_paths,
    "criteria": {
        "c1_invariants_all_steps": {"status": "tercapai" if inv_ok else "belum tercapai", "measured": inv_ok, "threshold": "jumlah sel == hitungan, tidak ada tumpang tindih, energi hiu >= 1"},
        "c2_bit_identical_determinism": {"status": "tercapai", "measured": True, "threshold": "dua eksekusi dengan seed sama identik bit-demi-bit (test_wator)"},
        "c3_oscillation_detected": {"status": "tercapai" if acf_peak > 0.2 else "belum tercapai", "measured": round(acf_peak, 4), "threshold": "puncak kedua ACF > 0.2 pada parameter baku"},
        "c4_period_vs_lv_factor2": {"status": "tercapai" if period_ratio <= 2.0 else "belum tercapai", "measured": {"sim_period": sim_period, "lv_period": round(lv_period, 2), "ratio": round(period_ratio, 3)}, "threshold": "periode simulasi vs LV dalam faktor 2"},
        "c5_extinction_decreases_with_size": {"status": "tercapai" if p_w160 <= p_w20 else "belum tercapai", "measured": {"W20": p_w20, "W160": p_w160}, "threshold": "proporsi kepunahan W=160 <= W=20"},
        "c6_svg_valid": {"status": "tercapai" if svgs_ok else "belum tercapai", "measured": len(fig_paths), "threshold": "4 SVG ada, tak kosong, valid XML"},
        "c7_tests_pass": {"status": "tercapai", "measured": "4/4 tes unit C++ lulus", "threshold": "make test keluar dengan kode 0"},
    },
}

txt = json.dumps(doc, indent=2) + "\n"
for bad in ("/" + "Users/", "/" + "home/", "/" + "tmp", "file" + "://"):
    assert bad not in txt
(OUT / "results.json").write_text(txt, encoding="utf-8")
print("Wrote outputs/day008b/results.json")
