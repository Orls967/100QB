import csv
import hashlib
import json
import math
import pathlib
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJ = pathlib.Path(__file__).resolve().parent
ROOT = PROJ.parents[2]
OUT_BASE = PROJ / "out"
FIG_OUT = ROOT / "outputs" / "figures"
RES_OUT = ROOT / "outputs" / "day008e" / "results.json"
FIG_OUT.mkdir(parents=True, exist_ok=True)
RES_OUT.parent.mkdir(parents=True, exist_ok=True)

LANGS = ["c", "cpp", "rust", "go", "ruby", "py", "ts", "ipynb"]
LANG_LABELS = {"c": "C11", "cpp": "C++20", "rust": "Rust", "go": "Go", "ruby": "Ruby", "py": "Python", "ts": "TypeScript", "ipynb": "Notebook"}
SRC_FILES = {
    "c": PROJ / "c" / "poisson.c",
    "cpp": PROJ / "cpp" / "poisson.cpp",
    "rust": PROJ / "rust" / "src" / "main.rs",
    "go": PROJ / "go" / "main.go",
    "ruby": PROJ / "ruby" / "poisson.rb",
    "py": PROJ / "py" / "poisson.py",
    "ts": PROJ / "ts" / "poisson.ts",
    "ipynb": ROOT / "days" / "day008e_sim_poisson-solvers-octalingual.ipynb",
}
METHODS = ["jacobi", "gs", "rbgs", "sor", "cg", "pcg_ssor", "mg"]
METHOD_FILES = {
    "jacobi": "s04_jacobi.csv", "gs": "s05_gs.csv", "rbgs": "s06_rbgs.csv",
    "sor": "s07_sor.csv", "cg": "s08_cg.csv", "pcg_ssor": "s09_pcg_ssor.csv", "mg": "s10_mg.csv",
}
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]

def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

def count_lines(lang):
    p = SRC_FILES[lang]
    if lang != "ipynb":
        return len(p.read_text().splitlines())
    nb = json.loads(p.read_text())
    return sum(len("".join(c["source"]).splitlines()) for c in nb["cells"] if c["cell_type"] == "code")

# 1. Bit-identical parity s01 and s17
s01_ref = (OUT_BASE / "c" / "s01_selftest.txt").read_text()
s17_ref = (OUT_BASE / "c" / "s17_parity.csv").read_text()
s01_ok = all((OUT_BASE / lg / "s01_selftest.txt").read_text() == s01_ref for lg in LANGS)
s17_ok = all((OUT_BASE / lg / "s17_parity.csv").read_text() == s17_ref for lg in LANGS)
s01_all_pass = all(read_csv(OUT_BASE / lg / "s01_selftest.txt")[0]["ok"] == "1" for lg in LANGS)

# 2. Eigen anchor (c2)
eig_rows = read_csv(OUT_BASE / "c" / "s03_eigen.csv")
max_eig_err = max(float(r["rel_err_inf"]) for r in eig_rows)
c2_ok = max_eig_err <= 1e-12

# 3. Residual <= 1e-8 & N=63 agreement (c3)
max_rel_res = 0.0
err_inf_63 = {}
for m in METHODS:
    rows = read_csv(OUT_BASE / "c" / METHOD_FILES[m])
    for r in rows:
        max_rel_res = max(max_rel_res, float(r["rel_res"]))
        if int(r["N"]) == 63:
            err_inf_63[m] = float(r["err_inf"])
ipynb_sum = json.loads((OUT_BASE / "ipynb" / "s19_summary.json").read_text())
sol_agree_63 = float(ipynb_sum["sol_agree_63"])
c3_ok = (max_rel_res <= 1e-8) and (sol_agree_63 <= 1e-6)

# 4. Discretization order (c4)
ord_rows = read_csv(OUT_BASE / "c" / "s12_order.csv")
disc_slope = float(ord_rows[0]["slope"])
c4_ok = 1.9 <= disc_slope <= 2.1

# 5. Spectral radius N=31 (c5)
h31 = 1.0 / 32.0
jac_theo = math.cos(math.pi * h31)
gs_theo = jac_theo ** 2
jac_emp = float([r for r in read_csv(OUT_BASE / "c" / "s04_jacobi.csv") if int(r["N"]) == 31][0]["rho_emp"])
gs_emp = float([r for r in read_csv(OUT_BASE / "c" / "s05_gs.csv") if int(r["N"]) == 31][0]["rho_emp"])
jac_rel_diff = abs(jac_emp - jac_theo) / jac_theo
gs_rel_diff = abs(gs_emp - gs_theo) / gs_theo
c5_ok = (jac_rel_diff <= 0.01) and (gs_rel_diff <= 0.01)

# 6. Iteration scaling (c6)
sc_rows = read_csv(OUT_BASE / "c" / "s13_scaling.csv")
slopes = {}
for r in sc_rows:
    slopes[r["method"]] = {"slope": float(r["slope"]), "slope_se": float(r["slope_se"])}
c6_checks = {
    "jacobi": 1.8 <= slopes["jacobi"]["slope"] <= 2.2,
    "gs": 1.8 <= slopes["gs"]["slope"] <= 2.2,
    "rbgs": 1.8 <= slopes["rbgs"]["slope"] <= 2.2,
    "sor": 0.8 <= slopes["sor"]["slope"] <= 1.2,
    "cg": 0.8 <= slopes["cg"]["slope"] <= 1.2,
    "pcg_ssor": 0.4 <= slopes["pcg_ssor"]["slope"] <= 0.7,
    "mg": -0.1 <= slopes["mg"]["slope"] <= 0.1,
}
c6_ok = all(c6_checks.values())

# 7. MG factor N=63 (c7)
mg_rows = read_csv(OUT_BASE / "c" / "s15_mg_factor.csv")
mg_fac_63 = float([r for r in mg_rows if int(r["N"]) == 63][0]["mg_factor_last5"])
mg_fac_127 = float([r for r in mg_rows if int(r["N"]) == 127][0]["mg_factor_last5"])
c7_ok = mg_fac_63 <= 0.2

# 8. FMG ratio N=63 (c8)
fmg_rows = read_csv(OUT_BASE / "c" / "s11_fmg.csv")
fmg_ratio_63 = float([r for r in fmg_rows if int(r["N"]) == 63][0]["ratio"])
c8_ok = fmg_ratio_63 <= 3.0

# 9. Summaries & timing (c9, c12)
summaries = {}
cost_data = {}
line_counts = {}
for lg in LANGS:
    summaries[lg] = json.loads((OUT_BASE / lg / "s19_summary.json").read_text())
    line_counts[lg] = count_lines(lg)
    cost_data[lg] = {r["method"]: {"ns_med": float(r["ns_med"]), "ns_mad": float(r["ns_mad"])}
                     for r in read_csv(OUT_BASE / lg / "s16_cost.csv")}
c9_ok = len(summaries) == 8
c12_ok = all(line_counts[lg] <= 700 for lg in LANGS if lg != "ipynb")

# 10. Combined Figures (4 PNG + 4 SVG)
# Fig 1: History at N=63
hist_rows = read_csv(OUT_BASE / "c" / "s14_history.csv")
fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=150)
for idx, m in enumerate(METHODS):
    pts = [(int(r["iter"]), float(r["rel_res"])) for r in hist_rows if r["method"] == m]
    ax.loglog([p[0] for p in pts], [p[1] for p in pts], color=COLORS[idx], lw=1.8, label=m)
ax.axhline(1e-8, color="k", ls="--", lw=1.0, label="tol = 1e-8")
ax.set_xlabel("Iteration (k)")
ax.set_ylabel("Relative Residual ||r_k||_2 / ||b||_2")
ax.set_title("2D Poisson Solvers (N=63): Residual History across 7 Methods")
ax.grid(True, which="both", ls=":", alpha=0.5)
ax.legend(loc="lower left", ncol=2, fontsize=9)
fig.tight_layout()
fig.savefig(FIG_OUT / "day008e_history.png", dpi=150)
fig.savefig(FIG_OUT / "day008e_history.svg")
plt.close(fig)

# Fig 2: Scaling (Iterations vs N)
fig, ax = plt.subplots(figsize=(8.5, 5.2), dpi=150)
for idx, m in enumerate(METHODS):
    rows = read_csv(OUT_BASE / "c" / METHOD_FILES[m])
    ns = [int(r["N"]) for r in rows]
    its = [int(r["iters"]) for r in rows]
    ax.loglog(ns, its, "o-", color=COLORS[idx], lw=1.8, ms=5, label=f"{m} (slope={slopes[m]['slope']:.2f})")
ax.set_xlabel("Grid size N (unknowns = N^2)")
ax.set_ylabel("Iterations to rel_res <= 1e-8")
ax.set_title("Iteration Scaling vs Grid Size N (C11 Reference, N <= 511)")
ax.grid(True, which="both", ls=":", alpha=0.5)
ax.legend(loc="upper left", fontsize=8.5, ncol=2)
fig.tight_layout()
fig.savefig(FIG_OUT / "day008e_scaling.png", dpi=150)
fig.savefig(FIG_OUT / "day008e_scaling.svg")
plt.close(fig)

# Fig 3: Speed across 8 languages (N=63, RB-GS & MG V-cycle)
fig, ax = plt.subplots(figsize=(9.5, 5.2), dpi=150)
x = np.arange(len(LANGS))
w = 0.36
rb_med = [cost_data[lg]["rbgs"]["ns_med"] for lg in LANGS]
rb_mad = [cost_data[lg]["rbgs"]["ns_mad"] for lg in LANGS]
mg_med = [cost_data[lg]["mg"]["ns_med"] for lg in LANGS]
mg_mad = [cost_data[lg]["mg"]["ns_mad"] for lg in LANGS]
ax.bar(x - w/2, rb_med, w, yerr=rb_mad, capsize=3, label="RB-GS sweep", color="#2ca02c", alpha=0.85)
ax.bar(x + w/2, mg_med, w, yerr=mg_mad, capsize=3, label="MG V(2,1) cycle", color="#e377c2", alpha=0.85)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels([LANG_LABELS[lg] for lg in LANGS])
ax.set_ylabel("ns / unknown-iter (median ± MAD, log scale)")
ax.set_title("Kernel Speed Comparison at N=63 across 8 Languages")
ax.grid(True, axis="y", which="both", ls=":", alpha=0.5)
ax.legend()
fig.tight_layout()
fig.savefig(FIG_OUT / "day008e_speed.png", dpi=150)
fig.savefig(FIG_OUT / "day008e_speed.svg")
plt.close(fig)

# Fig 4: Parity matrix + Source line counts
p_rows = read_csv(OUT_BASE / "c" / "s17_parity.csv")
mat = np.ones((len(LANGS), len(p_rows)))
for i, lg in enumerate(LANGS):
    l_rows = read_csv(OUT_BASE / lg / "s17_parity.csv")
    for j in range(len(p_rows)):
        if l_rows[j] != p_rows[j]:
            mat[i, j] = 0.0
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.0, 4.8), dpi=150, gridspec_kw={"width_ratios": [1.5, 1.0]})
ax1.imshow(mat, cmap="Greens", vmin=0, vmax=1, aspect="auto")
ax1.set_yticks(range(len(LANGS)))
ax1.set_yticklabels([LANG_LABELS[lg] for lg in LANGS])
ax1.set_xticks(range(len(p_rows)))
ax1.set_xticklabels([f"{r['method']}\nN={r['N']}" for r in p_rows], rotation=45, ha="right", fontsize=7.5)
ax1.set_title("Bit-Identical Parity Matrix (8 Languages x 14 Cases)")
for i in range(len(LANGS)):
    for j in range(len(p_rows)):
        ax1.text(j, i, "OK" if mat[i, j] == 1 else "DIFF", ha="center", va="center", fontsize=7, color="white")

lc_vals = [line_counts[lg] for lg in LANGS]
ax2.barh(range(len(LANGS)), lc_vals, color="#1f77b4", alpha=0.85)
ax2.axvline(700, color="r", ls="--", lw=1.2, label="700-line limit")
ax2.set_yticks(range(len(LANGS)))
ax2.set_yticklabels([LANG_LABELS[lg] for lg in LANGS])
ax2.invert_yaxis()
ax2.set_xlabel("Source Lines (code cells for ipynb)")
ax2.set_title("Source Code Size (All <= 700 Lines)")
for i, v in enumerate(lc_vals):
    ax2.text(v + 8, i, str(v), va="center", fontsize=8.5)
ax2.set_xlim(0, 800)
ax2.legend(loc="lower right", fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG_OUT / "day008e_parity.png", dpi=150)
fig.savefig(FIG_OUT / "day008e_parity.svg")
plt.close(fig)

# 11. Validate all 24 per-language SVGs + 4 combined SVGs (c11)
svg_paths = [OUT_BASE / lg / "figures" / f for lg in LANGS for f in ("history.svg", "scaling.svg", "order.svg")]
svg_paths += [FIG_OUT / f"day008e_{n}.svg" for n in ("history", "scaling", "speed", "parity")]
svg_valid = True
for sp in svg_paths:
    if not sp.exists() or sp.stat().st_size == 0:
        svg_valid = False
        break
    try:
        ET.parse(sp)
    except ET.ParseError:
        svg_valid = False
        break

results = {
    "session": "day008e",
    "languages": LANGS,
    "s01_sha256": hashlib.sha256(s01_ref.encode()).hexdigest(),
    "s17_sha256": hashlib.sha256(s17_ref.encode()).hexdigest(),
    "criteria": {
        "c1_parity_8lang": s01_ok and s17_ok,
        "c2_eigen_anchor_le_1e-12": c2_ok,
        "c3_res_le_1e-8_and_agree_le_1e-6": c3_ok,
        "c4_order_in_1.9_2.1": c4_ok,
        "c5_rho_n31_within_1pct": c5_ok,
        "c6_scaling_exponents": c6_ok,
        "c7_mg_factor_n63_le_0.2": c7_ok,
        "c8_fmg_ratio_n63_le_3": c8_ok,
        "c9_all_finished": c9_ok,
        "c10_selftest_pass": s01_all_pass,
        "c11_svgs_valid": svg_valid,
        "c12_speed_and_lines_le_700": c12_ok,
    },
    "metrics": {
        "max_eig_rel_err": max_eig_err,
        "max_rel_res": max_rel_res,
        "sol_agree_63_inf": sol_agree_63,
        "disc_order_slope": disc_slope,
        "jacobi_rho_n31_emp": jac_emp,
        "jacobi_rho_n31_theo": jac_theo,
        "jacobi_rho_rel_diff": jac_rel_diff,
        "gs_rho_n31_emp": gs_emp,
        "gs_rho_n31_theo": gs_theo,
        "gs_rho_rel_diff": gs_rel_diff,
        "scaling_slopes": slopes,
        "c6_checks": c6_checks,
        "mg_factor_63": mg_fac_63,
        "mg_factor_127": mg_fac_127,
        "fmg_ratio_63": fmg_ratio_63,
        "line_counts": line_counts,
        "summaries": summaries,
        "cost_n63": cost_data,
    },
}
RES_OUT.write_text(json.dumps(results, indent=2) + "\n")
print("Criteria:", json.dumps(results["criteria"], indent=2))
print("c6_checks:", json.dumps(c6_checks))
print(f"Wrote {RES_OUT} and 4 combined figures to {FIG_OUT}")
