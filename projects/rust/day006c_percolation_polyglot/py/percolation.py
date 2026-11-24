#!/usr/bin/env python3
"""Exact rational enumeration, Newman-Ziff canonical transform, FSS & plots (day006c)."""
from __future__ import annotations
import argparse, json, math, random, xml.etree.ElementTree as ET
from fractions import Fraction
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[4]
_LG_CACHE: dict[int, np.ndarray] = {}


def build_bonds(rows: int, cols: int) -> list[tuple[int, int]]:
    h = [(r * cols + c, r * cols + c + 1) for r in range(rows) for c in range(cols - 1)]
    return h + [(r * cols + c, (r + 1) * cols + c) for r in range(rows - 1) for c in range(cols)]


def _uf_crosses(nv: int, rows: int, cols: int, open_bonds: list[tuple[int, int]]) -> bool:
    p = list(range(nv + 2))
    vl, vr = nv, nv + 1

    def find(x: int) -> int:
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    for r in range(rows):
        p[find(r * cols)] = find(vl)
        p[find(r * cols + cols - 1)] = find(vr)
    for u, v in open_bonds:
        p[find(u)] = find(v)
    return find(vl) == find(vr)


def exact_enumeration(rows: int, cols: int) -> dict:
    bonds = build_bonds(rows, cols)
    eb, nv = len(bonds), rows * cols
    counts = [0] * (eb + 1)
    for mask in range(1 << eb):
        if _uf_crosses(nv, rows, cols, [bonds[i] for i in range(eb) if (mask >> i) & 1]):
            counts[mask.bit_count()] += 1
    return {"rows": rows, "cols": cols, "n_bonds": eb, "counts": counts, "p_half": Fraction(sum(counts), 1 << eb)}


def eval_exact_poly(counts: list[int], p: Fraction) -> Fraction:
    eb, q = len(counts) - 1, Fraction(1, 1) - p
    return sum((Fraction(counts[k]) * (p ** k) * (q ** (eb - k)) for k in range(eb + 1)), Fraction(0, 1))


def discover_self_dual_shapes() -> dict:
    out = {}
    for r, c in ((2, 3), (3, 2), (2, 2), (3, 4), (4, 3)):
        ex = exact_enumeration(r, c)
        out[f"{r}x{c}"] = {"n_bonds": ex["n_bonds"], "p_half": str(ex["p_half"]), "is_half": ex["p_half"] == Fraction(1, 2), "counts": ex["counts"]}
    return out


def python_newman_ziff(rows: int, cols: int, samples: int, seed: int) -> list[int]:
    bonds = build_bonds(rows, cols)
    eb, nv = len(bonds), rows * cols
    vl, vr, rng, idx, hist = nv, nv + 1, random.Random(seed), list(range(eb)), [0] * (eb + 1)
    for _ in range(samples):
        rng.shuffle(idx)
        p = list(range(nv + 2))

        def find(x: int) -> int:
            while p[x] != x:
                p[x] = p[p[x]]
                x = p[x]
            return x

        for r in range(rows):
            p[find(r * cols)] = find(vl)
            p[find(r * cols + cols - 1)] = find(vr)
        for step, bi in enumerate(idx, 1):
            u, v = bonds[bi]
            p[find(u)] = find(v)
            if find(vl) == find(vr):
                hist[step] += 1
                break
    return hist


def binom_weights(n: int, p: float) -> np.ndarray:
    if p <= 0.0 or p >= 1.0:
        w = np.zeros(n + 1)
        w[0 if p <= 0.0 else n] = 1.0
        return w
    if n not in _LG_CACHE:
        _LG_CACHE[n] = np.array([math.lgamma(k + 1) for k in range(n + 1)])
    lg, ks = _LG_CACHE[n], np.arange(n + 1, dtype=float)
    return np.exp(lg[n] - lg - lg[::-1] + ks * math.log(p) + (n - ks) * math.log(1.0 - p))


def canonical_eval(hist: list[int], ps: list[float]) -> tuple[list[float], list[float]]:
    arr = np.asarray(hist, dtype=float)
    n, q, ks = len(hist) - 1, np.cumsum(arr) / float(arr.sum()), np.arange(len(hist), dtype=float)
    vals, derivs = [], []
    for p in ps:
        w = binom_weights(n, p)
        vals.append(float(np.dot(w, q)))
        derivs.append(float(np.dot(w * ((ks - n * p) / (p * (1.0 - p))), q)) if 0.0 < p < 1.0 else 0.0)
    return vals, derivs


def ols_loglog(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    X = np.column_stack([np.ones_like(lx), lx])
    inv = np.linalg.inv(X.T @ X)
    b = inv @ (X.T @ ly)
    return float(b[0]), float(b[1]), math.sqrt(float(np.dot(ly - X @ b, ly - X @ b) / (lx.size - 2)) * inv[1, 1])


def save_fig(fig: plt.Figure, fig_dir: Path, stem: str) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    ps = [fig_dir / f"{stem}.png", fig_dir / f"{stem}.svg"]
    fig.savefig(ps[0], dpi=140, bbox_inches="tight")
    fig.savefig(ps[1], bbox_inches="tight")
    plt.close(fig)
    return [p.resolve().relative_to(REPO_ROOT).as_posix() if p.resolve().is_relative_to(REPO_ROOT) else p.as_posix() for p in ps]


def check_figures(figs: dict[str, list[str]]) -> bool:
    for paths in figs.values():
        for s in paths:
            p = Path(s) if Path(s).exists() else REPO_ROOT / s
            if not p.exists() or p.stat().st_size == 0:
                return False
            if p.suffix == ".svg":
                ET.parse(p)
    return True


def make_all_figures(sample_doc: dict, curves: list[dict], nu_fit: dict, d_fit: dict, fig_dir: Path, prefix: str = "day006c") -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    rows, cols = int(sample_doc["rows"]), int(sample_doc["cols"])
    nv, parent, sz = rows * cols, list(range(rows * cols)), [1] * (rows * cols)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edges = [(int(u), int(v)) for u, v in sample_doc["open_bonds"]]
    for u, v in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            if sz[ru] < sz[rv]:
                ru, rv = rv, ru
            parent[rv], sz[ru] = ru, sz[ru] + sz[rv]
    top = max(range(nv), key=lambda x: sz[x] if parent[x] == x else -1)
    bg, fg = [], []
    for u, v in edges:
        (fg if find(u) == top else bg).append([(u % cols, u // cols), (v % cols, v // cols)])
    fig, ax = plt.subplots(figsize=(6.0, 5.8))
    ax.add_collection(LineCollection(bg, colors="#94a3b8", linewidths=0.55, alpha=0.65))
    ax.add_collection(LineCollection(fg, colors="#dc2626", linewidths=1.15, alpha=0.95))
    ax.set(xlim=(-1, cols), ylim=(-1, rows), aspect="equal", title=f"Perkolasi Bond L={rows}, p=0.5 (S_max={sz[top]})")
    out["perc_cluster"] = save_fig(fig, fig_dir, f"{prefix}_perc_cluster")

    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    for c in curves:
        ax.plot(c["p_grid"], c["p_cross"], label=f"L={c['L']}")
    ax.axvline(0.5, ls="--", color="k", lw=0.8)
    ax.axhline(0.5, ls=":", color="k", lw=0.8)
    ax.set(xlabel="p", ylabel="P_cross(p)", title="Kurva Seberang P_cross(p)")
    ax.legend()
    out["perc_crossing_curves"] = save_fig(fig, fig_dir, f"{prefix}_perc_crossing_curves")

    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    for c in curves:
        ax.plot([(p - 0.5) * (c["L"] ** 0.75) for p in c["p_grid"]], c["p_cross"], label=f"L={c['L']}")
    ax.set(xlim=(-1.6, 1.6), xlabel="(p - 0.5) * L^(3/4)", ylabel="P_cross", title="Data Collapse (nu = 4/3)")
    ax.legend()
    out["perc_collapse"] = save_fig(fig, fig_dir, f"{prefix}_perc_collapse")

    ls, dps, smax = np.array([c["L"] for c in curves], float), np.array([c["dp_half"] for c in curves], float), np.array([c["mean_largest"] for c in curves], float)
    lg = np.geomspace(ls.min() * 0.9, ls.max() * 1.1, 50)
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.loglog(ls, dps, "o", label="dP/dp(1/2)")
    ax.loglog(lg, (dps[1] / (ls[1] ** 0.75)) * (lg ** 0.75), "--", label="Acuan L^(3/4)")
    ax.set(xlabel="L", ylabel="dP/dp(1/2)", title=f"Skala dP/dp (nu_hat={nu_fit['nu_hat']:.4f})")
    ax.legend()
    out["perc_derivative_scaling"] = save_fig(fig, fig_dir, f"{prefix}_perc_derivative_scaling")

    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.loglog(ls, smax, "s", label="<S_max(N/2)>")
    ax.loglog(lg, (smax[1] / (ls[1] ** (91.0 / 48.0))) * (lg ** (91.0 / 48.0)), "--", label="Acuan D=91/48")
    ax.set(xlabel="L", ylabel="<S_max>", title=f"Skala Klaster Terbesar (D_hat={d_fit['d_hat']:.4f})")
    ax.legend()
    out["perc_largest_cluster"] = save_fig(fig, fig_dir, f"{prefix}_perc_largest_cluster")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--go-json", type=Path, required=True)
    ap.add_argument("--sample-json", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, default=REPO_ROOT / "outputs" / "figures")
    args = ap.parse_args()

    shapes = discover_self_dual_shapes()
    c1_ok = shapes["2x3"]["is_half"] and shapes["3x4"]["is_half"] and (not shapes["3x2"]["is_half"]) and (not shapes["4x3"]["is_half"])
    ex34, nz_py, grid9 = exact_enumeration(3, 4), python_newman_ziff(3, 4, 15000, 20261124), [round(0.1 * k, 1) for k in range(1, 10)]
    nz_vals, _ = canonical_eval(nz_py, grid9)
    max_z_small = max(abs(pv - float(pe := eval_exact_poly(ex34["counts"], Fraction(k, 10)))) / math.sqrt(float(pe) * (1.0 - float(pe)) / 15000.0) for k, pv in zip(range(1, 10), nz_vals))

    go_data, sample_doc = json.loads(args.go_json.read_text()), json.loads(args.sample_json.read_text())
    p_grid = [round(0.35 + 0.005 * i, 4) for i in range(61)]
    curves, lat_rows, max_z_half = [], [], 0.0
    for item in go_data["lattices"]:
        l, m = int(item["L"]), int(item["samples"])
        p_vals, _ = canonical_eval(item["hist"], p_grid)
        (ph,), (dph,) = canonical_eval(item["hist"], [0.5])
        se_h = 0.5 / math.sqrt(m)
        zh = abs(ph - 0.5) / se_h
        max_z_half = max(max_z_half, zh)
        mean_s = float(item["sum_largest"]) / m
        var_s = max(float(item["sum_largest_sq"]) / m - mean_s * mean_s, 0.0)
        curves.append({"L": l, "p_grid": p_grid, "p_cross": p_vals, "dp_half": dph, "mean_largest": mean_s})
        lat_rows.append({"L": l, "n_bonds": item["n_bonds"], "samples": m, "p_half": ph, "se_half": se_h, "z_half": zh, "dp_half": dph, "mean_largest": mean_s, "se_largest": math.sqrt(var_s / m)})

    _, sl_dp, se_dp = ols_loglog([r["L"] for r in lat_rows], [r["dp_half"] for r in lat_rows])
    nu_hat, nu_se = 1.0 / sl_dp, se_dp / (sl_dp * sl_dp)
    _, d_hat, d_se = ols_loglog([r["L"] for r in lat_rows], [r["mean_largest"] for r in lat_rows])
    nu_fit, d_fit = {"slope_inv_nu": sl_dp, "slope_se": se_dp, "nu_hat": nu_hat, "nu_se": nu_se}, {"d_hat": d_hat, "d_se": d_se, "exact_d": 91.0 / 48.0}
    figs = make_all_figures(sample_doc, curves, nu_fit, d_fit, args.fig_dir)
    figs_ok = check_figures(figs)
    st = lambda ok: "tercapai" if ok else "belum tercapai"
    doc = {
        "session": "day006c", "date": "2026-11-24", "seed": int(go_data["seed"]), "run_completed_without_errors": True,
        "self_dual_discovery": {"chosen_orientation": "rows=L, cols=L+1 (left-right crossing)", "shapes": shapes},
        "small_lattice_check": {"shape": "3x4", "samples": 15000, "max_z_9pts": max_z_small},
        "lattices": lat_rows, "exponents": {**nu_fit, **d_fit}, "worker_invariance": {"workers": go_data["workers"], "invariant": bool(go_data["worker_invariant"])}, "figures": figs,
        "criteria": {
            "c1_exact_self_dual_half": {"status": st(c1_ok), "measured": {"2x3": shapes["2x3"]["p_half"], "3x4": shapes["3x4"]["p_half"], "3x2": shapes["3x2"]["p_half"]}, "threshold": "P(1/2) == 1/2 eksak hanya pada L x (L+1)"},
            "c2_nz_vs_exact_small_lattice": {"status": st(max_z_small <= 3.0), "measured": round(max_z_small, 4), "threshold": "maks |Z| <= 3 SE di 9 titik p kisi 3x4"},
            "c3_pcross_half_all_L": {"status": st(max_z_half <= 3.0), "measured": round(max_z_half, 4), "threshold": "|P_cross(1/2) - 0.5| <= 3 SE untuk tiap L"},
            "c4_nu_hat": {"status": st(abs(nu_hat - 4.0 / 3.0) < 0.15), "measured": {"nu_hat": round(nu_hat, 6), "se": round(nu_se, 6), "slope": round(sl_dp, 6)}, "threshold": "|nu_hat - 4/3| < 0.15"},
            "c5_fractal_dim_D_hat": {"status": st(abs(d_hat - 91.0 / 48.0) < 0.06), "measured": {"d_hat": round(d_hat, 6), "se": round(d_se, 6)}, "threshold": "|D_hat - 91/48| < 0.06"},
            "c6_figures_valid": {"status": st(figs_ok and len(figs) == 5), "measured": len(figs), "threshold": "5 pasang PNG+SVG ada, tak kosong, XML valid"},
            "c7_worker_invariance": {"status": st(bool(go_data["worker_invariant"])), "measured": go_data["workers"], "threshold": "1 vs N worker identik bit-demi-bit"},
            "c_wall_clock": {"status": "belum diketahui", "measured": None, "threshold": "elapsed_sec <= 570"},
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
