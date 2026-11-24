#!/usr/bin/env python3
"""Python independent brute-force reference, SAW validator, critical regressions & plots (day006b)."""

from __future__ import annotations

import argparse
import json
import math
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

A001411_18 = [
    4, 12, 36, 100, 284, 780, 2172, 5916, 16268,
    44100, 120292, 324932, 881500, 2374444, 6416596,
    17245332, 46466676, 124658732,
]


def brute_force_saw(nmax: int) -> list[int]:
    counts = [0] * (nmax + 1)
    counts[0] = 1
    if nmax <= 0:
        return counts
    vis = {(0, 0)}

    def dfs(x: int, y: int, s: int) -> None:
        ns = s + 1
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            pt = (x + dx, y + dy)
            if pt not in vis:
                counts[ns] += 1
                if ns < nmax:
                    vis.add(pt)
                    dfs(pt[0], pt[1], ns)
                    vis.remove(pt)

    dfs(0, 0, 0)
    return counts


def validate_saw(walk: list[list[int]], expected_n: int) -> bool:
    if len(walk) != expected_n + 1 or tuple(walk[0]) != (0, 0):
        return False
    seen = {(0, 0)}
    for i in range(1, len(walk)):
        px, py = walk[i - 1]
        cx, cy = walk[i]
        if abs(cx - px) + abs(cy - py) != 1 or (cx, cy) in seen:
            return False
        seen.add((cx, cy))
    return True


def ols_fit(x: list[float], y: list[float]) -> tuple[float, float, float, float, float]:
    xa, ya = np.asarray(x, float), np.asarray(y, float)
    X = np.column_stack([np.ones_like(xa), xa])
    inv = np.linalg.inv(X.T @ X)
    b = inv @ (X.T @ ya)
    res = ya - X @ b
    cov = float(np.dot(res, res) / (xa.size - 2)) * inv
    return float(b[0]), float(b[1]), math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1]), float(cov[0, 1])


def fit_ratio(counts: list[int], n_lo: int, n_hi: int) -> dict[str, float]:
    ns = list(range(n_lo, n_hi + 1))
    inv_n = [1.0 / n for n in ns]
    rn = [counts[n] / counts[n - 1] for n in ns]
    mu, sl, se_mu, se_sl, cov01 = ols_fit(inv_n, rn)
    gamma = 1.0 + sl / mu
    g0, g1 = -sl / (mu * mu), 1.0 / mu
    se_g = math.sqrt(max(g0 * g0 * se_mu * se_mu + 2 * g0 * g1 * cov01 + g1 * g1 * se_sl * se_sl, 0.0))
    return {"mu_hat": mu, "mu_se": se_mu, "gamma_hat": gamma, "gamma_se": se_g, "slope": sl}


def fit_scaling(ns: list[int], r2s: list[float]) -> dict[str, float]:
    b0, b1, se0, se1, _ = ols_fit([math.log(n) for n in ns], [math.log(r) for r in r2s])
    return {"nu_hat": 0.5 * b1, "nu_se": 0.5 * se1, "amp_hat": math.exp(b0), "amp_se": math.exp(b0) * se0}


def save_both(fig: plt.Figure, fig_dir: Path, stem: str) -> list[str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    paths = [fig_dir / f"{stem}.png", fig_dir / f"{stem}.svg"]
    fig.savefig(paths[0], dpi=150, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return [p.as_posix() for p in paths]


def make_figures(data: dict, rfit: dict, sfit: dict, fig_dir: Path, prefix: str = "day006b") -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    walks = data["samples_n300"]["samples"][:100]
    dists = [math.hypot(w[-1][0], w[-1][1]) for w in walks]
    dmin, dmax = min(dists), max(dists)

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    for i in sorted(range(len(walks)), key=lambda k: dists[k]):
        xs, ys = [p[0] for p in walks[i]], [p[1] for p in walks[i]]
        ax.plot(xs, ys, color=plt.cm.viridis((dists[i] - dmin) / max(dmax - dmin, 1e-9)), alpha=0.38, lw=0.85)
    ax.scatter([0], [0], color="#d62728", s=30, zorder=5)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(dmin, dmax))
    fig.colorbar(sm, ax=ax, shrink=0.85, label="Jarak ujung R")
    ax.set_aspect("equal")
    ax.set_title("100 Contoh SAW N = 300 (Pivot Madras-Sokal)")
    out["saw_samples"] = save_both(fig, fig_dir, f"{prefix}_saw_samples")

    ns = np.array([p["n"] for p in data["pivots"]], float)
    r2s = np.array([p["r2_mean"] for p in data["pivots"]], float)
    ng = np.geomspace(ns.min() * 0.9, ns.max() * 1.1, 80)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(ns, r2s, "o", label="Pivot MC <R^2>")
    ax.plot(ng, sfit["amp_hat"] * ng ** (2 * sfit["nu_hat"]), "-", label=f"Regresi nu={sfit['nu_hat']:.4f}")
    ax.plot(ng, (r2s[1] / ns[1] ** 1.5) * ng ** 1.5, "--", label="Acuan kemiringan 3/2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ax.set_ylabel("<R^2(N)>")
    ax.set_title("Penskalaan Log-Log <R^2> vs N")
    ax.legend()
    out["saw_scaling"] = save_both(fig, fig_dir, f"{prefix}_saw_scaling")

    cnts = data["go_counts"]
    nmax = data["nmax"]
    n_win = np.arange(nmax - 7, nmax + 1, dtype=float)
    inv_w = 1.0 / n_win
    rw = np.array([cnts[int(n)] / cnts[int(n) - 1] for n in n_win])
    xl = np.linspace(0.0, inv_w.max() * 1.1, 60)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.scatter(inv_w, rw, color="#1f77b4", label=f"r_n (n={nmax-7}..{nmax})")
    ax.plot(xl, rfit["mu_hat"] + rfit["slope"] * xl, "-r", label=f"mu_hat={rfit['mu_hat']:.5f}")
    ax.annotate(f"mu_hat={rfit['mu_hat']:.4f}+-{rfit['mu_se']:.4f}\ngamma_hat={rfit['gamma_hat']:.4f}+-{rfit['gamma_se']:.4f}",
                xy=(0.0, rfit["mu_hat"]), xytext=(0.02, rfit["mu_hat"] - 0.01))
    ax.set_xlabel("1/n")
    ax.set_ylabel("r_n = c_n / c_(n-1)")
    ax.set_title("Metode Rasio SAW: Ekstrapolasi mu & gamma")
    ax.legend(loc="lower right")
    out["saw_ratio"] = save_both(fig, fig_dir, f"{prefix}_saw_ratio")

    if data["rust_single"].get("heatmap"):
        hm = np.array(data["rust_single"]["heatmap"], float)
        hn = int(data["rust_single"]["heatmap_n"])
        hm_log = np.full_like(hm, np.nan)
        mask = hm > 0
        hm_log[mask] = np.log10(hm[mask])
        fig, ax = plt.subplots(figsize=(5.8, 5.4))
        im = ax.imshow(hm_log, origin="lower", extent=[-hn - 0.5, hn + 0.5, -hn - 0.5, hn + 0.5], cmap="magma")
        fig.colorbar(im, ax=ax, shrink=0.85, label="log10(kerapatan ujung)")
        ax.set_title(f"Peta Panas Titik Ujung SAW Eksak n = {hn} (Simetri D4)")
        out["saw_endpoint_heatmap"] = save_both(fig, fig_dir, f"{prefix}_saw_endpoint_heatmap")
    return out


def check_figures(figs: dict[str, list[str]]) -> bool:
    for paths in figs.values():
        for p_str in paths:
            p = Path(p_str)
            if not p.exists() or p.stat().st_size == 0:
                return False
            if p.suffix == ".svg":
                ET.parse(p)
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--go-json", type=Path, required=False)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, default=Path("outputs/figures"))
    ap.add_argument("--t0-file", type=Path, required=True)
    ap.add_argument("--stamp-final", action="store_true")
    args = ap.parse_args()

    t0 = int(args.t0_file.read_text().strip())
    wita = timezone(timedelta(hours=8))
    fmt_t = lambda ep: datetime.fromtimestamp(ep, tz=wita).strftime("%Y-%m-%d %H:%M:%S WITA")

    if args.stamp_final:
        doc = json.loads(args.out_json.read_text())
        now = int(time.time())
        el = now - t0
        doc["timing"]["finished_epoch"] = now
        doc["timing"]["finished_wita"] = fmt_t(now)
        doc["timing"]["elapsed_seconds"] = el
        doc["criteria"]["c10_ten_minute_wall_clock"]["status"] = "tercapai" if el <= 600 else "belum tercapai"
        doc["criteria"]["c10_ten_minute_wall_clock"]["measured"] = {"t0": t0, "finished": now, "elapsed_sec": el}
        args.out_json.write_text(json.dumps(doc, indent=2))
        return

    data = json.loads(args.go_json.read_text())
    nmax = int(data["nmax"])
    bf12 = brute_force_saw(12)
    rc = [int(x) for x in data["rust_single"]["counts"]]
    gc = [int(x) for x in data["go_counts"]]
    mism = [{"n": n, "got": gc[n], "ref": A001411_18[n - 1]} for n in range(1, 19) if gc[n] != A001411_18[n - 1]]
    rfit = fit_ratio(gc, nmax - 7, nmax)
    sfit = fit_scaling([p["n"] for p in data["pivots"]], [p["r2_mean"] for p in data["pivots"]])
    walks = data["samples_n300"]["samples"]
    valid_all = len(walks) == 100 and all(validate_saw(w, 300) for w in walks)
    figs = make_figures(data, rfit, sfit, args.fig_dir)
    figs_ok = check_figures(figs)
    now = int(time.time())
    el = now - t0

    st = lambda ok: "tercapai" if ok else "belum tercapai"
    doc = {
        "session": "day006b", "date": "2026-11-24", "seed": int(data["seed"]),
        "run_completed_without_errors": True,
        "timing": {"t0_epoch": t0, "t0_wita": fmt_t(t0), "finished_epoch": now, "finished_wita": fmt_t(now), "elapsed_seconds": el},
        "enumeration": {"nmax": nmax, "prefix_k": data["prefix_k"], "num_tasks": data["num_tasks"], "counts": gc, "py_bf_n12": bf12, "mismatches": mism},
        "exponents": {**rfit, **sfit},
        "pivots": [{k: p[k] for k in ("n", "seed", "r2_mean", "r2_se", "acc")} for p in data["pivots"]],
        "figures": figs,
        "criteria": {
            "c1_bruteforce_matches_rust_n12": {"status": st(bf12 == rc[:13] == gc[:13]), "measured": gc[12], "threshold": "persis n<=12 (c12=324932)"},
            "c2_go_matches_rust_single_all_n": {"status": st(rc == gc), "measured": gc[nmax], "threshold": f"Go == Rust n<= {nmax}"},
            "c3_memory_sequence_a001411": {"status": st(len(mism) == 0), "measured": len(mism), "threshold": "0 ketidakcocokan n=1..18"},
            "c4_mu_hat": {"status": st(abs(rfit["mu_hat"] - 2.638) < 0.03), "measured": {"val": round(rfit["mu_hat"], 6), "se": round(rfit["mu_se"], 6)}, "threshold": "|mu_hat - 2.638| < 0.03"},
            "c5_gamma_hat": {"status": st(abs(rfit["gamma_hat"] - 1.34375) < 0.08), "measured": {"val": round(rfit["gamma_hat"], 6), "se": round(rfit["gamma_se"], 6)}, "threshold": "|gamma_hat - 1.34375| < 0.08"},
            "c6_nu_hat": {"status": st(abs(sfit["nu_hat"] - 0.75) < 0.03), "measured": {"val": round(sfit["nu_hat"], 6), "se": round(sfit["nu_se"], 6)}, "threshold": "|nu_hat - 0.75| < 0.03"},
            "c7_sample_walks_valid": {"status": st(valid_all), "measured": len(walks), "threshold": "100 jalan N=300 valid (0 tabrakan)"},
            "c8_figures_valid": {"status": st(figs_ok), "measured": len(figs), "threshold": "PNG+SVG ada, tak kosong, XML valid"},
            "c9_worker_invariance": {"status": st(bool(data["worker_invariant"])), "measured": data["workers_n"], "threshold": "1 vs N worker identik"},
            "c10_ten_minute_wall_clock": {"status": st(el <= 600), "measured": {"t0": t0, "finished": now, "elapsed_sec": el}, "threshold": "<= 600 detik dari T0 dl.sh"},
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(doc, indent=2))


if __name__ == "__main__":
    main()
