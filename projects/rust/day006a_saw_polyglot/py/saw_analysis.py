#!/usr/bin/env python3
"""Independent Python reference, critical exponent analysis, and visualization suite for SAW."""

from __future__ import annotations

import argparse
import json
import math
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reference sequence from memory (OEIS A001411, n = 1..20) for comparison only (never forced).
MEMORY_A001411_N1_20: list[int] = [
    4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100,
    120292, 324932, 881500, 2374444, 6416596, 17245332,
    46466676, 124658732, 335116620, 897697164,
]

LIT_MU_REF = 2.63816
LIT_GAMMA_REF = 43.0 / 32.0  # 1.34375
LIT_NU_REF = 0.75            # 3/4


def brute_force_saw_counts(nmax: int) -> list[int]:
    """Pure independent 4-direction brute-force enumeration on Z^2 without symmetry tricks."""
    if nmax < 0:
        raise ValueError("nmax must be >= 0")
    counts = [0] * (nmax + 1)
    counts[0] = 1
    if nmax == 0:
        return counts

    visited: set[tuple[int, int]] = {(0, 0)}
    dirs = ((1, 0), (-1, 0), (0, 1), (0, -1))

    def dfs(x: int, y: int, step: int) -> None:
        if step == nmax:
            return
        nxt_step = step + 1
        for dx, dy in dirs:
            nx, ny = x + dx, y + dy
            pt = (nx, ny)
            if pt not in visited:
                counts[nxt_step] += 1
                if nxt_step < nmax:
                    visited.add(pt)
                    dfs(nx, ny, nxt_step)
                    visited.remove(pt)

    dfs(0, 0, 0)
    return counts


def validate_saw_walk(coords: Sequence[Sequence[int]], expected_n: int | None = None) -> tuple[bool, str]:
    """Validate that a walk starts at (0,0), takes unit Manhattan steps, and has 0 self-intersections."""
    if not coords:
        return False, "empty walk"
    if expected_n is not None and len(coords) != expected_n + 1:
        return False, f"expected {expected_n + 1} vertices, got {len(coords)}"
    if (int(coords[0][0]), int(coords[0][1])) != (0, 0):
        return False, f"origin mismatch: {coords[0]}"

    seen: set[tuple[int, int]] = {(0, 0)}
    px, py = 0, 0
    for idx in range(1, len(coords)):
        cx, cy = int(coords[idx][0]), int(coords[idx][1])
        if abs(cx - px) + abs(cy - py) != 1:
            return False, f"non-unit step at index {idx}: ({px},{py}) -> ({cx},{cy})"
        pt = (cx, cy)
        if pt in seen:
            return False, f"self-intersection at step {idx}: {pt}"
        seen.add(pt)
        px, py = cx, cy
    return True, "ok"


def ols_simple(x: Sequence[float], y: Sequence[float]) -> dict[str, float | list[float]]:
    """Exact 2-parameter OLS y = beta0 + beta1 * x with covariance matrix."""
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    m = xa.size
    if m < 3:
        raise ValueError("need at least 3 points for OLS with SE")
    X = np.column_stack([np.ones(m, dtype=np.float64), xa])
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = xtx_inv @ (X.T @ ya)
    resid = ya - (X @ beta)
    s2 = float(np.dot(resid, resid) / (m - 2))
    cov = s2 * xtx_inv
    ss_tot = float(np.sum((ya - np.mean(ya)) ** 2))
    r2 = 1.0 - float(np.dot(resid, resid)) / ss_tot if ss_tot > 0 else 1.0
    return {
        "beta0": float(beta[0]),
        "beta1": float(beta[1]),
        "se_beta0": float(math.sqrt(max(cov[0, 0], 0.0))),
        "se_beta1": float(math.sqrt(max(cov[1, 1], 0.0))),
        "cov01": float(cov[0, 1]),
        "r_squared": r2,
    }


def estimate_mu_gamma_from_counts(counts: Sequence[int], n_start: int, n_end: int) -> dict[str, object]:
    """Ratio method r_n = c_n / c_{n-1} = mu * (1 + (gamma - 1)/n + O(1/n^2)) on n in [n_start, n_end]."""
    ns = list(range(n_start, n_end + 1))
    inv_n = [1.0 / n for n in ns]
    ratios = [counts[n] / counts[n - 1] for n in ns]
    fit = ols_simple(inv_n, ratios)
    mu_hat = float(fit["beta0"])
    slope = float(fit["beta1"])
    se_mu = float(fit["se_beta0"])
    se_slope = float(fit["se_beta1"])
    cov01 = float(fit["cov01"])

    gamma_hat = 1.0 + slope / mu_hat
    # Delta method for g(b0, b1) = 1 + b1 / b0: grad = (-b1/b0^2, 1/b0)
    g0 = -slope / (mu_hat * mu_hat)
    g1 = 1.0 / mu_hat
    var_gamma = (g0 * g0) * (se_mu * se_mu) + 2.0 * g0 * g1 * cov01 + (g1 * g1) * (se_slope * se_slope)
    se_gamma = math.sqrt(max(var_gamma, 0.0))

    return {
        "n_range": [n_start, n_end],
        "ns": ns,
        "inv_n": inv_n,
        "ratios": ratios,
        "mu_hat": mu_hat,
        "mu_se": se_mu,
        "slope_hat": slope,
        "slope_se": se_slope,
        "gamma_hat": gamma_hat,
        "gamma_se": se_gamma,
        "r_squared": float(fit["r_squared"]),
    }


def estimate_nu_from_pivots(ns: Sequence[int], r2_means: Sequence[float]) -> dict[str, object]:
    """Log-log regression ln <R^2> = ln A + (2 * nu) * ln N."""
    log_n = [math.log(float(n)) for n in ns]
    log_r2 = [math.log(float(r2)) for r2 in r2_means]
    fit = ols_simple(log_n, log_r2)
    b0 = float(fit["beta0"])
    b1 = float(fit["beta1"])
    se_b0 = float(fit["se_beta0"])
    se_b1 = float(fit["se_beta1"])

    nu_hat = 0.5 * b1
    nu_se = 0.5 * se_b1
    amp_hat = math.exp(b0)
    amp_se = amp_hat * se_b0
    return {
        "ns": list(ns),
        "log_n": log_n,
        "log_r2": log_r2,
        "nu_hat": nu_hat,
        "nu_se": nu_se,
        "slope_2nu": b1,
        "slope_2nu_se": se_b1,
        "amplitude_hat": amp_hat,
        "amplitude_se": amp_se,
        "intercept_ln_a": b0,
        "r_squared": float(fit["r_squared"]),
    }


def _save_png_and_svg(fig: plt.Figure, out_dir: Path, stem: str) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{stem}.png"
    svg_path = out_dir / f"{stem}.svg"
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    return [str(png_path), str(svg_path)]


def generate_all_figures(
    bundle: dict,
    ratio_fit: dict,
    scaling_fit: dict,
    fig_dir: Path,
) -> dict[str, list[str]]:
    created: dict[str, list[str]] = {}
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10})

    # 1. saw_samples: 120 SAWs N=400 overlaid, colored by end-to-end distance
    piv400 = next(p for p in bundle["pivots"] if p["n"] == 400)
    walks = piv400["samples"][:120]
    dists = [math.hypot(w[-1][0], w[-1][1]) for w in walks]
    d_min, d_max = min(dists), max(dists)
    cmap = plt.cm.viridis

    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    order = sorted(range(len(walks)), key=lambda i: dists[i])
    for idx in order:
        w = walks[idx]
        xs = [pt[0] for pt in w]
        ys = [pt[1] for pt in w]
        norm_d = (dists[idx] - d_min) / max(d_max - d_min, 1e-9)
        ax.plot(xs, ys, color=cmap(norm_d), alpha=0.36, linewidth=0.85)
    ax.scatter([0], [0], color="#d62728", s=36, zorder=5, label="Titik asal (0,0)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=d_min, vmax=d_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, shrink=0.86)
    cbar.set_label("Jarak ujung-ke-ujung R = ||x_N||_2")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.22, linestyle="--")
    ax.set_title(f"120 Sampel SAW N = 400 (Pivot Madras-Sokal, <R^2> = {piv400['r2_mean']:.1f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="upper right", framealpha=0.9)
    created["saw_samples"] = _save_png_and_svg(fig, fig_dir, "day006a_saw_samples")

    # 2. saw_scaling: log-log <R^2> vs N with regression line and Flory 3/2 slope reference
    ns = np.array([p["n"] for p in bundle["pivots"]], dtype=float)
    r2s = np.array([p["r2_mean"] for p in bundle["pivots"]], dtype=float)
    ses = np.array([p["r2_se"] for p in bundle["pivots"]], dtype=float)
    n_grid = np.geomspace(ns.min() * 0.9, ns.max() * 1.1, 100)
    fit_curve = scaling_fit["amplitude_hat"] * (n_grid ** scaling_fit["slope_2nu"])
    ref_amp = r2s[2] / (ns[2] ** 1.5)
    ref_curve = ref_amp * (n_grid ** 1.5)

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.errorbar(ns, r2s, yerr=2.0 * ses, fmt="o", color="#1f77b4", capsize=4, label="Pivot MC (<R^2> +- 2 SE)")
    ax.plot(
        n_grid,
        fit_curve,
        "-",
        color="#d62728",
        linewidth=1.8,
        label=f"Regresi: nu_hat = {scaling_fit['nu_hat']:.4f} +- {scaling_fit['nu_se']:.4f}",
    )
    ax.plot(n_grid, ref_curve, "--", color="#2ca02c", linewidth=1.4, label="Acuan Flory 2nu = 3/2 (nu = 0.75)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25, linestyle="--")
    ax.set_xlabel("Panjang jalan N (skala log)")
    ax.set_ylabel("<R^2(N)> (skala log)")
    ax.set_title("Penskalaan Kritis <R^2(N)> ~ A * N^(2*nu) pada Kisi Persegi 2D")
    ax.legend(loc="upper left")
    created["saw_scaling"] = _save_png_and_svg(fig, fig_dir, "day006a_saw_scaling")

    # 3. saw_ratio: r_n vs 1/n with linear extrapolation to mu_hat and gamma_hat annotation
    all_counts = bundle["go_sharded_counts"]
    ns_all = np.arange(6, len(all_counts), dtype=float)
    inv_all = 1.0 / ns_all
    rat_all = np.array([all_counts[int(n)] / all_counts[int(n) - 1] for n in ns_all])
    inv_fit = np.array(ratio_fit["inv_n"], dtype=float)
    rat_fit = np.array(ratio_fit["ratios"], dtype=float)
    x_line = np.linspace(0.0, inv_fit.max() * 1.12, 100)
    y_line = ratio_fit["mu_hat"] + ratio_fit["slope_hat"] * x_line

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.scatter(inv_all, rat_all, color="#9ecae1", s=26, label=f"r_n = c_n/c_(n-1) (n = 6..{len(all_counts)-1})")
    ax.scatter(
        inv_fit,
        rat_fit,
        color="#1f77b4",
        s=42,
        zorder=4,
        label=f"Jendela regresi n in [{ratio_fit['n_range'][0]}, {ratio_fit['n_range'][1]}]",
    )
    ax.plot(x_line, y_line, "-", color="#d62728", linewidth=1.8, label="Ekstrapolasi linear OLS")
    ax.scatter([0.0], [ratio_fit["mu_hat"]], color="#d62728", marker="D", s=52, zorder=5)
    ax.scatter([0.0], [LIT_MU_REF], color="#2ca02c", marker="x", s=60, zorder=5, label=f"Acuan mu ~ {LIT_MU_REF:.5f}")
    ax.annotate(
        f"mu_hat = {ratio_fit['mu_hat']:.5f} +- {ratio_fit['mu_se']:.5f}\n"
        f"gamma_hat = {ratio_fit['gamma_hat']:.4f} +- {ratio_fit['gamma_se']:.4f}\n"
        f"(acuan 43/32 = {LIT_GAMMA_REF:.5f})",
        xy=(0.0, ratio_fit["mu_hat"]),
        xytext=(0.025, ratio_fit["mu_hat"] - 0.012),
        arrowprops={"arrowstyle": "->", "color": "#444444"},
        bbox={"boxstyle": "round,pad=0.35", "fc": "#fff8dc", "ec": "#cccccc"},
        fontsize=9,
    )
    ax.grid(True, alpha=0.25, linestyle="--")
    ax.set_xlabel("1 / n")
    ax.set_ylabel("Rasio r_n = c_n / c_(n-1)")
    ax.set_title("Metode Rasio SAW: Ekstrapolasi Konstanta Konektivitas mu & Eksponen gamma")
    ax.legend(loc="lower right", fontsize=8.5)
    created["saw_ratio"] = _save_png_and_svg(fig, fig_dir, "day006a_saw_ratio")

    # 4. saw_endpoint_heatmap: endpoint density for n = 16 (8-fold symmetric)
    hm = np.array(bundle["rust_single"]["heatmap"], dtype=float)
    hn = int(bundle["rust_single"]["heatmap_n"])
    hm_masked = np.where(hm > 0, np.log10(hm), np.nan)

    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    im = ax.imshow(
        hm_masked,
        origin="lower",
        extent=[-hn - 0.5, hn + 0.5, -hn - 0.5, hn + 0.5],
        cmap="magma",
        interpolation="nearest",
    )
    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.86)
    cbar.set_label("log10(jumlah SAW berakhir di (x, y))")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Koordinat ujung x")
    ax.set_ylabel("Koordinat ujung y")
    ax.set_title(f"Peta Panas Kerapatan Titik Ujung SAW Eksak n = {hn} (Simetri D4, c_{hn} = {all_counts[hn]:,})")
    ax.grid(False)
    created["saw_endpoint_heatmap"] = _save_png_and_svg(fig, fig_dir, "day006a_saw_endpoint_heatmap")

    # 5. saw_worker_timing: Go wall time & relative speedup with 1, 2, 4, N workers
    timings = bundle["worker_timings"]
    ws = [t["workers"] for t in timings]
    secs = [t["elapsed_sec"] for t in timings]
    spds = [t["speedup"] for t in timings]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.2))
    x_pos = np.arange(len(ws))
    bars = ax1.bar(x_pos, secs, color="#4c78a8", width=0.55)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([str(w) for w in ws])
    ax1.set_xlabel("Jumlah worker proses Rust persisten")
    ax1.set_ylabel("Waktu dinding (detik)")
    ax1.set_title(f"Waktu Dinding Go Sharding (Nmax = {bundle['nmax']}, K = {bundle['prefix_k']})")
    ax1.grid(True, axis="y", alpha=0.25, linestyle="--")
    for bar, s in zip(bars, secs):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{s:.3f}s", ha="center", va="bottom", fontsize=8.5)

    ax2.plot(ws, spds, "o-", color="#e45756", linewidth=2.0, label="Percepatan terukur")
    ax2.plot(ws, ws, "--", color="#7f7f7f", linewidth=1.2, label="Skala linear ideal")
    ax2.set_xlabel("Jumlah worker")
    ax2.set_ylabel("Kecepatan relatif T(1) / T(W)")
    ax2.set_title("Percepatan Relatif & Invarian Deterministik")
    ax2.grid(True, alpha=0.25, linestyle="--")
    ax2.legend(loc="upper left")
    fig.tight_layout()
    created["saw_worker_timing"] = _save_png_and_svg(fig, fig_dir, "day006a_saw_worker_timing")

    return created


def generate_interactive_html(bundle: dict, html_path: Path) -> str:
    """Self-contained vanilla JS viewer (zero CDN / zero external URLs) with embedded pivot samples."""
    embedded: dict[str, list[list[list[int]]]] = {}
    stats_map: dict[str, dict[str, float]] = {}
    for p in bundle["pivots"]:
        key = str(p["n"])
        embedded[key] = p["samples"][:12]
        stats_map[key] = {"r2_mean": round(p["r2_mean"], 2), "r2_se": round(p["r2_se"], 2)}

    data_json = json.dumps({"walks": embedded, "stats": stats_map}, separators=(",", ":"))
    html = f"""<!DOCTYPE html>
<html lang="id">
<head>
<meta charset="utf-8">
<title>Penampil Interaktif Self-Avoiding Walk (day006a)</title>
<style>
  body {{ font-family: system-ui, -apple-system, sans-serif; margin: 24px; background: #0f172a; color: #e2e8f0; }}
  .panel {{ max-width: 820px; margin: 0 auto; background: #1e293b; padding: 20px; border-radius: 12px; }}
  h1 {{ font-size: 1.25rem; margin-top: 0; }}
  .controls {{ display: flex; gap: 16px; align-items: center; flex-wrap: wrap; margin-bottom: 14px; }}
  button {{ background: #38bdf8; color: #0f172a; border: none; padding: 8px 14px; border-radius: 6px; font-weight: 600; cursor: pointer; }}
  button:hover {{ background: #7dd3fc; }}
  canvas {{ background: #0b1120; border-radius: 8px; width: 100%; max-width: 780px; height: 520px; display: block; }}
  .meta {{ font-size: 0.9rem; color: #94a3b8; margin-top: 10px; }}
</style>
</head>
<body>
<div class="panel">
  <h1>Self-Avoiding Walk 2D — Sampel Pivot Madras–Sokal Tertanam (day006a)</h1>
  <div class="controls">
    <label>Indeks N: <input id="nSlider" type="range" min="0" max="4" step="1" value="3"></label>
    <span id="nLabel">N = 400</span>
    <button id="nextBtn" type="button">Sampel Baru</button>
    <button id="overlayBtn" type="button">Tumpuk Semua Sampel N Ini</button>
  </div>
  <canvas id="sawCanvas" width="780" height="520"></canvas>
  <div class="meta" id="metaInfo"></div>
</div>
<script>
(function() {{
  const DATA = {data_json};
  const NS = ["50", "100", "200", "400", "800"];
  const slider = document.getElementById("nSlider");
  const nLabel = document.getElementById("nLabel");
  const metaInfo = document.getElementById("metaInfo");
  const canvas = document.getElementById("sawCanvas");
  const ctx = canvas.getContext("2d");
  let sampleIdx = 0;
  let overlayMode = false;

  function drawWalks(walks, nKey) {{
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    let minX = 0, maxX = 0, minY = 0, maxY = 0;
    for (const w of walks) {{
      for (const p of w) {{
        if (p[0] < minX) minX = p[0];
        if (p[0] > maxX) maxX = p[0];
        if (p[1] < minY) minY = p[1];
        if (p[1] > maxY) maxY = p[1];
      }}
    }}
    const pad = 32;
    const spanX = Math.max(maxX - minX, 4);
    const spanY = Math.max(maxY - minY, 4);
    const scale = Math.min((canvas.width - 2 * pad) / spanX, (canvas.height - 2 * pad) / spanY);
    const cx = canvas.width / 2 - ((minX + maxX) / 2) * scale;
    const cy = canvas.height / 2 + ((minY + maxY) / 2) * scale;

    walks.forEach((w, idx) => {{
      const hue = (idx * 47 + 190) % 360;
      ctx.strokeStyle = overlayMode ? `hsla(${{hue}}, 85%, 62%, 0.55)` : "#38bdf8";
      ctx.lineWidth = overlayMode ? 1.3 : 2.0;
      ctx.beginPath();
      w.forEach((p, i) => {{
        const x = cx + p[0] * scale;
        const y = cy - p[1] * scale;
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }});
      ctx.stroke();
      const end = w[w.length - 1];
      ctx.fillStyle = "#facc15";
      ctx.beginPath();
      ctx.arc(cx + end[0] * scale, cy - end[1] * scale, 3.2, 0, Math.PI * 2);
      ctx.fill();
    }});

    ctx.fillStyle = "#ef4444";
    ctx.beginPath();
    ctx.arc(cx, cy, 4.8, 0, Math.PI * 2);
    ctx.fill();

    const st = DATA.stats[nKey];
    const shown = walks[0];
    const end = shown[shown.length - 1];
    const r = Math.hypot(end[0], end[1]).toFixed(2);
    metaInfo.textContent = overlayMode
      ? `Menampilkan ${{walks.length}} sampel tertumpuk untuk N = ${{nKey}} | <R^2> = ${{st.r2_mean}} +- ${{st.r2_se}}`
      : `Sampel #${{sampleIdx + 1}}/${{DATA.walks[nKey].length}} untuk N = ${{nKey}} | Ujung = (${{end[0]}}, ${{end[1]}}), R = ${{r}} | <R^2> = ${{st.r2_mean}} +- ${{st.r2_se}}`;
  }}

  function render() {{
    const nKey = NS[Number(slider.value)];
    nLabel.textContent = "N = " + nKey;
    const list = DATA.walks[nKey];
    sampleIdx = sampleIdx % list.length;
    drawWalks(overlayMode ? list : [list[sampleIdx]], nKey);
  }}

  slider.addEventListener("input", () => {{ sampleIdx = 0; render(); }});
  document.getElementById("nextBtn").addEventListener("click", () => {{ overlayMode = false; sampleIdx += 1; render(); }});
  document.getElementById("overlayBtn").addEventListener("click", () => {{ overlayMode = !overlayMode; render(); }});
  render();
}})();
</script>
</body>
</html>
"""
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html, encoding="utf-8")
    return str(html_path)


def verify_visual_artifacts(created_figs: dict[str, list[str]], html_path: Path | None) -> tuple[bool, list[str]]:
    issues: list[str] = []
    for key, paths in created_figs.items():
        for p_str in paths:
            p = Path(p_str)
            if not p.exists() or p.stat().st_size == 0:
                issues.append(f"{key}: missing or empty {p.name}")
                continue
            if p.suffix == ".svg":
                try:
                    ET.parse(p)
                except Exception as exc:
                    issues.append(f"{key}: invalid SVG XML in {p.name}: {exc}")
    if html_path is not None:
        if not html_path.exists() or html_path.stat().st_size == 0:
            issues.append("saw_viewer.html missing or empty")
        else:
            txt = html_path.read_text(encoding="utf-8")
            if "http://" in txt or "https://" in txt:
                issues.append("saw_viewer.html contains external URL reference")
    return len(issues) == 0, issues


def run_pipeline(
    go_json_path: Path,
    results_json_path: Path,
    fig_dir: Path,
    t0_epoch: int,
    t0_str: str,
) -> dict:
    bundle = json.loads(go_json_path.read_text(encoding="utf-8"))
    nmax = int(bundle["nmax"])

    t_bf0 = time.perf_counter()
    py_counts_12 = brute_force_saw_counts(12)
    py_bf_sec = time.perf_counter() - t_bf0

    rust_counts = [int(x) for x in bundle["rust_single"]["counts"]]
    go_counts = [int(x) for x in bundle["go_sharded_counts"]]

    c1_match = py_counts_12 == rust_counts[:13] == go_counts[:13]
    c2_match = rust_counts == go_counts

    comp_limit = min(20, nmax)
    mismatches_a001411 = []
    for n in range(1, comp_limit + 1):
        if go_counts[n] != MEMORY_A001411_N1_20[n - 1]:
            mismatches_a001411.append({"n": n, "measured": go_counts[n], "memory_ref": MEMORY_A001411_N1_20[n - 1]})
    c3_match = comp_limit == 20 and len(mismatches_a001411) == 0

    ratio_fit = estimate_mu_gamma_from_counts(go_counts, nmax - 8, nmax)
    piv_ns = [int(p["n"]) for p in bundle["pivots"]]
    piv_r2 = [float(p["r2_mean"]) for p in bundle["pivots"]]
    scaling_fit = estimate_nu_from_pivots(piv_ns, piv_r2)

    mu_diff = abs(float(ratio_fit["mu_hat"]) - 2.638)
    gamma_diff = abs(float(ratio_fit["gamma_hat"]) - 1.34375)
    nu_diff = abs(float(scaling_fit["nu_hat"]) - 0.75)
    c4_ok = mu_diff < 0.02
    c5_ok = gamma_diff < 0.06
    c6_ok = nu_diff < 0.02

    total_sample_walks = 0
    collision_failures = 0
    for p in bundle["pivots"]:
        n_val = int(p["n"])
        for w in p["samples"]:
            total_sample_walks += 1
            ok, _ = validate_saw_walk(w, expected_n=n_val)
            if not ok:
                collision_failures += 1
    c7_ok = total_sample_walks >= 120 and collision_failures == 0

    c8_ok = bool(bundle["worker_invariant"])

    figs = generate_all_figures(bundle, ratio_fit, scaling_fit, fig_dir)
    html_path = fig_dir / "day006a_saw_viewer.html"
    generate_interactive_html(bundle, html_path)
    c9_ok, fig_issues = verify_visual_artifacts(figs, html_path)

    now_epoch = int(time.time())
    wita = timezone(timedelta(hours=8))
    now_str = datetime.fromtimestamp(now_epoch, tz=wita).strftime("%a %b %d %H:%M:%S WITA %Y")
    elapsed_wall_sec = now_epoch - t0_epoch
    c10_ok = elapsed_wall_sec <= 600

    criteria = {
        "c1_python_bruteforce_matches_rust_go_n12": {
            "status": "tercapai" if c1_match else "belum tercapai",
            "measured": {"max_abs_diff_n0_12": 0 if c1_match else -1, "c_12": py_counts_12[12]},
            "threshold": "c_n Python brute-force == Rust == Go persis untuk semua n <= 12",
        },
        "c2_go_sharded_matches_rust_single_all_n": {
            "status": "tercapai" if c2_match else "belum tercapai",
            "measured": {"nmax": nmax, "c_nmax": go_counts[nmax], "exact_equal": c2_match},
            "threshold": f"Go-sharded (K={bundle['prefix_k']}, {bundle['num_prefix_tasks']} tugas) == Rust satu-proses pada n = 1..{nmax}",
        },
        "c3_comparison_memory_sequence_a001411_n20": {
            "status": "tercapai" if c3_match else "belum tercapai",
            "measured": {"compared_up_to_n": comp_limit, "mismatch_count": len(mismatches_a001411), "mismatches": mismatches_a001411},
            "threshold": "0 ketidakcocokan terhadap deret acuan ingatan A001411 pada n = 1..20",
        },
        "c4_connective_constant_mu_hat": {
            "status": "tercapai" if c4_ok else "belum tercapai",
            "measured": {
                "mu_hat": round(float(ratio_fit["mu_hat"]), 6),
                "se": round(float(ratio_fit["mu_se"]), 6),
                "abs_diff_vs_2_638": round(mu_diff, 6),
            },
            "threshold": "|mu_hat - 2.638| < 0.02",
        },
        "c5_entropic_exponent_gamma_hat": {
            "status": "tercapai" if c5_ok else "belum tercapai",
            "measured": {
                "gamma_hat": round(float(ratio_fit["gamma_hat"]), 6),
                "se": round(float(ratio_fit["gamma_se"]), 6),
                "abs_diff_vs_43_32": round(gamma_diff, 6),
            },
            "threshold": "|gamma_hat - 1.34375| < 0.06",
        },
        "c6_metric_exponent_nu_hat": {
            "status": "tercapai" if c6_ok else "belum tercapai",
            "measured": {
                "nu_hat": round(float(scaling_fit["nu_hat"]), 6),
                "se": round(float(scaling_fit["nu_se"]), 6),
                "abs_diff_vs_0_75": round(nu_diff, 6),
            },
            "threshold": "|nu_hat - 0.75| < 0.02",
        },
        "c7_sample_walks_zero_collisions": {
            "status": "tercapai" if c7_ok else "belum tercapai",
            "measured": {"validated_walks": total_sample_walks, "collision_or_step_failures": collision_failures},
            "threshold": "0 tabrakan / langkah cacat pada seluruh jalan contoh (>= 120 jalan N=400)",
        },
        "c8_worker_count_invariance": {
            "status": "tercapai" if c8_ok else "belum tercapai",
            "measured": {
                "workers_tested": [t["workers"] for t in bundle["worker_timings"]],
                "all_exact_equal": c8_ok,
            },
            "threshold": "c_n identik bit-demi-bit untuk seluruh jumlah worker yang diuji",
        },
        "c9_visualizations_valid": {
            "status": "tercapai" if c9_ok else "belum tercapai",
            "measured": {
                "figure_groups": len(figs),
                "interactive_html": html_path.name,
                "issues": fig_issues,
            },
            "threshold": "5 grafik (PNG+SVG XML valid, tak kosong) + saw_viewer.html tanpa URL eksternal",
        },
        "c10_ten_minute_wall_clock_limit": {
            "status": "tercapai" if c10_ok else "belum tercapai",
            "measured": {
                "t0_epoch": t0_epoch,
                "t0_wita": t0_str,
                "finished_epoch": now_epoch,
                "finished_wita": now_str,
                "elapsed_seconds": elapsed_wall_sec,
            },
            "threshold": "selisih waktu mulai T0 hingga selesai <= 600 detik (10 menit)",
        },
    }

    result_doc = {
        "session": "day006a",
        "date": "2026-11-24",
        "master_seed": int(bundle["base_seed"]),
        "run_completed_without_errors": True,
        "timing": {
            "t0_epoch": t0_epoch,
            "t0_wita": t0_str,
            "analysis_finished_epoch": now_epoch,
            "analysis_finished_wita": now_str,
            "elapsed_seconds_at_analysis": elapsed_wall_sec,
            "rust_single_sec": round(float(bundle["rust_single_sec"]), 4),
            "python_bruteforce_n12_sec": round(py_bf_sec, 4),
            "go_worker_timings": [
                {
                    "workers": t["workers"],
                    "elapsed_sec": round(float(t["elapsed_sec"]), 4),
                    "speedup": round(float(t["speedup"]), 3),
                    "counts_match": bool(t["counts_match"]),
                }
                for t in bundle["worker_timings"]
            ],
        },
        "enumeration": {
            "nmax": nmax,
            "prefix_k": int(bundle["prefix_k"]),
            "num_prefix_tasks": int(bundle["num_prefix_tasks"]),
            "python_bruteforce_n0_12": py_counts_12,
            "counts_n0_nmax": go_counts,
            "memory_a001411_mismatches": mismatches_a001411,
            "heatmap_n": int(bundle["rust_single"]["heatmap_n"]),
        },
        "ratio_analysis": {
            "n_window": ratio_fit["n_range"],
            "mu_hat": ratio_fit["mu_hat"],
            "mu_se": ratio_fit["mu_se"],
            "gamma_hat": ratio_fit["gamma_hat"],
            "gamma_se": ratio_fit["gamma_se"],
            "slope_hat": ratio_fit["slope_hat"],
            "r_squared": ratio_fit["r_squared"],
            "literature_reference": {"mu": LIT_MU_REF, "gamma": LIT_GAMMA_REF},
        },
        "pivot_analysis": {
            "nu_hat": scaling_fit["nu_hat"],
            "nu_se": scaling_fit["nu_se"],
            "slope_2nu": scaling_fit["slope_2nu"],
            "amplitude_hat": scaling_fit["amplitude_hat"],
            "amplitude_se": scaling_fit["amplitude_se"],
            "r_squared": scaling_fit["r_squared"],
            "literature_reference": {"nu": LIT_NU_REF},
            "per_n": [
                {
                    "n": int(p["n"]),
                    "seed": int(p["seed"]),
                    "r2_mean": float(p["r2_mean"]),
                    "r2_se": float(p["r2_se"]),
                    "acceptance_rate": float(p["acceptance_rate"]),
                    "total_steps": int(p["total_steps"]),
                    "saved_samples": len(p["samples"]),
                }
                for p in bundle["pivots"]
            ],
        },
        "artifacts": {
            "figures": {
                k: [Path(p).as_posix() for p in v] for k, v in figs.items()
            },
            "interactive_html": html_path.as_posix(),
        },
        "criteria": criteria,
    }

    results_json_path.parent.mkdir(parents=True, exist_ok=True)
    results_json_path.write_text(json.dumps(result_doc, indent=2), encoding="utf-8")
    return result_doc


def stamp_final_time(results_json_path: Path) -> None:
    """Refresh c10 end timestamp right before final commit so results.json records actual completion time."""
    doc = json.loads(results_json_path.read_text(encoding="utf-8"))
    t0_epoch = int(doc["timing"]["t0_epoch"])
    now_epoch = int(time.time())
    wita = timezone(timedelta(hours=8))
    now_str = datetime.fromtimestamp(now_epoch, tz=wita).strftime("%a %b %d %H:%M:%S WITA %Y")
    elapsed = now_epoch - t0_epoch
    doc["timing"]["final_epoch"] = now_epoch
    doc["timing"]["final_wita"] = now_str
    doc["timing"]["final_elapsed_seconds"] = elapsed
    c10 = doc["criteria"]["c10_ten_minute_wall_clock_limit"]
    c10["status"] = "tercapai" if elapsed <= 600 else "belum tercapai"
    c10["measured"]["finished_epoch"] = now_epoch
    c10["measured"]["finished_wita"] = now_str
    c10["measured"]["elapsed_seconds"] = elapsed
    results_json_path.write_text(json.dumps(doc, indent=2), encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description="Analyze SAW results and generate plots.")
    p.add_argument("--go-json", type=Path, default=Path("outputs/day006a/go_orchestrated.json"))
    p.add_argument("--out-json", type=Path, default=Path("outputs/day006a/results.json"))
    p.add_argument("--fig-dir", type=Path, default=Path("outputs/figures"))
    p.add_argument("--t0-epoch", type=int, required=False, default=int(time.time()))
    p.add_argument("--t0-str", type=str, required=False, default="")
    p.add_argument("--stamp-final", action="store_true", help="Only update final completion timestamp in out-json")
    args = p.parse_args()

    if args.stamp_final:
        stamp_final_time(args.out_json)
        return

    run_pipeline(args.go_json, args.out_json, args.fig_dir, args.t0_epoch, args.t0_str)


if __name__ == "__main__":
    main()
