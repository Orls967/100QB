"""Visualisasi untuk day003a: permukaan SSVI, sigma_loc, g(y, T), konvergensi PDE, MC, barrier, vega-bucket, kalibrasi."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from qb100.plotting import PALETTE, save_figure, set_style


def make_all_figures(results: dict, folder: str | Path | None = None) -> list[Path]:
    """Hasilkan seluruh grafik day003a ke ``outputs/figures/``."""
    set_style()
    paths: list[Path] = []

    # 1. Permukaan implied vol & heatmap sigma_loc(y, T)
    surf = results["_plot_data"]["surfaces"]
    y_g = surf["y_grid"]
    T_g = surf["T_grid"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    im0 = axes[0].contourf(T_g, y_g, surf["iv_grid"] * 100.0, levels=20, cmap="viridis")
    fig.colorbar(im0, ax=axes[0], label="Implied Vol (%)")
    axes[0].set_title("Permukaan Implied Volatility SSVI")
    axes[0].set_xlabel("Maturitas T (tahun)")
    axes[0].set_ylabel("Log-moneyness y = ln(K/F)")

    im1 = axes[1].contourf(T_g, y_g, surf["sigma_loc_grid"] * 100.0, levels=20, cmap="magma")
    fig.colorbar(im1, ax=axes[1], label="Local Vol (%)")
    axes[1].set_title(r"Heatmap Volatilitas Lokal Dupire $\sigma_{\mathrm{loc}}(y, T)$")
    axes[1].set_xlabel("Maturitas T (tahun)")
    axes[1].set_ylabel("Log-moneyness y = ln(K/F)")
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_ssvi-implied-vol-and-local-vol", folder=folder))
    plt.close(fig)

    # 2. g(y, T) valid vs kontrafaktual
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    im_g = axes[0].contourf(T_g, y_g, surf["g_grid"], levels=20, cmap="coolwarm")
    fig.colorbar(im_g, ax=axes[0], label="g(y, T)")
    axes[0].set_title("Densitas Kupu-kupu g(y, T) — Parameter Valid")
    axes[0].set_xlabel("Maturitas T (tahun)")
    axes[0].set_ylabel("Log-moneyness y")

    for T_idx in [0, len(T_g) // 2, -1]:
        axes[1].plot(y_g, surf["g_grid"][:, T_idx], label=f"Valid T={T_g[T_idx]:.2f}")
    axes[1].plot(y_g, surf["g_bad_slice"], "--", color=PALETTE[1], label="Kontrafaktual (eta=4, rho=-0.9)")
    axes[1].axhline(0.0, color="black", lw=1.0, ls=":")
    axes[1].set_title("Irisan g(y, T) vs Uji Kontrafaktual")
    axes[1].set_xlabel("Log-moneyness y")
    axes[1].set_ylabel("g(y, T)")
    axes[1].legend()
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_butterfly-density-g", folder=folder))
    plt.close(fig)

    # 3. Kurva konvergensi PDE log-log
    conv = results["part3"]["flat_convergence"]
    ns = np.asarray(conv["resolutions"], dtype=float)
    eb = np.asarray(conv["backward_errors"], dtype=float)
    ef = np.asarray(conv["forward_errors"], dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(ns, eb, "o-", label=f"Backward PDE (slope={conv['backward_slope']:.3f})")
    ax.loglog(ns, ef, "s-", label=f"Forward Dupire PDE (slope={conv['forward_slope']:.3f})")
    ref_line = eb[0] * (ns / ns[0]) ** (-2.0)
    ax.loglog(ns, ref_line, "k--", alpha=0.6, label="Referensi orde 2 (-2.0)")
    ax.set_title("Konvergensi PDE Kasus Vol Datar (Log-Log)")
    ax.set_xlabel("Resolusi Grid N")
    ax.set_ylabel("Galat Harga Absolut")
    ax.legend()
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_pde-convergence", folder=folder))
    plt.close(fig)

    # 4. Heatmap galat vol implisit 5x4 (bp)
    ssvi_pde = results["part3"]["ssvi_pde"]
    m_vals = ssvi_pde["m_vals"]
    T_vals = ssvi_pde["T_vals"]
    err_b_bp = np.asarray(ssvi_pde["iv_err_backward_bp"], dtype=float)
    err_f_bp = np.asarray(ssvi_pde["iv_err_forward_bp"], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, mat, title in [
        (axes[0], err_b_bp, "Galat IV Backward PDE (bp)"),
        (axes[1], err_f_bp, "Galat IV Forward Dupire PDE (bp)"),
    ]:
        im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")
        ax.set_xticks(range(len(T_vals)), [f"{t:.2f}" for t in T_vals])
        ax.set_yticks(range(len(m_vals)), [f"{m:.1f}" for m in m_vals])
        ax.set_xlabel("Maturitas T")
        ax.set_ylabel("Moneyness K/F")
        ax.set_title(title)
        for i in range(len(m_vals)):
            for j in range(len(T_vals)):
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, label="bp")
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_iv-errors-5x4", folder=folder))
    plt.close(fig)

    # 5. MC vs PDE (z-score) & Weak Order
    mc_res = results["part3"]["mc_euler"]
    weak = results["part3"]["mc_weak_order"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].bar(range(len(m_vals)), mc_res["z_scores"], color=PALETTE[0])
    axes[0].axhline(3.0, color=PALETTE[1], ls="--", label="Batas 3 SE")
    axes[0].set_xticks(range(len(m_vals)), [f"{m:.1f}" for m in m_vals])
    axes[0].set_xlabel("Moneyness K/F (T=1)")
    axes[0].set_ylabel("|MC - PDE| / SE")
    axes[0].set_title("Monte Carlo Euler vs PDE (z-score)")
    axes[0].legend()

    steps = np.asarray(weak["steps"], dtype=float)
    werr = np.asarray(weak["abs_errors"], dtype=float)
    axes[1].loglog(steps, werr, "o-", label=f"Euler CRN (slope={weak['slope']:.3f})")
    axes[1].loglog(steps, werr[0] * (steps / steps[0]) ** (-1.0), "k--", alpha=0.6, label="Orde 1 (-1.0)")
    axes[1].set_xlabel("Langkah per Tahun")
    axes[1].set_ylabel("Galat Diskretisasi")
    axes[1].set_title("Weak Convergence Order Monte Carlo Euler")
    axes[1].legend()
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_mc-vs-pde-and-weak-order", folder=folder))
    plt.close(fig)

    # 6. Barrier & Vega-Bucket
    bar = results["part4"]["barrier"]
    vb = results["part4"]["vega_bucket"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    labels = ["Flat RR", "Flat PDE", "SSVI PDE", "SSVI MC Bridge"]
    vals = [
        bar["flat_rr_down_out"],
        bar["flat_pde_down_out"],
        bar["ssvi_pde_down_out"],
        bar["ssvi_mc_bridge_price"],
    ]
    errs = [0.0, 0.0, 0.0, 2.0 * bar["ssvi_mc_bridge_se"]]
    axes[0].bar(labels, vals, yerr=errs, capsize=5, color=PALETTE[:4])
    axes[0].set_ylabel("Harga Down-and-Out Call (H=85, K=100, T=1)")
    axes[0].set_title("Opsi Barrier: Reiner-Rubinstein, PDE & MC Bridge")

    b_labels = [f"T={t:.2f}" for t in vb["buckets"]] + ["Sum", "Parallel"]
    b_vals = vb["bucket_sensitivities"] + [vb["sum_bucket_sensitivities"], vb["parallel_sens"]]
    axes[1].bar(b_labels, b_vals, color=[PALETTE[0]] * 4 + [PALETTE[2], PALETTE[3]])
    axes[1].set_ylabel(r"$\partial C_{\mathrm{ATM}} / \partial \theta$")
    axes[1].set_title(f"Vega-Bucket vs Pergeseran Paralel (rel diff={vb['rel_diff']:.2e})")
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_barrier-and-vega-bucket", folder=folder))
    plt.close(fig)

    # 7. Distribusi estimasi kalibrasi R=50 & matriks korelasi
    cal = results["part4"]["calibration"]
    est = np.asarray(results["_plot_data"]["calib_estimates"], dtype=float)
    true_p = cal["true_params"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for idx, name in enumerate([r"$\rho$", r"$\eta$", r"$\gamma$"]):
        axes[idx].hist(est[:, idx], bins=12, color=PALETTE[idx], alpha=0.8, edgecolor="white")
        axes[idx].axvline(true_p[idx], color="black", ls="--", lw=1.5, label="Benar")
        axes[idx].axvline(np.mean(est[:, idx]), color=PALETTE[1], ls="-", lw=1.5, label="Mean")
        axes[idx].set_title(f"Distribusi {name} (R=50)")
        axes[idx].legend(fontsize=8)

    corr = np.asarray(cal["correlation_matrix"], dtype=float)
    im_c = axes[3].imshow(corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0)
    axes[3].set_title("Korelasi Parameter (Single Fit)")
    fig.colorbar(im_c, ax=axes[3], fraction=0.046)
    fig.tight_layout()
    paths.append(save_figure(fig, "day003a_calibration-mc-dist", folder=folder))
    plt.close(fig)

    return paths
