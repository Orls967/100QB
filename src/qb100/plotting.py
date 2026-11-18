"""Style grafik seragam dan helper untuk menyimpan figure ke ``outputs/figures``."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

from . import PROJECT_ROOT

FIGURES_DIR = PROJECT_ROOT / "outputs" / "figures"

#: Palet warna seragam untuk semua notebook.
PALETTE = ["#1f6feb", "#d1495b", "#2a9d8f", "#e9c46a", "#6a4c93", "#f4a261"]

_STYLE = {
    "figure.figsize": (10, 5),
    "figure.dpi": 100,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "legend.frameon": False,
    "font.size": 10,
}


def set_style() -> None:
    """Terapkan style grafik seragam 100QB ke matplotlib."""
    plt.rcParams.update(_STYLE)
    plt.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=PALETTE)


def save_figure(
    fig: plt.Figure,
    name: str,
    folder: str | Path | None = None,
    fmt: str = "png",
) -> Path:
    """Simpan figure ke ``outputs/figures`` (atau ``folder``) dan kembalikan path-nya.

    Parameters
    ----------
    fig : matplotlib Figure yang akan disimpan.
    name : nama file tanpa ekstensi, mis. ``"day001_returns"``.
    folder : folder tujuan; default ``outputs/figures`` di root repo.
    fmt : format file (``png``, ``svg``, ``pdf``, ...).
    """
    target_dir = Path(folder) if folder is not None else FIGURES_DIR
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / f"{name}.{fmt}"
    fig.savefig(path, format=fmt)
    return path
