#!/usr/bin/env python
"""Eksekusi semua notebook di ``days/`` untuk memastikan tidak ada error.

Contoh:
    python scripts/run_all.py                 # jalankan semua notebook
    python scripts/run_all.py --day 7         # hanya hari 007
    python scripts/run_all.py --template      # ikut jalankan template
    python scripts/run_all.py --inplace       # simpan output hasil eksekusi ke notebook

Exit code 0 bila semua lulus (termasuk saat folder days/ masih kosong), 1 bila ada yang gagal.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import nbformat
from nbconvert.preprocessors import CellExecutionError, ExecutePreprocessor

ROOT = Path(__file__).resolve().parents[1]
DAYS_DIR = ROOT / "days"
TEMPLATE = ROOT / "templates" / "day_template.ipynb"


def find_notebooks(day: int | None = None, include_template: bool = False) -> list[Path]:
    """Daftar notebook yang akan dijalankan (urut nama file)."""
    pattern = f"day{day:03d}_*.ipynb" if day else "day*.ipynb"
    paths = sorted(DAYS_DIR.glob(pattern)) if DAYS_DIR.exists() else []
    if include_template:
        paths.insert(0, TEMPLATE)
    return paths


def run_notebook(path: Path, timeout: int, kernel: str, inplace: bool) -> tuple[bool, str]:
    """Jalankan satu notebook dari atas ke bawah. Kembalikan (sukses, pesan)."""
    nb = nbformat.read(path, as_version=4)
    ep = ExecutePreprocessor(timeout=timeout, kernel_name=kernel)
    os.environ.setdefault("MPLBACKEND", "Agg")  # tidak membuka jendela grafik
    try:
        ep.preprocess(nb, {"metadata": {"path": str(path.parent)}})
    except CellExecutionError as exc:
        lines = str(exc).strip().splitlines()
        return False, lines[-1] if lines else "CellExecutionError"
    except Exception as exc:  # timeout, kernel mati, dll.
        return False, f"{type(exc).__name__}: {exc}"
    if inplace:
        nbformat.write(nb, path)
    return True, "ok"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Eksekusi semua notebook harian.")
    p.add_argument("--day", type=int, help="hanya jalankan notebook hari tertentu")
    p.add_argument("--template", action="store_true", help="ikut jalankan template")
    p.add_argument("--timeout", type=int, default=600, help="timeout per sel (detik)")
    p.add_argument("--kernel", default="python3", help="nama kernel Jupyter")
    p.add_argument("--inplace", action="store_true", help="simpan output ke notebook")
    args = p.parse_args(argv)

    notebooks = find_notebooks(args.day, args.template)
    if not notebooks:
        print("Tidak ada notebook di days/ — tidak ada yang perlu dijalankan.")
        return 0

    failures = []
    for path in notebooks:
        start = time.perf_counter()
        ok, msg = run_notebook(path, args.timeout, args.kernel, args.inplace)
        elapsed = time.perf_counter() - start
        status = "OK  " if ok else "GAGAL"
        print(f"[{status}] {path.relative_to(ROOT)} ({elapsed:.1f}s){'' if ok else ' — ' + msg}")
        if not ok:
            failures.append(path)

    print(f"\n{len(notebooks) - len(failures)}/{len(notebooks)} notebook lulus.")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
