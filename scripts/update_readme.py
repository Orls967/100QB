#!/usr/bin/env python
"""Regenerasi tabel progres di README.md dari isi folder ``days/``.

Metadata tiap notebook dibaca dari ``notebook.metadata["qb100"]`` (diisi oleh
``new_day.py``); bila tidak ada, jatuh balik ke informasi dari nama file.
Tabel ditulis di antara marker ``<!-- PROGRESS:START -->`` dan ``<!-- PROGRESS:END -->``.

Contoh:
    python scripts/update_readme.py           # tulis ke README.md
    python scripts/update_readme.py --check   # exit 1 bila README belum up to date
    python scripts/update_readme.py --dry-run # tampilkan tabel tanpa menulis
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import nbformat

ROOT = Path(__file__).resolve().parents[1]
DAYS_DIR = ROOT / "days"
README = ROOT / "README.md"
TOTAL_DAYS = 100

START = "<!-- PROGRESS:START -->"
END = "<!-- PROGRESS:END -->"
DAY_FILE_RE = re.compile(r"^day(\d{3})_(quant|sim)_(.+)\.ipynb$")


def collect_days(days_dir: Path = DAYS_DIR) -> list[dict]:
    """Kumpulkan metadata semua notebook harian, urut berdasarkan nomor hari."""
    rows = []
    if not days_dir.exists():
        return rows
    for path in sorted(days_dir.glob("day*.ipynb")):
        m = DAY_FILE_RE.match(path.name)
        if not m:
            continue
        meta = {}
        try:
            meta = nbformat.read(path, as_version=4).metadata.get("qb100", {}) or {}
        except Exception as exc:  # notebook rusak tidak boleh menggagalkan tabel
            print(f"Peringatan: gagal membaca {path.name}: {exc}", file=sys.stderr)
        rows.append(
            {
                "day": int(meta.get("day") or m.group(1)),
                "type": meta.get("type") or m.group(2),
                "topic": meta.get("topic") or m.group(3).replace("-", " ").title(),
                "date": meta.get("date") or "-",
                "file": path.relative_to(ROOT).as_posix(),
            }
        )
    return sorted(rows, key=lambda r: r["day"])
