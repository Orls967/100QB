#!/usr/bin/env python
"""Buat notebook hari berikutnya dari ``templates/day_template.ipynb``.

Contoh:
    python scripts/new_day.py --type quant --topic "Return harian & kumulatif"
    python scripts/new_day.py --type sim --topic "Estimasi pi Monte Carlo" --dry-run
    python scripts/new_day.py            # ambil jenis & topik dari ROADMAP.md

Nomor hari otomatis = nomor terbesar di ``days/`` + 1 (bisa dipaksa dengan ``--day``).
Nama file: ``days/dayNNN_<quant|sim>_<topik-singkat>.ipynb``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
import unicodedata
from pathlib import Path

import nbformat

ROOT = Path(__file__).resolve().parents[1]
DAYS_DIR = ROOT / "days"
TEMPLATE = ROOT / "templates" / "day_template.ipynb"
ROADMAP = ROOT / "ROADMAP.md"

DAY_FILE_RE = re.compile(r"^day(\d{3})_(quant|sim)_.+\.ipynb$")
ROADMAP_RE = re.compile(
    r"^- \[[ xX]\] \*\*(\d{3})\*\* · `(quant|sim)` · \*\*(.+?)\*\*(?: — (.*))?$"
)
TYPES = ("quant", "sim")


def slugify(text: str, max_len: int = 40) -> str:
    """Ubah topik menjadi slug pendek: huruf kecil, angka, dan tanda hubung."""
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    text = re.sub(r"[^a-zA-Z0-9]+", "-", text.lower()).strip("-")
    if len(text) > max_len:
        text = text[:max_len].rsplit("-", 1)[0] or text[:max_len]
    return text or "topik"


def existing_days(days_dir: Path = DAYS_DIR) -> list[int]:
    """Daftar nomor hari yang sudah ada di folder days/."""
    if not days_dir.exists():
        return []
    return sorted(int(m.group(1)) for p in days_dir.iterdir() if (m := DAY_FILE_RE.match(p.name)))


def next_day_number(days_dir: Path = DAYS_DIR) -> int:
    days = existing_days(days_dir)
    return (max(days) + 1) if days else 1


def roadmap_entry(day: int, roadmap: Path = ROADMAP) -> tuple[str, str] | None:
    """Cari (jenis, topik) untuk nomor hari tertentu di ROADMAP.md."""
    if not roadmap.exists():
        return None
    for line in roadmap.read_text(encoding="utf-8").splitlines():
        m = ROADMAP_RE.match(line.strip())
        if m and int(m.group(1)) == day:
            return m.group(2), m.group(3)
    return None
