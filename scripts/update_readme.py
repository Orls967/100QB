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
DAY_FILE_RE = re.compile(r"^day(\d{3})([a-z])?_(quant|sim)_(.+)\.ipynb$")


def _session_rank(session: str | None) -> int:
    if not session:
        return 0
    return ord(session.lower()) - ord("a")


def collect_days(days_dir: Path = DAYS_DIR) -> list[dict]:
    """Kumpulkan metadata semua notebook harian, urut berdasarkan (day, session)."""
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
        day = int(meta.get("day") or m.group(1))
        raw_sess = meta.get("session") if meta.get("session") is not None else (m.group(2) or "")
        session = "" if raw_sess in ("", "a") else str(raw_sess)
        if "topic_no" in meta:
            raw_tno = meta["topic_no"]
            topic_no = None if raw_tno in (None, "bonus") else int(raw_tno)
        else:
            topic_no = day
        rows.append(
            {
                "day": day,
                "session": session,
                "topic_no": topic_no,
                "type": meta.get("type") or m.group(3),
                "topic": meta.get("topic") or m.group(4).replace("-", " ").title(),
                "date": meta.get("date") or "-",
                "file": path.relative_to(ROOT).as_posix(),
            }
        )
    return sorted(rows, key=lambda r: (r["day"], _session_rank(r["session"])))


def render_table(rows: list[dict], total: int = TOTAL_DAYS) -> str:
    """Render ringkasan + tabel markdown progres (string kosong bila belum ada notebook)."""
    if not rows:
        return ""
    n_quant = sum(r["type"] == "quant" for r in rows)
    n_sim = sum(r["type"] == "sim" for r in rows)
    done = len(rows)
    pct = round(100 * done / total) if total else 0
    lines = [
        f"**Progres: {done}/{total} topik ({pct}%)** · quant: {n_quant} · sim: {n_sim}",
        "",
        "| Hari | Topik No | Jenis | Topik | Tanggal | Notebook |",
        "|:----:|:--------:|:-----:|-------|---------|----------|",
    ]
    for r in rows:
        topic = str(r["topic"]).replace("|", "\\|")
        day_label = f"{r['day']}{r.get('session', '')}"
        raw_tno = r.get("topic_no") if "topic_no" in r else r["day"]
        topic_no = "bonus" if raw_tno is None else f"{int(raw_tno):03d}"
        lines.append(
            f"| {day_label} | {topic_no} | `{r['type']}` | {topic} | {r['date']} | "
            f"[{Path(r['file']).name}]({r['file']}) |"
        )
    return "\n".join(lines)


def replace_block(text: str, block: str) -> str:
    """Ganti isi di antara marker START/END dengan ``block``."""
    pattern = re.compile(re.escape(START) + r".*?" + re.escape(END), re.DOTALL)
    if not pattern.search(text):
        raise ValueError(f"Marker {START} / {END} tidak ditemukan di README.md")
    body = f"\n{block}\n" if block else "\n"
    return pattern.sub(lambda _: f"{START}{body}{END}", text, count=1)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Regenerasi tabel progres di README.md.")
    p.add_argument("--readme", type=Path, default=README, help="path README (default: README.md)")
    p.add_argument("--check", action="store_true", help="exit 1 bila README belum up to date")
    p.add_argument("--dry-run", action="store_true", help="cetak tabel tanpa menulis file")
    args = p.parse_args(argv)

    rows = collect_days()
    block = render_table(rows)
    if args.dry_run:
        print(block)
        return 0

    old = args.readme.read_text(encoding="utf-8")
    try:
        new = replace_block(old, block)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 2

    if args.check:
        if new != old:
            print("README.md belum up to date. Jalankan: python scripts/update_readme.py")
            return 1
        print("README.md sudah up to date.")
        return 0

    if new != old:
        args.readme.write_text(new, encoding="utf-8")
        print(f"README.md diperbarui ({len(rows)} notebook).")
    else:
        print(f"README.md sudah up to date ({len(rows)} notebook).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
