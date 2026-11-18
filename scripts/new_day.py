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
from zoneinfo import ZoneInfo

import nbformat

ROOT = Path(__file__).resolve().parents[1]
DAYS_DIR = ROOT / "days"
TEMPLATE = ROOT / "templates" / "day_template.ipynb"
ROADMAP = ROOT / "ROADMAP.md"
WITA_TZ = ZoneInfo("Asia/Makassar")


def wita_today(now: dt.datetime | None = None) -> str:
    """Tanggal hari ini dalam zona waktu WITA (Asia/Makassar) format YYYY-MM-DD."""
    current = now.astimezone(WITA_TZ) if now is not None else dt.datetime.now(WITA_TZ)
    return current.date().isoformat()

DAY_FILE_RE = re.compile(r"^day(\d{3})([a-z])?_(quant|sim)_.+\.ipynb$")
ROADMAP_RE = re.compile(
    r"^- \[([ xX])\] \*\*(\d{3})\*\* · `(quant|sim)` · \*\*(.+?)\*\*(?: — (.*))?$"
)
TYPES = ("quant", "sim")


def slugify(text: str, max_len: int = 40) -> str:
    """Ubah topik menjadi slug pendek: huruf kecil, angka, dan tanda hubung."""
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    text = re.sub(r"[^a-zA-Z0-9]+", "-", text.lower()).strip("-")
    if len(text) > max_len:
        text = text[:max_len].rsplit("-", 1)[0] or text[:max_len]
    return text or "topik"


def session_rank(session: str | None) -> int:
    """Urutan numerik huruf sesi: '' atau 'a' -> 0, 'b' -> 1, dst."""
    if not session:
        return 0
    return ord(session.lower()) - ord("a")


def next_session_letter(session: str | None) -> str:
    """Huruf sesi berikutnya: ''/'a' -> 'b', 'b' -> 'c', dst."""
    rank = session_rank(session)
    return chr(ord("a") + rank + 1)


def existing_sessions(days_dir: Path = DAYS_DIR) -> set[tuple[int, str]]:
    """Kumpulan (day, session) yang sudah ada di folder days/."""
    if not days_dir.exists():
        return set()
    found: set[tuple[int, str]] = set()
    for p in days_dir.iterdir():
        if m := DAY_FILE_RE.match(p.name):
            sess = "" if (m.group(2) in (None, "a")) else m.group(2)
            found.add((int(m.group(1)), sess))
    return found


def latest_notebook_meta(days_dir: Path = DAYS_DIR) -> dict | None:
    """Ambil metadata notebook terakhir di days/ berdasarkan (day, session)."""
    if not days_dir.exists():
        return None
    entries: list[dict] = []
    for p in days_dir.iterdir():
        m = DAY_FILE_RE.match(p.name)
        if not m:
            continue
        meta = {}
        try:
            meta = nbformat.read(p, as_version=4).metadata.get("qb100", {}) or {}
        except Exception:
            pass
        day = int(meta.get("day") or m.group(1))
        raw_sess = meta.get("session") if meta.get("session") is not None else (m.group(2) or "")
        sess = "" if raw_sess in ("", "a") else str(raw_sess)
        entries.append(
            {
                "day": day,
                "session": sess,
                "date": meta.get("date") or "",
                "topic_no": meta.get("topic_no"),
                "path": p,
            }
        )
    if not entries:
        return None
    entries.sort(key=lambda e: (e["day"], session_rank(e["session"])))
    return entries[-1]


def existing_days(days_dir: Path = DAYS_DIR) -> list[int]:
    """Daftar nomor hari yang sudah ada di folder days/."""
    if not days_dir.exists():
        return []
    return sorted({int(m.group(1)) for p in days_dir.iterdir() if (m := DAY_FILE_RE.match(p.name))})


def resolve_day_and_session(
    days_dir: Path = DAYS_DIR, today: str | None = None
) -> tuple[int, str]:
    """Tentukan (nomor_hari, huruf_sesi) dari notebook terakhir & tanggal WITA."""
    today = today or wita_today()
    last = latest_notebook_meta(days_dir)
    if last is None:
        return 1, ""
    if last["date"] == today:
        return last["day"], next_session_letter(last["session"])
    return last["day"] + 1, ""


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


def build_notebook(day: int, nb_type: str, topic: str, date: str, slug: str):
    """Isi template dengan metadata hari dan kembalikan objek notebook."""
    nb = nbformat.read(TEMPLATE, as_version=4)
    repl = {
        "{{DAY}}": f"{day:03d}",
        "{{DATE}}": date,
        "{{TYPE}}": nb_type,
        "{{TOPIC}}": topic,
        "{{SLUG}}": slug,
    }
    for cell in nb.cells:
        for key, val in repl.items():
            cell.source = cell.source.replace(key, val)
    nb.metadata["qb100"] = {
        "template": False,
        "day": day,
        "date": date,
        "type": nb_type,
        "topic": topic,
        "slug": slug,
    }
    return nb


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Buat notebook hari berikutnya dari template.")
    p.add_argument("--type", choices=TYPES, help="jenis notebook: quant atau sim")
    p.add_argument("--topic", help="topik hari ini (default: ambil dari ROADMAP.md)")
    p.add_argument("--slug", help="nama singkat untuk file (default: dari topik)")
    p.add_argument("--day", type=int, help="paksa nomor hari (default: otomatis)")
    p.add_argument("--date", default=dt.date.today().isoformat(), help="tanggal YYYY-MM-DD")
    p.add_argument("--dry-run", action="store_true", help="tampilkan rencana tanpa membuat file")
    p.add_argument("--force", action="store_true", help="timpa file bila sudah ada")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    day = args.day or next_day_number()
    if not 1 <= day <= 999:
        print(f"Nomor hari tidak valid: {day}", file=sys.stderr)
        return 2
    if day in existing_days() and not args.force:
        print(f"Hari {day:03d} sudah ada di days/ (pakai --force untuk menimpa).", file=sys.stderr)
        return 1

    nb_type, topic = args.type, args.topic
    if nb_type is None or topic is None:
        entry = roadmap_entry(day)
        if entry is None:
            print("Butuh --type dan --topic (topik tidak ditemukan di ROADMAP.md).", file=sys.stderr)
            return 2
        nb_type = nb_type or entry[0]
        topic = topic or entry[1]

    slug = slugify(args.slug or topic)
    path = DAYS_DIR / f"day{day:03d}_{nb_type}_{slug}.ipynb"
    nb = build_notebook(day, nb_type, topic, args.date, slug)

    if args.dry_run:
        print("[dry-run] Tidak ada file yang dibuat.")
        print(f"  hari   : {day:03d}")
        print(f"  jenis  : {nb_type}")
        print(f"  topik  : {topic}")
        print(f"  tanggal: {args.date}")
        print(f"  file   : {path.relative_to(ROOT)}")
        print(f"  sel    : {len(nb.cells)}")
        return 0

    DAYS_DIR.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, path)
    print(f"Notebook dibuat: {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
