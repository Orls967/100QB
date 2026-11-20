#!/usr/bin/env python
"""Pemeriksa tautan relatif Markdown dan struktur dokumen penjelasan proyek 100QB."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VALID_STATUSES = {"tercapai", "sebagian", "belum tercapai", "belum diketahui"}
LINK_RE = re.compile(r"(?<!\!)\[[^\]]*\]\(([^)\s]+)\)|!\[[^\]]*\]\(([^)\s]+)\)")
SKIP_PARTS = {"build", "build-asan", "build-ubsan", "build_asan", "build_ubsan", "build_release", "node_modules", "_deps", "target", ".venv"}


def _iter_markdown_files(root: Path) -> list[Path]:
    files: list[Path] = []
    readme = root / "README.md"
    if readme.exists():
        files.append(readme)
    for base in (root / "docs", root / "projects"):
        if not base.exists():
            continue
        for p in sorted(base.rglob("*.md")):
            if p.name == "PROMPT.md":
                continue
            if any(part in SKIP_PARTS or part.startswith("build") for part in p.parts):
                continue
            if base.name == "projects" and p.name != "README.md":
                continue
            files.append(p)
    return files


def check_links(files: list[Path]) -> list[str]:
    """Periksa tautan relatif dan larangan path absolut /Users/ atau file://."""
    errors: list[str] = []
    for f in files:
        text = f.read_text(encoding="utf-8")
        if "/Users/" in text or "file://" in text:
            errors.append(f"{f}: mengandung '/Users/' atau 'file://'")
        stripped = re.sub(r"```.*?```", "", text, flags=re.S)
        for m in LINK_RE.finditer(stripped):
            target = m.group(1) or m.group(2)
            if re.match(r"^(https?:|mailto:|#)", target):
                continue
            rel_path = target.split("#", 1)[0]
            if not rel_path:
                continue
            resolved = (f.parent / rel_path).resolve()
            if not resolved.exists():
                errors.append(f"{f}: tautan rusak -> {target}")
    return errors


def check_explanation_doc(path: Path) -> list[str]:
    """Periksa 11 bagian bernomor, >=10 poin teori, dan baris status sah pada dokumen penjelasan."""
    errors: list[str] = []
    text = path.read_text(encoding="utf-8")

    for sec in range(1, 12):
        if not re.search(rf"^##\s+{sec}\.\s+", text, flags=re.M):
            errors.append(f"{path.name}: bagian '## {sec}.' tidak ditemukan")

    # Cek >= 10 poin teori di Bagian 9
    sec9_match = re.search(r"^##\s+9\.\s+.*?(?=^##\s+10\.)", text, flags=re.M | re.S)
    if not sec9_match:
        errors.append(f"{path.name}: isi Bagian 9 tidak terbaca")
    else:
        theory_items = re.findall(r"^\d+\.\s+\*\*[^*]+\*\*", sec9_match.group(0), flags=re.M)
        if len(theory_items) < 10:
            errors.append(f"{path.name}: jumlah poin landasan teori {len(theory_items)} < 10")

    # Cek baris tabel status di Bagian 7
    sec7_match = re.search(r"^##\s+7\.\s+.*?(?=^##\s+8\.)", text, flags=re.M | re.S)
    if not sec7_match:
        errors.append(f"{path.name}: isi Bagian 7 tidak terbaca")
    else:
        status_rows = re.findall(r"^\|\s*[^|]+\|\s*`([^`]+)`\s*\|", sec7_match.group(0), flags=re.M)
        if not status_rows:
            errors.append(f"{path.name}: tidak ada baris status di tabel Bagian 7")
        for st in status_rows:
            if st not in VALID_STATUSES:
                errors.append(f"{path.name}: status tidak sah '{st}'")

    return errors


def run_checks(root: Path = ROOT) -> list[str]:
    errors = check_links(_iter_markdown_files(root))
    exp_dir = root / "docs" / "project-explanations"
    if exp_dir.exists():
        for doc in sorted(exp_dir.glob("day*.md")):
            errors.extend(check_explanation_doc(doc))
    return errors


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Periksa tautan dan struktur dokumen 100QB.")
    p.add_argument("--root", type=Path, default=ROOT, help="root direktori repo")
    args = p.parse_args(argv)
    errors = run_checks(args.root)
    if errors:
        for err in errors:
            print(f"[ERROR] {err}", file=sys.stderr)
        print(f"Total pelanggaran: {len(errors)}", file=sys.stderr)
        return 1
    print("check_docs: 0 pelanggaran.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
