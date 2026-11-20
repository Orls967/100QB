"""Tes untuk scripts/check_docs.py."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import check_docs  # noqa: E402


def test_check_docs_detects_broken_link(tmp_path):
    doc_dir = tmp_path / "docs" / "project-explanations"
    doc_dir.mkdir(parents=True)
    (tmp_path / "README.md").write_text("[ok](docs/project-explanations/README.md)\n", encoding="utf-8")
    (doc_dir / "README.md").write_text("[hilang](tidak_ada.md)\n", encoding="utf-8")
    errors = check_docs.run_checks(tmp_path)
    assert any("tidak_ada.md" in e for e in errors)
    assert check_docs.main(["--root", str(tmp_path)]) == 1


def test_check_docs_passes_current_repo():
    errors = check_docs.run_checks(ROOT)
    assert errors == [], f"Pelanggaran ditemukan: {errors[:5]}"
