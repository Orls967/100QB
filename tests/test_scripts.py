"""Tes ringan untuk scripts/ (tanpa membuat file di days/)."""

import datetime as dt
import sys
from pathlib import Path
from zoneinfo import ZoneInfo

import nbformat
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import new_day  # noqa: E402
import update_readme  # noqa: E402


def test_slugify():
    assert new_day.slugify("Monte Carlo Harga Saham (GBM)") == "monte-carlo-harga-saham-gbm"
    assert new_day.slugify("Estimasi π") == "estimasi"
    assert len(new_day.slugify("a" * 30 + " " + "b" * 30)) <= 40


def test_next_day_number(tmp_path):
    assert new_day.next_day_number(tmp_path) == 1
    (tmp_path / "day003_sim_x.ipynb").touch()
    (tmp_path / "notes.txt").touch()
    assert new_day.next_day_number(tmp_path) == 4


def test_build_notebook_fills_metadata_and_placeholders():
    nb = new_day.build_notebook(7, "quant", "GBM", "2026-01-01", "gbm")
    assert nb.metadata["qb100"]["day"] == 7
    assert nb.metadata["qb100"]["type"] == "quant"
    assert all("{{" not in c.source for c in nb.cells)
    assert "Hari 007" in nb.cells[0].source


def test_roadmap_has_100_parseable_entries():
    for day in (1, 50, 100):
        entry = new_day.roadmap_entry(day)
        assert entry is not None and entry[0] in ("quant", "sim")


def test_dry_run_creates_no_file(capsys):
    before = set(new_day.DAYS_DIR.iterdir())
    assert new_day.main(["--type", "sim", "--topic", "Coba", "--dry-run"]) == 0
    assert set(new_day.DAYS_DIR.iterdir()) == before
    assert "dry-run" in capsys.readouterr().out


def test_update_readme_collect_and_render(tmp_path, monkeypatch):
    nb = new_day.build_notebook(1, "sim", "Pi | MC", "2026-01-01", "pi")
    nbformat.write(nb, tmp_path / "day001_sim_pi.ipynb")
    monkeypatch.setattr(update_readme, "ROOT", tmp_path.parent)
    rows = update_readme.collect_days(tmp_path)
    assert rows[0]["day"] == 1 and rows[0]["type"] == "sim"
    table = update_readme.render_table(rows)
    assert "1/100" in table and "Pi \\| MC" in table


def test_update_readme_replace_block():
    text = f"x\n{update_readme.START}\nold\n{update_readme.END}\ny"
    out = update_readme.replace_block(text, "new")
    assert "new" in out and "old" not in out
    with pytest.raises(ValueError):
        update_readme.replace_block("tanpa marker", "new")


def test_readme_progress_is_up_to_date():
    assert update_readme.main(["--check"]) == 0


def test_resolve_day_and_session_empty_dir(tmp_path):
    assert new_day.resolve_day_and_session(tmp_path, today="2026-11-19") == (1, "")
