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


def test_resolve_day_and_session_same_date(tmp_path):
    nb1 = new_day.build_notebook(1, "quant", "Return", "2026-11-19", "return", session="", topic_no=1)
    nbformat.write(nb1, tmp_path / "day001_quant_return.ipynb")
    assert new_day.resolve_day_and_session(tmp_path, today="2026-11-19") == (1, "b")

    nb2 = new_day.build_notebook(1, "sim", "Pi", "2026-11-19", "pi", session="b", topic_no=2)
    nbformat.write(nb2, tmp_path / "day001b_sim_pi.ipynb")
    assert new_day.resolve_day_and_session(tmp_path, today="2026-11-19") == (1, "c")


def test_resolve_day_and_session_new_date(tmp_path):
    nb = new_day.build_notebook(1, "sim", "Pi", "2026-11-19", "pi", session="b", topic_no=2)
    nbformat.write(nb, tmp_path / "day001b_sim_pi.ipynb")
    assert new_day.resolve_day_and_session(tmp_path, today="2026-11-20") == (2, "")


def test_resolve_day_and_session_midnight_wita(tmp_path):
    before_midnight = dt.datetime(2026, 11, 19, 15, 59, tzinfo=ZoneInfo("UTC"))
    after_midnight = dt.datetime(2026, 11, 19, 16, 0, tzinfo=ZoneInfo("UTC"))
    date_before = new_day.wita_today(before_midnight)
    date_after = new_day.wita_today(after_midnight)
    assert date_before == "2026-11-19"
    assert date_after == "2026-11-20"

    nb = new_day.build_notebook(1, "quant", "Return", date_before, "return", topic_no=1)
    nbformat.write(nb, tmp_path / "day001_quant_return.ipynb")
    assert new_day.resolve_day_and_session(tmp_path, today=date_before) == (1, "b")
    assert new_day.resolve_day_and_session(tmp_path, today=date_after) == (2, "")


def test_next_unchecked_topic_no(tmp_path):
    roadmap = tmp_path / "ROADMAP.md"
    roadmap.write_text(
        "- [x] **001** · `quant` · **Satu** — a.\n"
        "- [ ] **002** · `sim` · **Dua** — b.\n",
        encoding="utf-8",
    )
    assert new_day.next_unchecked_topic_no(roadmap) == 2
    assert new_day.roadmap_entry(2, roadmap) == ("sim", "Dua")


def test_update_readme_collect_and_render_sessions(tmp_path, monkeypatch):
    nb1b = new_day.build_notebook(1, "sim", "Pi", "2026-11-19", "pi", session="b", topic_no=2)
    nb1 = new_day.build_notebook(1, "quant", "Return", "2026-11-19", "ret", session="", topic_no=1)
    nbformat.write(nb1b, tmp_path / "day001b_sim_pi.ipynb")
    nbformat.write(nb1, tmp_path / "day001_quant_ret.ipynb")
    monkeypatch.setattr(update_readme, "ROOT", tmp_path.parent)
    rows = update_readme.collect_days(tmp_path)
    assert [(r["day"], r["session"], r["topic_no"]) for r in rows] == [(1, "", 1), (1, "b", 2)]
    table = update_readme.render_table(rows)
    assert "| 1 | 001 | `quant` | `python` |" in table
    assert "| 1b | 002 | `sim` | `python` |" in table


def test_update_readme_bonus_topic_no(tmp_path, monkeypatch):
    nb = new_day.build_notebook(
        1, "quant", "Opsi Asia", "2026-11-19", "opsi-asia", session="c", bonus=True
    )
    assert nb.metadata["qb100"]["topic_no"] is None
    nbformat.write(nb, tmp_path / "day001c_quant_opsi-asia.ipynb")
    monkeypatch.setattr(update_readme, "ROOT", tmp_path.parent)
    rows = update_readme.collect_days(tmp_path)
    assert rows[0]["topic_no"] is None
    table = update_readme.render_table(rows)
    assert "| 1c | bonus | `quant` | `python` |" in table


def test_update_readme_collects_multi_language_meta(tmp_path, monkeypatch):
    days = tmp_path / "days"
    rust = tmp_path / "rust" / "day001d_asian"
    go_dir = tmp_path / "go" / "day001e_asian"
    ts_dir = tmp_path / "typescript" / "day001f_asian"
    rb_dir = tmp_path / "ruby" / "day001i_asian"
    days.mkdir()
    rust.mkdir(parents=True)
    go_dir.mkdir(parents=True)
    ts_dir.mkdir(parents=True)
    rb_dir.mkdir(parents=True)
    (rust / "meta.json").write_text(
        '{"day": 1, "session": "d", "topic_no": null, "type": "quant", '
        '"topic": "Opsi Asia (Rust)", "date": "2026-11-19", "link": "projects/rust/day001d_asian"}',
        encoding="utf-8",
    )
    (go_dir / "meta.json").write_text(
        '{"day": 1, "session": "e", "topic_no": null, "type": "quant", "language": "go", '
        '"topic": "Opsi Asia (Go)", "date": "2026-11-19", "link": "projects/go/day001e_asian"}',
        encoding="utf-8",
    )
    (ts_dir / "meta.json").write_text(
        '{"day": 1, "session": "f", "topic_no": null, "type": "quant", "language": "typescript", '
        '"topic": "Opsi Asia (TypeScript)", "date": "2026-11-19", "link": "projects/typescript/day001f_asian"}',
        encoding="utf-8",
    )
    (rb_dir / "meta.json").write_text(
        '{"day": 1, "session": "i", "topic_no": null, "type": "quant", "language": "ruby", '
        '"topic": "Opsi Asia (Ruby)", "date": "2026-11-19", "link": "projects/ruby/day001i_asian"}',
        encoding="utf-8",
    )
    assert "ruby" in update_readme.LANG_DIRS
    monkeypatch.setattr(update_readme, "ROOT", tmp_path)
    rows = update_readme.collect_days(
        days, lang_roots=[tmp_path / "rust", tmp_path / "go", tmp_path / "typescript", tmp_path / "ruby"]
    )
    assert len(rows) == 4
    assert [r["language"] for r in rows] == ["rust", "go", "typescript", "ruby"]
    table = update_readme.render_table(rows)
    assert "| 1d | bonus | `quant` | `rust` | Opsi Asia (Rust) | 2026-11-19 | [day001d_asian](projects/rust/day001d_asian) |" in table
    assert "| 1e | bonus | `quant` | `go` | Opsi Asia (Go) | 2026-11-19 | [day001e_asian](projects/go/day001e_asian) |" in table
    assert "| 1f | bonus | `quant` | `typescript` | Opsi Asia (TypeScript) | 2026-11-19 | [day001f_asian](projects/typescript/day001f_asian) |" in table
    assert "| 1i | bonus | `quant` | `ruby` | Opsi Asia (Ruby) | 2026-11-19 | [day001i_asian](projects/ruby/day001i_asian) |" in table
