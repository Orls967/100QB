# Changelog

Semua perubahan penting pada project **100QB** dicatat di sini, per hari.
Format mengikuti [Keep a Changelog](https://keepachangelog.com/id/1.1.0/).

## [Hari 001] — 2026-11-19 · Return harian & kumulatif

### Ditambahkan
- Aturan dan dukungan multi-sesi per hari (`day001`, `day001b`, ...) berbasis tanggal WITA (`Asia/Makassar`) pada `scripts/new_day.py`, `scripts/update_readme.py`, `scripts/run_all.py`, dan `templates/day_template.ipynb`.
- Fungsi `cumulative_returns` di `src/qb100/finance.py` beserta unit test di `tests/test_finance.py`.
- Notebook `days/day001_quant_return-harian-kumulatif.ipynb` (topik 001) beserta grafik `outputs/figures/day001_return-harian-kumulatif.png`.

### Diperbarui
- `README.md` (tabel progres dengan kolom Hari & Topik No, bagian scripts, dan konvensi penamaan).
- `ROADMAP.md` (centang topik 001).

## [Hari 000] — 2026-11-19 · Setup baseline

### Ditambahkan
- Struktur repo: `days/`, `templates/`, `src/qb100/`, `scripts/`, `tests/`, `data/`, `outputs/figures/`.
- Package `qb100`: `plotting` (style seragam, `save_figure`), `data` (yfinance + cache CSV),
  `finance` (returns, log returns, volatilitas, Sharpe, drawdown), `simulation` (`get_rng`,
  random walk, Monte Carlo, ringkasan MC).
- Template notebook harian `templates/day_template.ipynb` dengan metadata `qb100`.
- Scripts `new_day.py`, `update_readme.py`, `run_all.py`.
- Pytest untuk `finance`, `simulation`, `data`/`plotting`, dan scripts.
- `ROADMAP.md` berisi 100 topik (50 quant, 50 sim) dalam 4 fase.
- `README.md`, `LICENSE` (MIT), `requirements.txt`, `pyproject.toml`, `.env.example`.
