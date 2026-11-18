# Changelog

Semua perubahan penting pada project **100QB** dicatat di sini, per hari.
Format mengikuti [Keep a Changelog](https://keepachangelog.com/id/1.1.0/).

## [Hari 001c] — 2026-11-19 · Opsi Asia: Monte Carlo, variance reduction & Greeks (Bonus)

### Ditambahkan
- Modul `src/qb100/options.py` (Black-Scholes call/put, geometric Asian call diskret Kemna–Vorst, simulasi GBM tervektorisasi, estimator Naive/Antithetic/Control Variate, serta Delta *pathwise* & *bump-and-revalue* CRN) beserta unit test di `tests/test_options.py`.
- Dukungan `--bonus` (`topic_no: null` / tampilan `"bonus"` di tabel progres) pada `scripts/new_day.py`, `scripts/update_readme.py`, dan `tests/test_scripts.py`.
- Notebook off-roadmap `days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb` beserta grafik `day001c_asian-convergence-loglog.png` dan `day001c_asian-estimators-ci.png`.

### Diperbarui
- `README.md` (tabel progres baris `1c` dengan `Topik No = bonus`).

## [Hari 001b] — 2026-11-19 · Estimasi π dengan Monte Carlo

### Ditambahkan
- Helper `running_mean` dan `running_se` di `src/qb100/simulation.py` beserta unit test di `tests/test_simulation.py`.
- Notebook `days/day001b_sim_estimasi-dengan-monte-carlo.ipynb` (topik 002) beserta 4 grafik (`day001b_pi-scatter.png`, `day001b_pi-convergence.png`, `day001b_pi-error-scaling.png`, `day001b_pi-histogram.png`) di `outputs/figures/`.

### Diperbarui
- `README.md` (tabel progres baris `1b`) dan `ROADMAP.md` (centang topik 002).

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
