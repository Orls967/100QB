# 100QB

![Python](https://img.shields.io/badge/python-3.11%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)
![Progres](https://img.shields.io/badge/progres-100%20hari-orange)

**100QB** adalah project belajar Python iseng: gabungan **100 Days of Quant** dan **Kebun Simulasi**.
Setiap hari dibuat **satu Jupyter notebook kecil yang selesai**, tentang satu topik keuangan
kuantitatif atau satu simulasi seru.

## 💡 Konsep

Setiap notebook harian masuk ke salah satu dari dua jenis:

| Jenis | Isi | Contoh topik |
|-------|-----|--------------|
| `quant` | Konsep keuangan/kuantitatif, memakai data nyata (yfinance) atau data simulasi | return, volatilitas, Sharpe, Monte Carlo harga saham, VaR, GBM, Black-Scholes, korelasi, backtest sederhana |
| `sim` | Simulasi atau perhitungan iseng ("Kebun Simulasi") | estimasi π Monte Carlo, Game of Life, fraktal Mandelbrot, random walk, antrean, gerak peluru, peluruhan radioaktif |

Daftar lengkap 100 topik (dari dasar ke lanjut) ada di **[docs/ROADMAP.md](docs/ROADMAP.md)**.

## ⚙️ Setup

Butuh **Python 3.11+**.

```bash
git clone https://github.com/Orls967/100QB.git
cd 100QB

# buat & aktifkan virtualenv
python3.11 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

# install dependensi + package lokal qb100 (mode editable)
pip install -r requirements.txt
pip install -e .

# (opsional) salin contoh variabel environment
cp .env.example .env

# cek semuanya berjalan
pytest
```

> Tidak punya Python 3.11? Bisa pakai [uv](https://docs.astral.sh/uv/):
> `uv venv --python 3.11 .venv && uv pip install -r requirements.txt -e .`

Lalu buka notebook dengan `jupyter lab` atau `jupyter notebook`.

## 🛠️ Cara pakai scripts

### `scripts/new_day.py` — buat notebook hari berikutnya

Nomor hari dan huruf sesi (`b`, `c`, dst.) ditentukan otomatis dari metadata notebook terakhir
dan tanggal hari ini (zona waktu WITA). Bila masih di tanggal yang sama, nomor hari tetap sama
dan sesi naik ke huruf berikutnya. Bila `--type`/`--topic` tidak diisi, topik berikutnya yang
belum dicentang di docs/ROADMAP.md diambil otomatis.

```bash
python scripts/new_day.py --type quant --topic "Return harian & kumulatif"
python scripts/new_day.py --type sim --topic "Estimasi pi Monte Carlo" --dry-run   # cek tanpa membuat file
python scripts/new_day.py                                                           # ikuti docs/ROADMAP.md
```

Opsi lain: `--slug` (nama file singkat), `--day` (paksa nomor), `--session` (paksa huruf sesi), `--topic-no` (paksa nomor topik), `--date YYYY-MM-DD`, `--force`.

### `scripts/update_readme.py` — perbarui tabel progres

Membaca metadata semua notebook di `days/` lalu menulis ulang tabel di bagian **Progres** README ini.

```bash
python scripts/update_readme.py            # tulis ke README.md
python scripts/update_readme.py --dry-run  # tampilkan tabel saja
python scripts/update_readme.py --check    # exit 1 bila README belum up to date
```

### `scripts/run_all.py` — jalankan semua notebook

Mengeksekusi setiap notebook dari atas ke bawah dan melaporkan yang gagal.

```bash
python scripts/run_all.py              # semua notebook di days/
python scripts/run_all.py --day 7      # hanya hari 007
python scripts/run_all.py --template   # ikut uji template
python scripts/run_all.py --inplace    # simpan output ke notebook
```

### `scripts/check_docs.py` — periksa tautan dan struktur dokumen

```bash
python scripts/check_docs.py           # verifikasi tautan relatif, path absolut, dan struktur docs/project-explanations/
```

## 📁 Struktur folder

```text
100QB/
├── README.md              # dokumentasi ini
├── docs/                  # ROADMAP, CHANGELOG, dan penjelasan proyek
├── projects/              # implementasi bonus per bahasa (C, C++, Go, Ruby, Rust, TypeScript)
├── LICENSE                # MIT
├── requirements.txt       # dependensi Python
├── pyproject.toml         # metadata package qb100 + konfigurasi pytest
├── .env.example           # contoh variabel environment
├── days/                  # notebook harian (dayNNN_<jenis>_<topik>.ipynb)
├── templates/
│   └── day_template.ipynb # template notebook harian
├── src/qb100/             # package helper bersama
│   ├── plotting.py        # style grafik seragam + save_figure
│   ├── data.py            # ambil & cache harga yfinance ke data/cache
│   ├── finance.py         # returns, log returns, volatilitas, Sharpe, drawdown
│   └── simulation.py      # seed (get_rng), random walk, Monte Carlo generik
├── scripts/
│   ├── new_day.py         # buat notebook hari berikutnya
│   ├── update_readme.py   # perbarui tabel progres README
│   ├── run_all.py         # eksekusi semua notebook
│   └── check_docs.py      # periksa tautan & struktur dokumen penjelasan
├── tests/                 # pytest untuk src/ dan scripts/
├── data/                  # data lokal (raw/ & cache/ tidak di-commit)
└── outputs/               # hasil proyek dan grafik
```

> Nama package Python adalah `qb100` karena identifier Python tidak boleh diawali angka.

## 🏷️ Konvensi penamaan notebook

```text
days/dayNNN<sesi>_<quant|sim>_<topik-singkat>.ipynb
```

- `NNN` = nomor hari 3 digit (`001` … `100`).
- `[sesi]` = `a` untuk sesi pertama di hari itu, lalu `b`, `c`, dst. bila ada lebih dari satu notebook pada tanggal WITA yang sama (`day001a_...`, `day001b_...`).
- `<quant|sim>` = jenis notebook.
- `<topik-singkat>` = slug huruf kecil dengan tanda hubung, mis. `monte-carlo-harga-saham`.

Contoh: `days/day001a_quant_return-harian-kumulatif.ipynb`, `days/day001b_sim_estimasi-pi-monte-carlo.ipynb`.

Grafik yang disimpan memakai awalan yang sama, mis. `outputs/figures/day027_gbm-paths.png`.

## 📓 Isi standar tiap notebook

Setiap notebook (dibuat dari `templates/day_template.ipynb`) wajib memiliki:

1. **Judul + metadata** — nomor hari, tanggal, jenis, dan topik (juga disimpan di metadata JSON notebook `metadata.qb100`).
2. **Konsep** — penjelasan singkat dalam bahasa Indonesia.
3. **Setup** — import library dan helper dari package `qb100`.
4. **Implementasi** — kode Python utama.
5. **Visualisasi** — minimal satu grafik.
6. **Kesimpulan (3 poin)** — tiga hal yang dipelajari.
7. **Referensi** — sumber bacaan.

## 📈 Progres

Tabel di bawah dibuat otomatis oleh `python scripts/update_readme.py`. Rencana lengkap: [docs/ROADMAP.md](docs/ROADMAP.md) · Riwayat: [docs/CHANGELOG.md](docs/CHANGELOG.md).

<!-- PROGRESS:START -->
**Progres: 56/100 topik (56%)** · quant: 45 · sim: 11

| Hari | Topik No | Jenis | Bahasa | Topik | Tanggal | Notebook / Kode |
|:----:|:--------:|:-----:|:------:|-------|---------|-----------------|
| 1a | 001 | `quant` | `python` | Return harian & kumulatif | 2026-11-19 | [day001a_quant_return-harian-kumulatif.ipynb](days/day001a_quant_return-harian-kumulatif.ipynb) |
| 1b | 002 | `sim` | `python` | Estimasi π dengan Monte Carlo | 2026-11-19 | [day001b_sim_estimasi-dengan-monte-carlo.ipynb](days/day001b_sim_estimasi-dengan-monte-carlo.ipynb) |
| 1c | bonus | `quant` | `python` | Opsi Asia: Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001c_quant_opsi-asia-monte-carlo-variance.ipynb](days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb) |
| 1d | bonus | `quant` | `rust` | Opsi Asia (Rust): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001d_asian_options](projects/rust/day001d_asian_options) |
| 1e | bonus | `quant` | `go` | Opsi Asia (Go): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001e_asian_options](projects/go/day001e_asian_options) |
| 1f | bonus | `quant` | `typescript` | Opsi Asia (TypeScript): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001f_asian_options](projects/typescript/day001f_asian_options) |
| 1g | bonus | `quant` | `cpp` | Opsi Asia (C++): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001g_asian_options](projects/cpp/day001g_asian_options) |
| 1h | bonus | `quant` | `c` | Opsi Asia (C): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001h_asian_options](projects/c/day001h_asian_options) |
| 1i | bonus | `quant` | `ruby` | Opsi Asia (Ruby): Monte Carlo, variance reduction & Greeks | 2026-11-19 | [day001i_asian_options](projects/ruby/day001i_asian_options) |
| 2b | bonus | `quant` | `python` | Estimasi kovarians & portofolio varians minimum: Wishart, Marchenko-Pastur, Ledoit-Wolf | 2026-11-20 | [day002b_quant_estimasi-kovarians-portofolio-varians.ipynb](days/day002b_quant_estimasi-kovarians-portofolio-varians.ipynb) |
| 2c | bonus | `quant` | `python` | HMM dua rezim: forward-backward, Baum-Welch & Viterbi | 2026-11-20 | [day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb](days/day002c_quant_hmm-dua-rezim-forward-backward-baum.ipynb) |
| 2d | bonus | `quant` | `python` | Risiko kredit Vasicek: distribusi kerugian eksak, LHP & importance sampling | 2026-11-20 | [day002d_quant_risiko-kredit-vasicek-distribusi.ipynb](days/day002d_quant_risiko-kredit-vasicek-distribusi.ipynb) |
| 2e | bonus | `quant` | `python` | Model Heston: pricing Fourier, kalibrasi & Monte Carlo QE | 2026-11-20 | [day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb](days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb) |
| 2f | bonus | `quant` | `rust` | Inferensi Bayesian state-space: Kalman, particle filter & PMMH | 2026-11-20 | [day002f_state_space_pmmh](projects/rust/day002f_state_space_pmmh) |
| 2g | bonus | `quant` | `cpp` | AAD untuk Greeks Monte Carlo: tape mode-balik, pathwise & keranjang 10 aset | 2026-11-20 | [day002g_aad_greeks](projects/cpp/day002g_aad_greeks) |
| 2h | bonus | `quant` | `cpp` | Hawkes + limit order book: simulasi eksak, MLE, mesin pencocokan & fakta stilistik | 2026-11-20 | [day002h_hawkes_orderbook](projects/cpp/day002h_hawkes_orderbook) |
| 2i | bonus | `quant` | `cpp` | Hull-White: pohon trinomial, Jamshidian, kalibrasi & Bermudan | 2026-11-20 | [day002i_hull_white](projects/cpp/day002i_hull_white) |
| 2j | bonus | `quant` | `ruby` | Longevity risk: Lee–Carter, anuitas stokastik & modal solvabilitas | 2026-11-20 | [day002j_longevity_risk](projects/ruby/day002j_longevity_risk) |
| 2k | bonus | `quant` | `ruby` | Opsi Put Amerika (Ruby): CRR Binomial, LSM Monte Carlo & PSOR PDE | 2026-11-20 | [day002k_american_put](projects/ruby/day002k_american_put) |
| 3a | bonus | `quant` | `python` | Volatilitas lokal Dupire dari permukaan SSVI | 2026-11-21 | [day003a_quant_volatilitas-lokal-dupire-ssvi.ipynb](days/day003a_quant_volatilitas-lokal-dupire-ssvi.ipynb) |
| 3b | bonus | `quant` | `rust` | Eksekusi optimal: Almgren–Chriss, dampak transien, likuiditas stokastik | 2026-11-21 | [day003b_optimal_execution](projects/rust/day003b_optimal_execution) |
| 3c | bonus | `quant` | `go` | Risiko sistemik jaringan: kliring Eisenberg–Noe, kaskade Watts & biaya kebangkrutan | 2026-11-21 | [day003c_systemic_risk](projects/go/day003c_systemic_risk) |
| 3d | bonus | `quant` | `typescript` | Risiko ekor: GARCH–EVT, VaR/ES dan backtesting | 2026-11-21 | [day003d_tail_risk](projects/typescript/day003d_tail_risk) |
| 3e | bonus | `quant` | `go` | CVA swap suku bunga: Hull–White, intensitas CIR dan wrong-way risk | 2026-11-21 | [day003e_cva_wrong_way](projects/go/day003e_cva_wrong_way) |
| 3f | bonus | `quant` | `go` | Teori ruin Cramér–Lundberg dan simulasi kejadian langka | 2026-11-21 | [day003f_ruin_theory](projects/go/day003f_ruin_theory) |
| 4a | bonus | `quant` | `go` | Struktur suku bunga: Nelson–Siegel, AFNS dan filter Kalman | 2026-11-22 | [day004a_term_structure](projects/go/day004a_term_structure) |
| 4b | bonus | `quant` | `go` | AFNS — bias sampel terbatas, rentang kalender (long-span), dan frekuensi observasi (infill) | 2026-11-22 | [day004b_afns_inference](projects/go/day004b_afns_inference) |
| 4c | bonus | `quant` | `go` | AFNS — koreksi bias mean-reversion dan inferensi sampel terbatas dengan jackknife serta bootstrap parametrik | 2026-11-22 | [day004c_afns_bias_correction](projects/go/day004c_afns_bias_correction) |
| 4d | bonus | `quant` | `go` | AFNS — interval bootstrap yang peka bias dan batas parameter | 2026-11-22 | [day004d_afns_interval_calibration](projects/go/day004d_afns_interval_calibration) |
| 4e | bonus | `quant` | `go` | Harga listrik negatif dan arbitrase baterai (DP stokastik) | 2026-11-22 | [day004e_power_storage](projects/go/day004e_power_storage) |
| 4f | bonus | `quant` | `rust` | Port lengkap AFNS term structure, kalibrasi interval, harga listrik negatif, dan valuasi proyek baterai terpadu | 2026-11-22 | [day004f_energy_transition](projects/rust/day004f_energy_transition) |
| 4g | bonus | `quant` | `rust` | Market making optimal: Guéant–Lehalle–Tapia, HJB dan salah spesifikasi | 2026-11-22 | [day004g_market_making](projects/rust/day004g_market_making) |
| 4h | bonus | `quant` | `rust` | Volatilitas kasar: rough Bergomi, skema hybrid, skew ATM dan identifiabilitas Hurst | 2026-11-22 | [day004h_rough_bergomi](projects/rust/day004h_rough_bergomi) |
| 4i | bonus | `quant` | `rust` | Cointegration dan trading OU optimal: Dickey-Fuller, Engle-Granger, Johansen dan ambang Bertram | 2026-11-22 | [day004i_pairs_trading](projects/rust/day004i_pairs_trading) |
| 4j | bonus | `quant` | `rust` | Deep hedging dari nol: jaringan saraf, risiko entropik, biaya transaksi dan pasar tak lengkap | 2026-11-22 | [day004j_deep_hedging](projects/rust/day004j_deep_hedging) |
| 4k | bonus | `quant` | `rust` | AMM: impermanent loss, LVR, biaya dan likuiditas terkonsentrasi | 2026-11-22 | [day004k_amm_lvr](projects/rust/day004k_amm_lvr) |
| 5a | bonus | `quant` | `rust` | SABR: CEV eksak, ekspansi Hagan, kepadatan negatif, kalibrasi dan delta Bartlett | 2026-11-23 | [day005a_sabr](projects/rust/day005a_sabr) |
| 5b | bonus | `quant` | `rust` | Harga opsi Eropa Black-Scholes & put-call parity | 2026-11-23 | [day005b_black_scholes](projects/rust/day005b_black_scholes) |
| 5c | bonus | `quant` | `rust` | Kriteria Kelly dan risiko drawdown: peluang sentuh eksak, Kelly fraksional dan galat estimasi | 2026-11-23 | [day005c_kelly_drawdown](projects/rust/day005c_kelly_drawdown) |
| 5d | bonus | `quant` | `rust` | Jump-diffusion Merton: seri eksak, PIDE Crank-Nicolson, kalibrasi smile dan risiko ekor | 2026-11-23 | [day005d_merton_jump_diffusion](projects/rust/day005d_merton_jump_diffusion) |
| 5e | bonus | `quant` | `ruby` | Kopula, ketergantungan ekor dan alokasi modal Euler | 2026-11-23 | [day005e_copula_tail_dependence](projects/ruby/day005e_copula_tail_dependence) |
| 5f | bonus | `quant` | `ruby` | Valuasi tranche CDO: LHP Vasicek, rekursi ASB, base correlation & exponential twisting | 2026-11-23 | [day005f_cdo_tranche](projects/ruby/day005f_cdo_tranche) |
| 5g | bonus | `quant` | `ruby` | Multilevel Monte Carlo (Giles) untuk CIR dan Heston: laju kuat, hukum kompleksitas & antithetic | 2026-11-23 | [day005g_mlmc_heston](projects/ruby/day005g_mlmc_heston) |
| 5h | bonus | `quant` | `ruby` | Backtest dan overfitting: Sharpe, bias seleksi, pengujian ganda dan PBO | 2026-11-23 | [day005h_backtest_overfitting](projects/ruby/day005h_backtest_overfitting) |
| 5i | bonus | `quant` | `ruby` | Keluarga GARCH, QMLE Sandwich, MLE Student-t, VaR/ES, DCC, HAR-RV & salah spesifikasi | 2026-11-23 | [day005i_garch_family](projects/ruby/day005i_garch_family) |
| 5j | bonus | `quant` | `ruby` | Model Faktor, PCA BBP, GRS, Fama-MacBeth Shanken, Risk Parity, Black-Litterman, Markowitz QP & Backtest Portofolio | 2026-11-23 | [day005j_factor_portfolio_lab](projects/ruby/day005j_factor_portfolio_lab) |
| 6a | bonus | `sim` | `python+rust+go` | Self-Avoiding Walk (SAW) Poliglot: Enumerasi Eksak, Eksponen Kritis & Pivot Monte Carlo | 2026-11-24 | [day006a_saw_polyglot](projects/rust/day006a_saw_polyglot) |
| 6b | bonus | `sim` | `python+rust+go` | Self-Avoiding Walk (SAW) Poliglot: Enumerasi Eksak, Pivot Monte Carlo & Eksponen Kritis | 2026-11-24 | [day006b_saw_polyglot](projects/rust/day006b_saw_polyglot) |
| 6c | bonus | `sim` | `python+rust+go` | Perkolasi Bond Kisi Persegi Poliglot: Newman-Ziff, Dualitas-Diri Eksak & Eksponen Kritis | 2026-11-24 | [day006c_percolation_polyglot](projects/rust/day006c_percolation_polyglot) |
| 6d | bonus | `sim` | `cpp` | Solver PDE Heston 2D ADI Craig–Sneyd (C++): Tangga Optimasi Terukur & Roofline | 2026-11-24 | [day006d_heston_adi_cpp](projects/cpp/day006d_heston_adi_cpp) |
| 6e | bonus | `sim` | `python` | HMC/NUTS dari Nol untuk Volatilitas Stokastik Bayesian & Kalibrasi SBC | 2026-11-24 | [day006e_sim_bayesian-sv-nuts.ipynb](days/day006e_sim_bayesian-sv-nuts.ipynb) |
| 7a | bonus | `sim` | `cpp` | Kekacauan Deterministik (C++): Eksponen Lyapunov Peta Logistik & Lorenz, Konstanta Feigenbaum, dan Dimensi Kaplan–Yorke | 2026-11-25 | [day007a_chaos_lyapunov_cpp](projects/cpp/day007a_chaos_lyapunov_cpp) |
| 8a | bonus | `sim` | `python` | Model Ising 2D Lintas Tiga Bahasa (C, C++, Python): Metropolis, Wolff, Swendsen–Wang, FSS & Paritas Bit-Identik | 2026-11-26 | [day008a_sim_ising-trilingual.ipynb](days/day008a_sim_ising-trilingual.ipynb) |
| 8a | bonus | `sim` | `c+cpp+python` | Model Ising 2D Lintas Tiga Bahasa (C, C++, Python): Metropolis, Wolff, Swendsen–Wang, FSS & Paritas Bit-Identik | 2026-11-26 | [day008a_ising_trilingual](projects/c/day008a_ising_trilingual) |
| 8b | bonus | `sim` | `cpp` | Simulasi Wa-Tor di C++20 (Ikan vs Hiu), Fit Lotka–Volterra RK4 & Efek Ukuran Hingga | 2026-11-26 | [day008b_wator_cpp](projects/cpp/day008b_wator_cpp) |
| 8c | bonus | `sim` | `cpp` | Model Vicsek (Aktif-Materi) di C++20 dengan Cell List O(N) vs Brute Force O(N^2) | 2026-11-26 | [day008c_vicsek_cpp](projects/cpp/day008c_vicsek_cpp) |
<!-- PROGRESS:END -->

## 📄 Lisensi

Dirilis di bawah [lisensi MIT](LICENSE).
