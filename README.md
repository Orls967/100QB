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

Daftar lengkap 100 topik (dari dasar ke lanjut) ada di **[ROADMAP.md](ROADMAP.md)**.

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
belum dicentang di ROADMAP.md diambil otomatis.

```bash
python scripts/new_day.py --type quant --topic "Return harian & kumulatif"
python scripts/new_day.py --type sim --topic "Estimasi pi Monte Carlo" --dry-run   # cek tanpa membuat file
python scripts/new_day.py                                                           # ikuti ROADMAP.md
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

## 📁 Struktur folder

```text
100QB/
├── README.md              # dokumentasi ini
├── ROADMAP.md             # daftar 100 topik (quant & sim)
├── CHANGELOG.md           # catatan perubahan per hari
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
│   └── run_all.py         # eksekusi semua notebook
├── tests/                 # pytest untuk src/ dan scripts/
├── data/                  # data lokal (raw/ & cache/ tidak di-commit)
└── outputs/figures/       # grafik hasil notebook
```

> Nama package Python adalah `qb100` karena identifier Python tidak boleh diawali angka.

## 🏷️ Konvensi penamaan notebook

```text
days/dayNNN[sesi]_<quant|sim>_<topik-singkat>.ipynb
```

- `NNN` = nomor hari 3 digit (`001` … `100`).
- `[sesi]` = kosong untuk sesi pertama di hari itu, lalu `b`, `c`, dst. bila ada lebih dari satu notebook pada tanggal WITA yang sama (`day001_...`, `day001b_...`).
- `<quant|sim>` = jenis notebook.
- `<topik-singkat>` = slug huruf kecil dengan tanda hubung, mis. `monte-carlo-harga-saham`.

Contoh: `days/day001_quant_return-harian-kumulatif.ipynb`, `days/day001b_sim_estimasi-pi-monte-carlo.ipynb`.

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

Tabel di bawah dibuat otomatis oleh `python scripts/update_readme.py`. Rencana lengkap: [ROADMAP.md](ROADMAP.md) · Riwayat: [CHANGELOG.md](CHANGELOG.md).

<!-- PROGRESS:START -->
<!-- PROGRESS:END -->

## 📄 Lisensi

Dirilis di bawah [lisensi MIT](LICENSE).
