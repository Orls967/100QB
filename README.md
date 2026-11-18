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

Nomor hari otomatis (nomor terbesar di `days/` + 1). Bila `--type`/`--topic` tidak diisi,
jenis dan topik diambil dari ROADMAP.md.

```bash
python scripts/new_day.py --type quant --topic "Return harian & kumulatif"
python scripts/new_day.py --type sim --topic "Estimasi pi Monte Carlo" --dry-run   # cek tanpa membuat file
python scripts/new_day.py                                                           # ikuti ROADMAP.md
```

Opsi lain: `--slug` (nama file singkat), `--day` (paksa nomor), `--date YYYY-MM-DD`, `--force`.

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
