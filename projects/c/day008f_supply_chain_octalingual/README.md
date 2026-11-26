# Day 8f — Jaringan Rantai Pasok Oktalingual (40 Tahap)

> **Catatan Duplikasi**: Pemeriksaan duplikasi Langkah 0 menemukan topik bursa/order-book sudah ada di `day002h` dan `day004g`, sehingga `day008f` mengimplementasikan topik cadangan resmi: **jaringan rantai pasok 3 pabrik / 4 gudang / 10 toko dengan kebijakan `(s,S)`, peramalan EWMA, pengiriman, keterlambatan, kit 4-SKU, dan biaya**.

Implementasi alur **40 tahap (`s01`–`s40`, `all`)** secara **independen di 8 bahasa** (`C11`, `C++20`, `Rust`, `Go`, `Ruby`, `Python`, `TypeScript`, dan `Jupyter Notebook`), dengan paritas bilangan bulat `int64` bit-identik pada `S01`, `S04`, `S05`, `S06`, dan `S34`.

## Cara Menjalankan
```bash
make test    # kompilasi 4 bahasa + selftest (s01) 7 CLI
make run     # jalankan seluruh 8 bahasa (40 tahap) + compare.py
```

Keluaran per bahasa ditulis ke `out/<lang>/` (`103` berkas per bahasa; log besar S10 diabaikan oleh `.gitignore`), ringkasan ke `outputs/day008f/results.json`, dan 6 pasang grafik gabungan ke `outputs/figures/day008f_*.{png,svg}`.
