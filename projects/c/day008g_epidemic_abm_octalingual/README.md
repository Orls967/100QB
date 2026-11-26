# Day 008g — Simulasi Epidemi Berbasis Agen SEIR Berlapis (39 Tahap × 8 Bahasa)

Implementasi independen alur 39 tahap (`s01`–`s39`) simulasi epidemi berbasis agen (SEIR + `Ia/Is/H/ICU/R/D` + *waning*, 3 lapisan kontak, 2 varian, intervensi berlapis, kapasitas rumah sakit/ICU, kontrol negatif, dan anti-hampa) di 8 bahasa: **C11, C++20, Rust, Go, Ruby, Python, TypeScript, dan Jupyter Notebook**.

## Cara Menjalankan
```bash
make test
./run_all.sh
../../../.venv/bin/python compare.py
```
