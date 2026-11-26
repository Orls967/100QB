# day008h — Jaringan Antrean Terminal Bandara Oktalingual (48 Tahap)

Proyek simulasi jaringan antrean terminal bandara berbasis agen (penumpang, bagasi, jadwal penerbangan, staf per jam, gangguan/perbaikan) dengan validasi eksak menggunakan rantai Markov diskret (`Geo/Geo/c/K` dan seri 2 stasiun), uji metamorfik (`S26`–`S28`), kontrol negatif (`S25`), reduksi varians CRN (`S36`), dan aturan anti-hampa (`S38`, `compare.py`).

## Menjalankan
```bash
make test        # selftest & paritas 7 CLI
./run_all.sh     # jalankan 7 bahasa CLI
../../../.venv/bin/python compare.py
```
