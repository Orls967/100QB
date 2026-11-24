# day006d — Solver PDE Heston 2D ADI Craig–Sneyd & Tangga Optimasi Terukur (C++20)

Implementasi solver persamaan diferensial parsial (PDE) 2D **Heston** untuk opsi *call* Eropa menggunakan skema **Alternating Direction Implicit (ADI) Craig–Sneyd** ($\theta = 1/2$) dalam **C++20 murni** (pustaka standar saja, `899` baris total `src/` + `tests/`), dilengkapi kuadratur semi-analitik referensi (**Gil–Pelaez** & **Lewis** dengan fungsi karakteristik *little trap* Albrecher et al.) serta **tangga optimasi 6 tingkat (`V0`–`V5`)** yang diukur ketat terhadap model *roofline* Apple M4 (`arm64`).

## Intuisi Tangga Optimasi

1. **`V0` (Naif)**: `std::vector<std::vector<double>>` (`U[i][j]`), alokasi buffer kerja tiap langkah waktu, dan eliminasi Thomas dengan pembagian (`fdiv`) di dalam loop.
2. **`V1` (Kontigu + Urutan Loop)**: Array datar rata 128-byte (`j * ld + i`, `ld` kelipatan 16 `double`), prealokasi di luar loop waktu, dan anotasi `__restrict`.
3. **`V2` (Prekomputasi LU Thomas)**: Koefisien PDE tidak bergantung waktu sehingga faktorisasi LU tridiagonal (`xinv`, `xcp`, `vinv`, `vcp`) dihitung sekali; pembagian dihapus dari loop dalam (menguntungkan saat seluruh larik muat di cache L2, namun menambah kapasitas kerja `+4,99 MiB` pada grid `801×401`).
4. **`V3` (Batched Thomas Arah-$v$)**: Memproses blok kolom selebar `B` secara serempak pada tata letak *row-major* yang sudah tersisip antar-garis-$v$, mengubah akses ber-*stride* buruk menjadi akses kontigu dan mengaktifkan auto-vektorisasi NEON (`width=2, interleave=4`).
5. **`V4` (Fusi Operator & Sapuan)**: Menggabungkan evaluasi turunan campuran eksplisit + RHS + `x-solve` per baris di cache L1d serta menggabungkan pengurangan $-\frac{1}{2}A_2$ ke dalam sapuan maju `v-solve` (menurunkan lalu lintas memori model dari `280 B` ke `152 B` per sel-langkah).
6. **`V5` (`std::thread` + `std::barrier` Tanpa False Sharing)**: Partisi statis baris (tahap-$x$) dan blok 16-kolom / 128-byte (tahap-$v$) tanpa *false sharing*, menghasilkan solusi **bit-identik** dengan `V2`–`V4` untuk sembarang jumlah *thread* (`T ∈ {1, 2, 4, 10}`).

## Ringkasan Hasil (`801×401`, `M = 200`, Set A: Feller Dilanggar)

| Varian | Median ± MAD (s) | ns / sel-langkah | GFLOP/s | GB/s (model) | Speedup vs `V0` | Speedup vs Sebelumnya | Catatan |
|---|---:|---:|---:|---:|---:|---:|---|
| `V0` | `3,2755 ± 0,0282` | `50,99` | `1,04` | `4,86` | `1,00×` | `1,00×` | Baseline naif |
| `V1` | `1,7929 ± 0,0016` | `27,91` | `1,90` | `8,89` | `1,83×` | `1,83×` | `12,48 MiB` (< L2 P-core `16 MiB`) |
| `V2` | `2,0952 ± 0,0034` | `32,62` | `1,62` | `8,58` | `1,56×` | `0,86×` | **Memburuk di 801×401** (`17,48 MiB` > L2 `16 MiB`); naik `1,10–1,15×` di `≤ 401×201` |
| `V3` | `1,2098 ± 0,0053` | `18,83` | `2,81` | `14,87` | `2,71×` | `1,73×` | Batched Thomas (`B = 64`), NEON aktif |
| `V4` | `1,1602 ± 0,0064` | `18,06` | `2,93` | `8,42` | `2,82×` | `1,04×` | Fusi lintasan memori (`152 B`/sel-langkah) |
| `V5` (`T=10`) | `0,2622 ± 0,0214` | `4,08` | `12,98` | `37,23` | **`12,49×`** | `4,42×` | Bit-identik `V2`–`V4` (`T=8` puncak `0,2422 s`) |

Verifikasi numerik Set A ($S_0=K=100, T=1, r=0{,}03, \kappa=1{,}5, \theta=0{,}04, \sigma=0{,}5, \rho=-0{,}7, v_0=0{,}04$, $2\kappa\theta = 0{,}12 < \sigma^2 = 0{,}25$):
- Referensi semi-analitik Gil–Pelaez = `8,802660962859` (selisih paritas terhadap Lewis = `4,55e-13`; cocok dengan `qb100` Python `8,802660962859`).
- Galat PDE `401×201, M=200` = `2,54e-3` (`<= 5e-3`), `801×401, M=400` = `6,07e-4`; orde konvergensi empiris 3-level = `2,28` (Set A) dan `2,01` (Set B, Feller terpenuhi).

## Cara Build, Uji & Benchmark

```bash
cd projects/cpp/day006d_heston_adi_cpp

# Kompilasi & jalankan 8 unit test
make test

# Laporan auto-vektorisasi kompiler
make vecreport

# Jalankan seluruh protokol benchmark (portabel + -mcpu=native)
make bench

# Hasilkan grafik PNG+SVG dan outputs/day006d/results.json
../../../.venv/bin/python py/analyze.py
```

Penjelasan lengkap metode, penurunan rumus, analisis *roofline*, dan jejak bukti tersedia di [day006d_sim_heston-adi-cpp.md](../../../docs/project-explanations/day006d_sim_heston-adi-cpp.md).
