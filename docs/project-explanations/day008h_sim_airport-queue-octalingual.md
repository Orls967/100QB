# Hari 8h (`day008h`) — Jaringan Antrean Terminal Bandara Berbasis Agen dengan Validasi DTMC Eksak (48 Tahap, 8 Bahasa)

> **Catatan Konteks**: Dikerjakan pada 2026-11-26 (`day008h`) dalam mode santai tanpa batas waktu sebagai proyek oktalingual (`c`, `cpp`, `rust`, `go`, `ruby`, `py`, `ts`, `ipynb`) dengan inti simulasi bilangan bulat murni (`int64`/basis points `BP = 10000`), validasi terhadap solusi eksak Rantai Markov Diskret (*Discrete-Time Markov Chain*, DTMC), uji metamorfik sebagai tahap eksplisit, dan pemeriksaan anti-hampa otomatis (`vacuous_flags = []`).

---

## 1. Ringkasan Eksekutif

Proyek `day008h` membangun simulator **jaringan antrean terminal bandara berbasis agen** dengan **48 tahap identik (`s01`–`s48`, `all`)** di **8 bahasa** (`C11`, `C++20`, `Rust`, `Go`, `Ruby`, `Python`, `TypeScript`, dan `Jupyter Notebook`), menghasilkan **101 berkas keluaran bernama sama** per bahasa (`808` berkas keluaran terverifikasi di `out/<lang>/`, ditambah `96` SVG bahasa + `8` pasang grafik gabungan PNG+SVG di `outputs/figures/day008h_*.{png,svg}`).

Hasil utama berdasarkan `outputs/day008h/results.json`:
- **Paritas Bit-Identik (`c1`)**: `8/8` bahasa identik **0-byte** pada `s01_selftest.txt`, `s05_queue_script.csv`, `s06_bag_script.csv`, dan `s40_parity.csv` (`3` hari × `8 640` slot = `25 920` slot).
- **Validasi DTMC Eksak (`c2`–`c4`)**: Selisih relatif `s07`–`s09` antar-8 bahasa `0.00e+00` (`< 1e-10`), selisih bentuk tertutup `Geo/Geo/1` vs `Geo/Geo/1/60` sebesar `1.74e-11` (`< 1e-9`), validasi stasiun tunggal `9/9` dalam 3 SE di C (`>= 8/9` di semua 8 bahasa), dan validasi seri 2-stasiun `4/4` dalam 3 SE di semua 8 bahasa.
- **Kalibrasi & Invarian (`c5`, `c8`)**: Pada simulasi utama C (`365` hari), stasiun tersibuk adalah `SECURITY` dengan utilisasi `0.6284` (`∈ [0.55, 0.92]`) dan rasio kehilangan penumpang `0.1133` (`< 0.15`), dengan `0` pelanggaran konservasi aliran maupun identitas Little (`little_diff = 0`).
- **Kejujuran Ilmiah & Anti-Hampa (`c9`, `c11`, `c12`, `c13`, `c15`, `c20`)**: `vacuous_flags = []`. Empat kriteria (`c9`, `c11`, `c12`, `c13`) dilaporkan jujur sebagai `belum tercapai` dan satu (`c15`) sebagai `sebagian` tanpa penyetelan ambang pasca-eksekusi.

---

## 2. Pertanyaan, Tujuan & Kriteria Keberhasilan A Priori

1. Dapatkah simulator jaringan antrean diskret multi-stasiun (`CHECKIN -> SECURITY -> IMIGRASI -> GATE -> CUSTOMS` + `BAGGAGE_BELT`) diimplementasikan dengan **paritas bit-identik lintas 8 bahasa** melalui inti `int64` dan RNG `SplitMix64 -> xoshiro256++`?
2. Apakah hasil simulasi stokastik cocok secara kuantitatif terhadap **solusi eksak DTMC** (`Geo/Geo/c/K` tunggal dan seri `Geo/Geo/1/K1 -> Geo/Geo/1/K2`) dalam batas galat statistik (`<= 3 SE`)?
3. Apakah mekanisme prioritas `BISNIS`, gangguan/perbaikan server, roster per jam, kesabaran (*abandonment*), dan penutupan gerbang menghasilkan efek terukur (`>= 3 SE`, non-hampa) serta lolos kontrol negatif dan uji metamorfik?

---

## 3. Arsitektur Model & Rancangan 48 Tahap

Setiap bahasa dipecah menjadi **11 modul sejajar**: `rng`, `config`, `dtmc`, `station`, `passenger`, `network`, `engine`, `stats`, `analysis`, `svg`, dan `main` (`projects/c/day008h_airport_queue_octalingual/<lang>/`), serta notebook mandiri `days/day008h_sim_airport-queue-octalingual.ipynb`.

### Fase Slot Tetap (1 Slot = 10 Detik, 1 Hari = 8 640 Slot)
Dalam setiap slot `t`, keadaan diperbarui secara deterministik dalam **8 fase berurutan**:
1. **Kedatangan**: Profil per jam (`arr_prof[hr]`) dikalikan `arr_scale_bp / 10000`; menentukan kelas (`EKONOMI=0`, `BISNIS=1`, `TRANSIT=2`), tujuan (`DOM=0`, `INTL=1`), jumlah bagasi (`0..2`), dan jadwal keberangkatan `dep_slot = t + rand_int(110, 220)`. Penumpang `TRANSIT` langsung menuju `SECURITY`.
2. **Penyelesaian layanan**: Iterasi stasiun `0..5`, server `0..c_max-1`; uji Bernoulli dengan `mu_bp` (atau `security_extra_mu_bp` jika terkena pemeriksaan tambahan).
3. **Rutinisasi & transisi**: `CHECKIN -> SECURITY -> (IMIGRASI jika INTL) -> GATE -> selesai` (atau `CUSTOMS` bila diaktifkan). Di `GATE`, jika `t >= dep_slot - gate_close_lead`, penumpang mengalami `MISSED_FLT` dengan peluang `rebook_prob_bp`. Bagasi menuju `BAGGAGE_BELT` dengan peluang hilang `bag_miss_bp`.
4. **Penugasan server**: Server `IDLE` & aktif mengambil dari kepala antrean; di `CHECKIN` dan `SECURITY`, kelas `BISNIS` mendapat prioritas non-preemptive (`pop_next`).
5. **Abandonment**: Penumpang yang menunggu `>= patience_slots` di stasiun penumpang keluar sebagai `ABANDON`.
6. **Gangguan/perbaikan**: Server aktif dapat rusak (`SRV_DOWN`) dengan peluang `break_prob_bp` (pelanggan dikembalikan ke kepala antrean `push_head`); server `DOWN` pulih dengan peluang `repair_prob_bp`.
7. **Pergantian roster**: Setiap `360` slot (`1` jam), jumlah server aktif disesuaikan dengan tabel `roster[hr][st]`.
8. **Pencatatan & invarian**: Mencatat `sum_l`, `busy_slots`, `hr_util_sum`, dan memeriksa invarian kapasitas/konservasi.

---

## 4. Spesifikasi Eksperimen & Beban Kerja

| Bahasa | Toolchain | LOC | Beban Utama (`S14`) | Beban Sapuan (`S29`–`S33`) | Berkas/Bahasa |
|---|---|---:|---:|---:|---:|
| `c` | `Apple clang 16.0.0 (-std=c11)` | `1587` | `365` hari | `3` hari | `101` |
| `cpp` | `Apple clang 16.0.0 (-std=c++20)` | `1229` | `365` hari | `3` hari | `101` |
| `rust` | `rustc 1.93.1` | `1066` | `365` hari | `3` hari | `101` |
| `go` | `go1.27.4` | `2032` | `365` hari | `3` hari | `101` |
| `ts` | `node v24.21.0 (--experimental-strip-types)` | `1699` | `60` hari | `3` hari | `101` |
| `py` | `python 3.11.18` | `1535` | `7` hari | `1` hari | `101` |
| `ruby` | `ruby 4.0.9` | `988` | `7` hari | `1` hari | `101` |
| `ipynb` | `python 3.11.18 (ipynb)` | `915` | `3` hari | `1` hari | `101` |

Paritas (`S39`/`S40`) dijalankan seragam **`3` hari (`25 920` slot)** di semua 8 bahasa.

---

## 5. Hasil Kuantitatif & Perbandingan Lintas Bahasa

### 5.1 Paritas Bit-Identik & Kecepatan (`S40`, `S42`)
- `s01_selftest.txt`: `ok=1`, `h_rng=3df66645d51ba189`, `inv_ok=1`, `tot_arr=490`, `tot_done=490` — identik 0-byte di 8/8 bahasa.
- `s40_parity.csv`: identik 0-byte di 8/8 bahasa pada `3` hari paritas.
- **Kecepatan (`slot-pax/detik` median $\pm$ MAD)**:
  - `cpp`: `61.27M ± 0.98M`
  - `c`: `58.49M ± 0.27M`
  - `rust`: `57.71M ± 0.59M`
  - `go`: `50.61M ± 1.48M`
  - `ts`: `3.42M ± 17.7k` (`BigInt` pada RNG/hash + `number` untuk penghitung)
  - `ipynb`: `566.9k ± 2.4k`
  - `py`: `542.1k ± 4.9k`
  - `ruby`: `331.2k ± 1.4k`

### 5.2 Validasi DTMC, Invarian & Sapuan (`S07`–`S11`, `S15`, `S29`–`S36`)
- `S09`: `rel_err_L = 1.74e-11` antara rumus tertutup `Geo/Geo/1` dan solusi DTMC terpotong `K=60`.
- `S10` (50 replika, C): `9/9` metrik (`L`, `U`, `p_block` untuk `c ∈ {1,2,3}`) berada dalam `|z| <= 3.0`.
- `S11` (50 replika, C): `4/4` metrik (`L1, L2, U1, U2`) berada dalam `|z| <= 3.0` di seluruh 8 bahasa.
- `S15` & `S41`: `violations = 0`, `little_diff = 0`, `max_diff = 0` di 8/8 bahasa.
- `S29`–`S32`: `z_srv = 35.5`, `z_biz = 7.0`, `z_econ = 0.8`, `z_brk = 3.0`, `z_blt = 7.4`.
- `S36` (CRN): `var_ratio = 1.2264` (`< 2.0` karena desinkronisasi aliran acak akibat *abandonment* dan pengalihan rute).

---

## 6. Visualisasi Gabungan

Delapan grafik gabungan disimpan dalam format PNG dan SVG di `outputs/figures/`:
- `day008h_dtmc_validation.{png,svg}` — Z-score validasi DTMC `S10`
- `day008h_hourly_utilization.{png,svg}` — Utilisasi `SECURITY` per jam lintas beban (`C`, `TS`, `Py`, `ipynb`)
- `day008h_station_queues.{png,svg}` — Rata-rata waktu tunggu per stasiun (`C`, 365 hari)
- `day008h_sweeps_effect.{png,svg}` — Ukuran efek (`z-score`) sapuan arah `S29`–`S32`
- `day008h_ensemble_ci.{png,svg}` — Pita ensambel 30 replika `S34` (`mean ± 1.96 SE`)
- `day008h_tornado_sensitivity.{png,svg}` — Analisis sensitivitas tornado `S45` (`±20%`)
- `day008h_parity_matrix.{png,svg}` — Matriks paritas bit-identik `8 bahasa × 4 berkas`
- `day008h_benchmark_loc.{png,svg}` — Perbandingan *throughput* (`slot-pax/s`, skala log)

---

## 7. Evaluasi Kriteria A Priori

| Kriteria | Status | Nilai Terukur (`outputs/day008h/results.json`) | Ambang A Priori |
|---|---|---|---|
| `c1_bit_parity` | `tercapai` | `8/8 bahasa identik 0-byte pada s01,s05,s06,s40` | `8/8 identik` |
| `c2_dtmc_exact` | `tercapai` | `max_rel_diff=0.00e+00, s09_rel_err=1.74e-11` | `< 1e-10 & < 1e-9` |
| `c3_single_val` | `tercapai` | `C: 9/9 dalam 3 SE (semua 8 bahasa >= 8/9)` | `>= 90% dalam 3 SE` |
| `c4_tandem_val` | `tercapai` | `C: 4/4 dalam 3 SE (8/8 bahasa 4/4)` | `>= 90% dalam 3 SE` |
| `c5_invariants` | `tercapai` | `0 pelanggaran & little_diff=0 pada 8/8 bahasa` | `0 pelanggaran` |
| `c6_replay_det` | `tercapai` | `8/8 bahasa match=1` | `match=1` |
| `c7_replay_indep` | `tercapai` | `8/8 bahasa max_diff=0` | `max_diff=0` |
| `c8_calibration` | `tercapai` | `C max_util=0.6284, loss=0.1133` | `util in [0.55,0.92], loss < 0.15` |
| `c9_util_little` | `belum tercapai` | `max_u_diff=6.72e-04, max_little_rel=1.29e-16` | `< 1e-4` |
| `c10_neg_controls` | `tercapai` | `8/8 bahasa 3/3 kontrol negatif ok=1` | `semua ok=1` |
| `c11_metamorphic` | `belum tercapai` | `S26 match=1, S27 max_diff=0.025043, S28 min_mono=1.00` | `S26=1, S27<0.02, S28>=0.95` |
| `c12_sweeps_dir` | `belum tercapai` | `z_srv=35.5, z_biz=7.0, z_econ=0.8, z_brk=3.0, z_blt=7.4` | `semua >= 3 SE & non-hampa` |
| `c13_crn_ratio` | `belum tercapai` | `C var_ratio=1.2264` | `>= 2.0x` |
| `c14_ci_overlap` | `tercapai` | `8/8 bahasa CI tumpang tindih` | `overlap` |
| `c15_branches` | `sebagian` | `C s_main zero_branches=['BELT_BLOCK'] (BELT_BLOCK terpicu di s32)` | `15/15 cabang >= 1` |
| `c16_100_files` | `tercapai` | `101 berkas/bahasa identik di 8/8 bahasa` | `>= 100 berkas` |
| `c17_all_ran` | `tercapai` | `8/8 bahasa selesai` | `8/8 selesai` |
| `c18_selftest` | `tercapai` | `8/8 bahasa ok=1` | `8/8 ok=1` |
| `c19_svgs` | `tercapai` | `8/8 bahasa 12 SVG valid XML` | `>= 12 SVG/bahasa` |
| `c20_anti_vacuous` | `tercapai` | `vacuous_flags=0` | `0 flags` |
| `c21_compress` | `tercapai` | `8/8 bahasa roundtrip_ok=1` | `roundtrip_ok=1` |
| `c22_bench_loc` | `tercapai` | `C 58489802+-266295 slot-pax/s (1587 LOC)` | `dilaporkan dengan MAD & LOC` |

---

## 8. FAKTA vs HIPOTESIS & Deviasi

### FAKTA
- Seluruh 8 bahasa menyelesaikan 48 tahap tanpa galat dan menghasilkan paritas 0-byte pada `s01`, `s05`, `s06`, `s40`.
- `vacuous_flags` berjumlah `0` setelah pemeriksaan kolom konstan di `compare.py`.
- `c9` gagal tipis pada hukum utilisasi (`6.72e-04 > 1e-4`), sementara identitas Little `s19` terpenuhi hingga `1.29e-16`.
- `c11` gagal tipis pada `S27` (`max_diff = 0.025043` vs `< 0.02`) di beban singkat `1` hari.
- `c12` memiliki `z_econ = 0.8` karena hanya sebagian kecil penumpang `BISNIS` yang mendahului `EKONOMI` (`z_biz = 7.0`).
- `c13` menghasilkan `var_ratio = 1.2264` (`< 2.0`).

### HIPOTESIS
- Selisih `6.72e-04` pada `U - U_law` berasal dari agen yang masih berada dalam sistem pada akhir simulasi (`dep_cnt < arr_cnt`) dan pengembalian pelanggan ke antrean saat gangguan server (`requeue_head`).
- Rendahnya rasio reduksi varians CRN (`1.2264`) disebabkan oleh *desynchronization* aliran RNG bersama ketika perubahan jumlah server mengubah urutan kedatangan/penyelesaian dan kejadian *abandonment*.

### Deviasi
1. **`tsc` tidak tersedia**: TypeScript dijalankan langsung dengan `node v24.21.0 --experimental-strip-types` tanpa `tsc --strict`.
2. **Kriteria gagal dilaporkan apa adanya**: `c9`, `c11`, `c12`, `c13` berstatus `belum tercapai` dan `c15` berstatus `sebagian` tanpa mengubah ambang atau menyetel parameter pasca-eksekusi.

---

## 9. Landasan Teori (Terhubung ke Implementasi)

1. **Antrean Geo/Geo/c/K Diskret** — Model waktu terslot dengan kedatangan Bernoulli $\lambda$ dan keberangkatan binomial $\text{Bin}(\min(n,c), \mu)$ diimplementasikan di `dtmc.{c,cpp,rs,go,rb,py,ts}` (`solve_dtmc_single`).
2. **Rantai Markov Diskret (DTMC) & Matriks Transisi** — Ruang keadaan hingga $\{0,\dots,K\}$ dengan transisi fase awal-slot $P_{ij}$ membangun kebenaran dasar eksak `s07_dtmc_single.csv`.
3. **Distribusi Stasioner & Iterasi Gauss–Seidel** — Pemecahan $\pi P = \pi$ secara iteratif di tempat hingga $\max_j |\pi^{\text{baru}}_j - \pi_j| < 10^{-14}$ memastikan konvergensi numerik lintas 8 bahasa hingga `< 1e-10`.
4. **Rasio Geometrik & Rumus Tertutup Geo/Geo/1** — Rasio $\rho_g = \frac{\lambda(1-\mu)}{\mu(1-\lambda)}$ menghasilkan bentuk tertutup $L = \frac{\lambda(1-\lambda)}{\mu - \lambda}$ pada `solve_geo_closed_form` (`s09_geo_closed_form.csv`).
5. **Properti BASTA (Bernoulli Arrivals See Time Averages)** — Peluang pemblokiran pada `Geo/Geo/c/K` sama dengan probabilitas stasioner keadaan penuh $\pi_K$, divalidasi oleh `s10_single_validation.csv`.
6. **Jaringan Seri Terbatas & Blocking** — Sistem 2-stasiun `Geo/Geo/1/K1 -> Geo/Geo/1/K2` (`49` keadaan) dengan *blocking* saat $n_2 - d_2 = K_2$ dimodelkan di `solve_dtmc_tandem` (`s08`, `s11`).
7. **Batas Validitas Jaringan Jackson** — Jaringan bandara nyata melanggar asumsi Jackson (ada kapasitas terbatas, prioritas non-preemptive, gangguan, dan *abandonment*), sehingga validasi eksak dibatasi pada subsistem `S07`–`S11`.
8. **Hukum Little & Identitas Jalur Sampel** — Pada horizon hingga, $\sum_t L(t) \equiv \sum_i W_i$ berlaku sebagai identitas aljabar eksak jika agen tersisa dihitung hingga akhir simulasi (`s15_invariants.csv`).
9. **Hukum Utilisasi** — $U = \lambda_{\text{eff}} E[S] / c$ diuji di `s16_utilization.csv`; penyimpangan kecil muncul akibat agen tersisa di akhir horizon dan *requeue* akibat gangguan.
10. **Disiplin Prioritas Non-Preemptive** — Penumpang `CLS_BIZ` diprioritaskan di antrean `CHECKIN` dan `SECURITY` tanpa memutus layanan aktif melalui `Station::pop_next`.
11. **Abandonment & Kesabaran Terbatas** — Pelanggan yang menunggu $\ge \text{patience\_slots}$ keluar di Fase 5, mencegah ledakan panjang antrean saat beban puncak.
12. **Blocking & Buffer Terbatas** — Sabuk bagasi `STN_BELT` berkapasitas `belt_cap` menolak bagasi berlebih (`BELT_BLOCK`), diuji secara parametrik pada `s32_sweep_belt.csv`.
13. **Proses Gangguan & Perbaikan Server** — Rantai 2-keadaan (`IDLE/BUSY <-> DOWN`) dengan peluang `break_prob_bp` dan `repair_prob_bp` mengembalikan agen aktif ke kepala antrean (`push_head`).
14. **Deteksi Pemanasan MSER-5** — *Marginal Standard Error Rule* blok-5 (`s12_warmup.csv`) memilih panjang pemanasan $d^*$ yang meminimalkan $\text{Var}(\bar{Y}_{n,d})$.
15. **Replikasi Independen vs Batch Means** — Perbandingan lebar selang kepercayaan `95%` antara 20 replika independen dan rata-rata simulasi utama dilakukan pada `s35_ci_methods.csv`.
16. **Common Random Numbers (CRN) & Desinkronisasi** — Pemakaian aliran RNG yang sama antar konfigurasi (`s36_crn.csv`) memberi reduksi varians moderat (`1.2264x`) karena cabang stokastik menggeser sinkronisasi.
17. **Uji Metamorfik & Kontrol Negatif** — Invariansi permutasi server (`s26`), penskalaan waktu (`s27`), monotonisitas server (`s28`), dan 3 kontrol negatif (`s25`) memverifikasi perilaku kualitatif model.
18. **Determinisme Bilangan Bulat & FNV-1a 64-bit** — Penghindaran *floating point* di inti simulasi dan pemakaian `uint64`/`BigInt` pada `SplitMix64`, `xoshiro256++`, dan FNV-1a menjamin paritas 0-byte lintas 8 bahasa.
19. **Replay Independen dari Log Kejadian** — Rekonstruksi deret antrean hanya dari `station_events_*.csv` (`+1`/`-1`) di `S41` memverifikasi konsistensi pencatatan (`max_diff = 0`).
20. **Kriteria Anti-Hampa (`vacuous_flags`)** — Pemeriksaan otomatis di `compare.py` menolak kriteria yang lolos karena kolom konstan atau kondisi tanpa variasi nyata.

---

## 10. Kesimpulan

Implementasi `day008h` berhasil menyatukan **48 tahap simulasi jaringan antrean bandara** dengan **paritas bit-identik 0-byte di 8 bahasa**, validasi eksak terhadap **DTMC tunggal dan seri**, **101 berkas keluaran per bahasa**, serta pelaporan jujur atas empat kriteria yang tidak terpenuhi secara penuh (`c9`, `c11`, `c12`, `c13`) dan satu kriteria parsial (`c15`).

---

## 11. Artefak & Reproduksibilitas

- Proyek utama: `projects/c/day008h_airport_queue_octalingual/` (`SPEC.md`, `config.kv`, `Makefile`, `run_all.sh`, `compare.py`, `c/`, `cpp/`, `rust/`, `go/`, `ruby/`, `py/`, `ts/`, `out/`)
- Notebook: `days/day008h_sim_airport-queue-octalingual.ipynb`
- Ringkasan & grafik: `outputs/day008h/results.json`, `outputs/figures/day008h_*.{png,svg}`
- Reproduksi: `cd projects/c/day008h_airport_queue_octalingual && ./run_all.sh && ../../../.venv/bin/python compare.py`
