# Penjelasan Proyek — Hari 8g (`day008g`): Simulasi Epidemi Berbasis Agen SEIR Berlapis Lintas 8 Bahasa (39 Tahap)

## 1. Identitas Proyek
- **Hari / Sesi**: `8g` (`day008g`), **Tanggal**: `2026-11-26` (WITA)
- **Topik**: Simulasi Epidemi Berbasis Agen (ABM) SEIR Berlapis (Rumah Tangga, Sekolah/Tempat Kerja, Komunitas) dengan Rumah Sakit/ICU, 2 Varian, dan Intervensi Multi-Lapis — **39 Tahap Identik** di **8 Bahasa** (`c`, `cpp`, `rust`, `go`, `ruby`, `py`, `ts`, `ipynb`)
- **Direktori Proyek**: `projects/c/day008g_epidemic_abm_octalingual/`
- **Notebook**: `days/day008g_sim_epidemic-abm-octalingual.ipynb`
- **Artefak Keluaran**: `outputs/day008g/results.json`, `outputs/figures/day008g_*.{png,svg}`, dan `out/<lang>/` (`72` berkas/bahasa)
- **Catatan Duplikasi**: Topik epidemi berbasis agen belum ada di `docs/project-explanations/`; item `034` pada `docs/ROADMAP.md` (SIR sederhana) tidak disentuh dan tidak dicentang.

## 2. Masalah yang Dipecahkan
Model kompartemen ODE klasik ($dS/dt = -\beta S I / N$) mengasumsikan pencampuran homogen (*well-mixed*) dan mengabaikan struktur kontak rumah tangga, sekolah/tempat kerja, heterogenitas *superspreader*, keterlambatan pelaporan kasus, serta kapasitas rumah sakit/ICU yang terbatas. Sebaliknya, implementasi ABM lintas bahasa sering kali gagal mempertahankan determinisme atau menghasilkan kriteria evaluasi yang lolos secara hampa (*vacuous*). Proyek ini memecahkan kedua masalah tersebut sekaligus: (1) membangun ABM SEIR berlapis dengan inti aritmetika bilangan bulat murni (`int64`, basis points `1/10000`) sehingga 8 bahasa menghasilkan jejak kejadian dan keadaan akhir **bit-identik**, dan (2) menegakkan aturan **Anti-Hampa** (`vacuous_flags`, kontrol negatif, ukuran efek $\ge 3\text{ SE}$, dan replay independen `S37`).

## 3. Fungsi & Tujuan
1. Menyediakan **39 tahap identik** (`s01`–`s39`, `all`) di 8 bahasa (`c`, `cpp`, `rust`, `go`, `ruby`, `py`, `ts`, `ipynb`) yang menghasilkan **72 berkas keluaran bernama sama** per bahasa.
2. Membuktikan paritas bit-identik (`0` selisih byte pada `s01_selftest.txt`, `s05_contact_script.csv`, `s06_progression_script.csv`, `s34_parity.csv`) pada beban paritas `N = 2 000`, `60` hari.
3. Mengukur dampak intervensi (tutup sekolah, *lockdown*, masker, tes-isolasi-penelusuran, vaksinasi bocor memudar) dan varian ke-2 terhadap serangan akhir, $R_t$, okupansi RS/ICU, dan IFR menurut umur.
4. Memvalidasi seluruh klaim menggunakan kontrol negatif (`S25`), 4 sapuan parameter berulang (`S26`–`S29`), ensambel 20 replika (`S30`), detektor anti-hampa otomatis (`S32` & `compare.py`), dan rekonstruksi log independen (`S37`).

## 4. Uraian Kode & Struktur Berkas
Setiap implementasi bahasa dipecah menjadi 9 modul sejajar:
- `rng`: `SplitMix64` (konstanta `0x9e3779b97f4a7c15`, `0xbf58476d1ce4e5b9`, `0x94d049bb133111eb`), `xoshiro256++`, `fnv1a_bytes`, `fnv1a_u32pair`, `fnv1a_i64`, dan `rand_int(lo, hi) = lo + ((next_u64() >> 11) % (hi - lo + 1))`.
- `population`: pemuat `config.kv`, struktur `Person` (umur `0..2`, `hh_id`, `sw_id`, `status`, `variant`, `vax_doses`, `compliance`, `quar_until`, `super_flag`, `offspring`), serta sintesis CSR rumah tangga (`1..6` orang) dan sekolah/tempat kerja.
- `contact`: langkah kontak harian 3 lapisan (`HH`, `SW`, `COMM`), penegakan masker (`mask_eff_bp * compliance`), karantina (`quar_until >= day`), dan efikasi vaksin bocor yang memudar setelah 30 hari (`vax_eff_bp`).
- `disease`: progresi `S -> E -> (Ia | Is) -> (R | H -> (ICU -> (R | D)) | R | D)` serta waning `R -> S`, ditambah penerapan sinkron kandidat infeksi harian (`cand_tgt`).
- `interventions`: impor varian 1 (`variant_day = 45`), tes-isolasi-penelusuran kontak rumah tangga (`test_prob_bp`, `trace_prob_bp`), dan vaksinasi prioritas umur (`60+ -> 18–59 -> 0–17`).
- `engine`: `analytical_r0` (matriks generasi berikutnya), `solve_final_size` (iterasi titik tetap $z = 1 - e^{-R_0 z}$), `check_invariants` (`S12`), dan loop harian `sim_run`.
- `analysis`: penulisan keluaran `s01`–`s37`, kontrol negatif (`S25`), 4 sapuan (`S26`–`S29`), ensambel (`S30`), dan replay independen (`S37`).
- `svg`: penulisan 10 SVG tulisan tangan (`S38`) ke `out/<lang>/figures/*.svg`.
- `main`: penggerak CLI (`s01` dan `all`), benchmark `s35_bench.csv` (5 ulang median $\pm$ MAD pada `120 000` orang-hari), `s39_summary.json`, dan `files_manifest.csv`.

## 5. Metode & Implementasi
- **Paritas Lintas Bahasa**: Seluruh inti simulasi memakai `int64` (`uint64` untuk RNG/FNV-1a; `BigInt` di TypeScript; masker `0xffffffffffffffff` di Ruby/Python) tanpa `float`/`exp`/`log`. Lima aliran RNG independen (`0: pop`, `1: contact`, `2: disease`, `3: interv`, `4: hosp`) diinisialisasi dengan `master_seed ^ fnv1a_u32pair(100, i)`.
- **Beban Kerja Nyata**:
  - Paritas (`S33`/`S34`): `N = 2 000`, `60` hari di seluruh 8 bahasa.
  - Utama (`S10`/`S11`): C, C++, Rust, Go = `N = 100 000, 180d`; TypeScript = `N = 20 000, 180d`; Ruby & Python = `N = 3 000, 180d`; Notebook = `N = 1 500, 180d`.
  - Sapuan & Ensambel (`S26`–`S30`): C, C++, Rust, Go, TS = `N_sw = 1 000, 90d`; Ruby, Python, Notebook = `N_sw = 500, 90d`.
- **Kalibrasi Baseline (`config.kv`)**: `beta_bp_init = 400`, `vax_daily_bp = 200` menghasilkan baseline tanpa intervensi tidak jenuh ($R_{0,\text{ana}} = 1.7839$, serangan akhir `0.6820`).

## 6. Target & Kriteria Keberhasilan
Kriteria `c1`–`c18` ditetapkan a priori di prompt dan dievaluasi secara otomatis oleh `compare.py` tanpa penyetelan ambang setelah hasil keluar:
- Paritas bit-identik (`c1`), nol pelanggaran invarian (`c2`), replay deterministik (`c3`), replay log independen (`c4`), baseline tidak jenuh (`c5`), $R_0$ empiris vs analitik (`c6`), ukuran akhir vs $z$ (`c7`), efek intervensi $\ge 3\text{ SE}$ & kontrol negatif (`c8`), arah sapuan $\ge 3\text{ SE}$ (`c9`), kapasitas RS (`c10`), dispersi $k < 1$ (`c11`), cakupan cabang (`c12`), `>= 70` berkas/bahasa (`c13`), `all` selesai (`c14`), `selftest` lulus (`c15`), SVG valid (`c16`), `vacuous_flags` kosong (`c17`), dan kecepatan + LOC (`c18`).

## 7. Status Pencapaian Berbasis Bukti

| Kriteria | Status | Nilai Terukur (Bukti `outputs/day008g/results.json`) | Ambang / Keterangan |
|:---:|:---:|---|---|
| `c1` | `tercapai` | `8/8 bahasa bit-identik (S01,S05,S06,S34)` | `8/8` bahasa `0` selisih byte |
| `c2` | `tercapai` | `0 pelanggaran di 8/8 bahasa` | `violations == 0` (`S12`) |
| `c3` | `tercapai` | `8/8 bahasa replay identik` | `match == 1` (`S33`) |
| `c4` | `tercapai` | `8/8 bahasa max_diff == 0` | `max_diff == 0` (`S37`) |
| `c5` | `tercapai` | `attack_nointerv=0.6820, R0_ana=1.7839` | Serangan $\in [0.20, 0.80]$, $R_0 \in [1.5, 3.5]$ |
| `c6` | `belum tercapai` | `R0_emp=1.2398, R0_ana=1.7839, ratio=0.6950` | Rasio di luar $[0.80, 1.25]$ (temuan) |
| `c7` | `tercapai` | `\|0.6820 - 0.7263\| = 0.0443` | Selisih $\le 0.10$ dari $z_{\text{analytic}}$ |
| `c8` | `tercapai` | `no=0.6820 -> iv=0.0056 (z=453.5 SE), neg_ctrl=3/3 OK` | Efek $\ge 3\text{ SE}$ & kontrol negatif `0` |
| `c9` | `belum tercapai` | `comp z=3.31 SE, vax z=3.73 SE, var z=0.17 SE` | Arah varian lemah (`0.17 SE`) saat intervensi aktif |
| `c10` | `tercapai` | `peak_h_main=7/250, peak_h_no=250/250 (rej_no=396)` | RS tidak penuh saat intervensi, penuh tanpa intervensi |
| `c11` | `belum tercapai` | `k_disp=4.0321 (mean=0.9998, var=1.2477)` | Overdispersi moderat, bukan $k < 1$ (temuan) |
| `c12` | `sebagian` | `19/21 cabang di s_main; HOSP_REJ & ICU_REJ di s_no (396)` | Penolakan RS hanya terjadi tanpa intervensi |
| `c13` | `tercapai` | `72 berkas/bahasa di 8/8 bahasa, mismatch=0` | $\ge 70$ berkas bernama identik |
| `c14` | `tercapai` | `8/8 bahasa selesai tanpa galat` | `all` selesai di seluruh bahasa |
| `c15` | `tercapai` | `8/8 selftest ok=1` | `S01` lulus semua bahasa |
| `c16` | `tercapai` | `80 SVG per-bahasa + 6 SVG gabungan valid XML (galat=0)` | Semua SVG lolos parser XML |
| `c17` | `tercapai` | `vacuous_flags=[]` | Tidak ada kolom konstan/identik/bulat semu |
| `c18` | `tercapai` | `C=211.6M pd/s (1284 LOC), py=3.20M pd/s (1002 LOC)` | Dilaporkan dengan MAD di `results.json` |

## 8. Fakta vs Hipotesis
- **Fakta**:
  1. Seluruh 8 bahasa menghasilkan hash event `ev_hash` dan hash akhir `fin_hash` yang identik byte demi byte pada `S34` (`N = 2 000, 60d`), serta `0` pelanggaran invarian (`S12`) dan `0` selisih replay independen (`S37`).
  2. Pada `N = 100 000` (C/C++/Rust/Go), intervensi menurunkan serangan akhir dari `0.6820` menjadi `0.0056` ($z = 453.5\text{ SE}$) dan memangkas puncak okupansi RS dari `250/250` (`396` penolakan) menjadi `7/250` (`0` penolakan).
  3. Rasio $R_{0,\text{emp}} / R_{0,\text{ana}} = 0.6950$ (`c6`) dan $k_{\text{disp}} = 4.0321$ (`c11`) tidak memenuhi ambang a priori, dan sapuan varian (`s29`) hanya berbeda `0.17 SE` saat masker dan vaksin sudah aktif (`c9`).
- **Hipotesis**:
  1. $R_{0,\text{emp}}$ lebih rendah daripada $R_{0,\text{ana}}$ karena saturasi lokal di dalam rumah tangga berukuran `1..6` serta periode laten diskrit (`2` hari) mengurangi laju pertumbuhan awal dibanding asumsi *well-mixed* NGM.
  2. Varian 1 memiliki efek kecil pada `s29` (`0.17 SE`) karena intervensi baseline (masker hari 15, vaksinasi hari 20) sudah menekan $R_t < 1$ sebelum varian masuk dengan 5 kasus impor.

## 9. Landasan Teori
1. **Model Kompartemen SEIR & Perluasan Klinis**: Struktur kompartemen `S -> E -> (Ia | Is) -> (R | H -> (ICU -> (R | D)) | R | D)` memperluas model Kermack–McKendrick dengan periode laten `E`, jalur asimtomatik `Ia`, serta rawat inap dan kematian bertingkat.
2. **Simulasi Berbasis Agen (ABM) vs Persamaan Diferensial**: Berbeda dari ODE deterministik yang mengasumsikan populasi tak-hingga dan pencampuran homogen, ABM melacak keadaan diskrit tiap individu sehingga efek stokastik kepunahan dini dan klaster keluarga muncul alami.
3. **Jaringan Kontak Berlapis (HH, SW, COMM)**: Pemisahan kontak ke rumah tangga (klaster padat ukuran `1..6`), sekolah/tempat kerja (`SW`), dan komunitas acak (`COMM`) merepresentasikan struktur sosial nyata dan memengaruhi ambang epidemi.
4. **Matriks Generasi Berikutnya & $R_0$ Analitik**: Nilai $R_0$ dihitung dari rata-rata derajat kontak per lapisan dikali peluang tular $\beta$, durasi infeksius $D_{\text{inf}}$, dan bobot infektivitas asimtomatik $(p_{\text{Is}} + p_{\text{Ia}} r_{\text{Ia}})$.
5. **Hubungan Ukuran Akhir Epidemi**: Untuk model homogen, fraksi populasi terinfeksi $z$ memenuhi persamaan titik tetap $z = 1 - \exp(-R_0 z)$; selisih `0.0443` pada `c7` menunjukkan kedekatan skala makro dengan deviasi struktur jaringan.
6. **Ambang Kekebalan Kawanan (*Herd Immunity*)**: Fraksi imun kritis $p_c = 1 - 1/R_0$ (`0.4394` untuk $R_0 = 1.7839$) menentukan titik ketika pertumbuhan bersih berbalik negatif.
7. **Estimasi $R_t$ dan *Generation Interval***: Rasio insidensi tersaring `s13_rt.csv` (jendela 5 hari kini vs 5 hari tertinggal 4 hari) memantau perubahan transmisi seiring intervensi, sementara `generation_interval.csv` mencatat jeda waktu infeksi sumber ke target.
8. **Dispersi Keturunan & *Superspreading***: Estimator momen $k = \mu^2 / (\sigma^2 - \mu)$ pada `s18_offspring.csv` mengukur heterogenitas penularan; nilai $k = 4.0321$ mencerminkan overdispersi moderat dari `10%` agen *superspreader*.
9. **Kapasitas Rumah Sakit & Mortalitas Bergantung Beban**: Keterbatasan tempat tidur biasa dan ICU (`250` dan `60` per `100 000`) menggandakan peluang kematian (`2x`) saat pasien ditolak (`HOSP_REJ`/`ICU_REJ`), memodelkan kolaps layanan kesehatan pada skenario tanpa intervensi.
10. **Tes, Isolasi & Penelusuran Kontak**: Pengujian simtomatik yang diskalakan oleh kepatuhan individu diikuti karantina anggota rumah tangga memutus rantai penularan `HH` dan `SW` secara selektif.
11. **Vaksin Bocor & Kekebalan Memudar (*Waning*)**: Efikasi vaksin dimodelkan secara probabilistik (*leaky*) terhadap infeksi dan keparahan, disertai peluruhan linear setelah hari ke-30 dan transisi `R -> S` pasca `waning_min_days`.
12. **Dinamika Varian Lebih Menular**: Masuknya varian 1 (`variant_mult_bp = 15000`, `1.5x`) pada `variant_day` menguji apakah keunggulan transmisi mampu memicu gelombang kedua di bawah tekanan intervensi.
13. **Kontrol Negatif & Kriteria Anti-Hampa**: Tiga kontrol negatif (`zero_contact`, `zero_seed`, `perfect_vax` pada `S25`) serta pemeriksaan varians nol dan keidentikan kelompok (`vacuous_flags`) menjamin metrik tidak lolos secara artifisial.
14. **Determinisme Bilangan Bulat, FNV-1a & Replay Independen**: Penggunaan `SplitMix64`, `xoshiro256++`, `int64` basis points, dan hash FNV-1a memungkinkan paritas bit-identik 8 bahasa serta verifikasi `S37` yang merekonstruksi seluruh kurva harian dari `progression_log.csv` dengan `max_diff = 0`.

## 10. Kesimpulan & Rencana Berikutnya
- **Kesimpulan**: Seluruh 8 bahasa (`c`, `cpp`, `rust`, `go`, `ruby`, `py`, `ts`, `ipynb`) menjalankan **39 tahap** secara deterministik dengan **paritas bit-identik 0-byte**, `0` pelanggaran invarian, `0` selisih replay independen, dan `vacuous_flags = []`. Dari 18 kriteria a priori, **14 `tercapai`**, **1 `sebagian`** (`c12`), dan **3 `belum tercapai`** (`c6`, `c9` varian, `c11` dispersi $k$) yang dicatat jujur sebagai temuan ilmiah model.
- **Rencana Berikutnya**: Jika ingin menguji dominasi varian (`c9` varian) atau *superspreading* ekstrem (`c11` $k < 1$), dapat ditambahkan skenario sapuan tanpa vaksin dini dan distribusi derajat *heavy-tailed* (mis. Pareto) pada lapisan komunitas.

## 11. Jejak Bukti
- **Kode Sumber**: `projects/c/day008g_epidemic_abm_octalingual/{c,cpp,rust,go,ruby,py,ts}/` dan `days/day008g_sim_epidemic-abm-octalingual.ipynb`.
- **Pembanding & Hasil**: `projects/c/day008g_epidemic_abm_octalingual/compare.py`, `outputs/day008g/results.json`.
- **Grafik Gabungan**: `outputs/figures/day008g_{epi_curves,r0_finalsize,sweeps,hospital_ifr,parity_perf,branches}.{png,svg}`.
