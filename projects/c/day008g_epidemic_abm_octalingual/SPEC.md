# SPEC — Hari 8g (`day008g`): Simulasi Epidemi Berbasis Agen SEIR Berlapis Oktalingual (39 Tahap)

> **Catatan Duplikasi (Langkah 0)**: Topik epidemi berbasis agen belum ada di `docs/project-explanations/` (hanya `docs/ROADMAP.md` item 034 SIR sederhana yang belum dicentang; proyek ini tidak menyentuh atau mencentang item 034).

## 1. Model Kompartemen & Jaringan Kontak Berlapis
- **Status Penyakit (`0..8`)**:
  - `0: S` (Susceptible), `1: E` (Exposed), `2: Ia` (Infectious Asymptomatic), `3: Is` (Infectious Symptomatic), `4: H` (Hospitalized), `5: ICU` (Intensive Care), `6: R` (Recovered), `7: D` (Dead).
  - Transisi: `S -> E -> (Ia | Is) -> (R | H -> (ICU -> (R | D)) | R | D)` serta waning `R -> S` (peluang `waning_bp` per hari setelah `waning_min_days`).
- **3 Lapisan Kontak**:
  1. **Rumah Tangga (`HH`)**: ukuran `1..6` orang, derajat `size - 1`, kontak penuh setiap hari.
  2. **Sekolah / Tempat Kerja (`SW`)**: anak umur `0–17` masuk sekolah (~25/kelas, `sw_contacts_child = 6`/hari); dewasa `18–59` masuk tempat kerja (~15/unit, `sw_contacts_adult = 5`/hari); lansia `60+` tidak punya kontak kerja/sekolah (`0`).
  3. **Komunitas (`COMM`)**: acak seluruh populasi; `90%` agen standar (`comm_contacts = 3`/hari), `10%` *superspreader* (`comm_contacts_super = 12`/hari) agar distribusi keturunan heterogen (`k < 1`).
- **Rumah Sakit & ICU**:
  - Kapasitas: `hosp_beds = max(1, (N * hosp_beds_per_10k) / 10000)`, `icu_beds = max(1, (N * icu_beds_per_10k) / 10000)`.
  - Jika penuh saat masuk `H` atau `ICU`, pasien **ditolak** (cabang `HOSP_REJECT` / `ICU_REJECT`) dan peluang kematian berlipat `2x` selama fase akut.
- **2 Varian & Intervensi**:
  - **Varian 0** (asal) dan **Varian 1** (masuk pada `variant_day` dengan pengali transmisi `variant_mult_bp / 10000`, default `1.5x`).
  - **Intervensi**: tutup sekolah (`school_close_day`), lockdown (`lockdown_day`, kurangi kontak SW & COMM), masker (`mask_day`, kurangi peluang tular `mask_eff_bp`), tes-isolasi-penelusuran (`test_day`, simtomatik dites dengan `test_prob_bp` * `compliance_bp`, positif -> isolasi, anggota HH dikarantina dengan `trace_prob_bp`), vaksinasi bocor memudar (`vax_day`, dosis harian `vax_daily_bp` * `N / 10000`, prioritas `60+ -> 18-59 -> 0-17`, efikasi infeksi `vax_eff_bp` yang memudar linear hingga `waning`).

## 2. Prinsip Paritas Bilangan Bulat & RNG
- **Inti simulasi murni `int64`**: Seluruh peluang dalam **basis points** (`1/10000`, `BP = 10000`), durasi dalam **hari**, hitungan dalam `int64`. Tanpa `float`/`exp`/`log` pada inti simulasi (`S05`, `S06`, `S09`, `S10`, `S11`, `S33`, `S34`).
- **RNG**: `SplitMix64` menginisialisasi `xoshiro256++`:
  - Aliran per tugas: `seed = splitmix64(master_seed ^ fnv1a_u32pair(etype, eid))` (`0: pop, 1: contact, 2: disease, 3: interv, 4: hosp`).
  - `rand_int(lo, hi) = lo + int64((next_u64() >> 11) % uint64(hi - lo + 1))`.
- **Sinkron & Deterministik**: Dalam satu hari, infeksi baru dikumpulkan pada iterasi penular berurutan `id` naik (`0..N-1`), duplikat target `S` dibuang menurut kemunculan pertama, lalu diterapkan serentak di akhir hari.

## 3. Invarian Eksak (`S12`) & Anchor Analitik (`S15`, `S16`)
- **Invarian Eksak**:
  1. `S + E + Ia + Is + H + ICU + R + D == N` setiap hari.
  2. `D[t] >= D[t-1]` (kematian kumulatif monoton tak-turun).
  3. `hosp_occ <= hosp_cap` dan `icu_occ <= icu_cap` setiap hari.
  4. Jumlah entri di `transmissions_log.csv` == jumlah transisi `S -> E` kumulatif.
- **Anchor Analitik (Matriks Generasi Berikutnya)**:
  - $R_0 = \frac{\beta_{\text{bp}}}{10000} \cdot \left(\bar{c}_{\text{HH}} + \bar{c}_{\text{SW}} + \bar{c}_{\text{COMM}}\right) \cdot D_{\text{inf}} \cdot \left(p_{\text{Is}} + p_{\text{Ia}} \cdot r_{\text{Ia}}\right)$, dengan $D_{\text{inf}}$ rata-rata durasi infeksius (`4.5` hari), $p_{\text{Ia}} = 0.35$, $p_{\text{Is}} = 0.65$, $r_{\text{Ia}} = 0.60$.
  - Ukuran akhir teoretis $z = 1 - \exp(-R_0 z)$ (iterasi titik tetap 80 langkah) dan ambang kekebalan kawanan $p_c = 1 - 1/R_0$.

## 4. Alur 39 Tahap (`s01`–`s39`, `all`) & 77 Berkas Keluaran
Setiap bahasa menulis berkas identik ke `out/<lang>/`:
- `s01`: `s01_selftest.txt`
- `s02`: `s02_config.csv`
- `s03`: `s03_population.csv`, `s03_households.csv`, `s03_schools_work.csv`
- `s04`: `s04_network.csv`
- `s05`: `s05_contact_script.csv` (>= 10 kasus)
- `s06`: `s06_progression_script.csv` (>= 10 kasus)
- `s07`: `s07_interventions.csv`
- `s08`: `s08_calibration.csv` (cari `beta_bp` sehingga $R_0 \in [1.5, 3.5]$ & serangan tanpa intervensi $\in [20\%, 80\%]$)
- `s09`: `s09_baseline_small.csv`
- `s10`: `daily_counts_nointerv.csv`, `transmissions_log_nointerv.csv`
- `s11`: `daily_counts.csv`, `incidence_by_age.csv`, `hospital_daily.csv`, `deaths_by_age.csv`, `tests.csv`, `quarantine.csv`, `vaccination.csv`, `variant_share.csv`, `transmissions_log.csv`, `generation_interval.csv`, `agents_final.csv`, `progression_log.csv`
- `s12`–`s37`: `s12_invariants.csv`, `s13_rt.csv`, `s14_growth.csv`, `s15_r0.csv`, `s16_final_size.csv`, `s17_attack.csv`, `s18_offspring.csv`, `s19_hospital.csv`, `s20_ifr.csv`, `s21_testing.csv`, `s22_vaccine.csv`, `s23_variant.csv`, `s24_waning.csv`, `s25_negative_controls.csv`, `s26_sweep_comp_test.csv`, `s27_sweep_vax.csv`, `s28_sweep_school_lock.csv`, `s29_sweep_variant_mask.csv`, `s30_ensemble.csv`, `s31_branches.csv`, `s32_vacuous_check.csv`, `s33_replay.csv`, `s34_parity.csv`, `s35_bench.csv`, `s36_sizes.csv`, `s37_replay_counts.csv`
- `s38`: 10 SVG tulisan tangan di `out/<lang>/figures/{epi_curves,incidence_age,rt_series,hospital_cap,layer_attack,offspring_dist,sweep_heatmap,vaccine_effect,variant_share,ensemble_band}.svg`
- `s39`: `s39_summary.json`, `files_manifest.csv`

## 5. Beban Kerja
- **Paritas (`S33`/`S34`)**: `N = 2 000`, `60` hari (identik 8 bahasa).
- **Utama (`S10`/`S11`)**: C, C++, Rust, Go = `N = 100 000`, `180` hari; TypeScript = `N = 20 000`, `180` hari; Ruby, Python = `N = 3 000`, `180` hari; Notebook = `N = 1 500`, `180` hari.
- **Sapuan & Ensambel (`S26`–`S30`)**: C, C++, Rust, Go, TS = `N_sw = 1 000`, `90` hari; Ruby, Python, Notebook = `N_sw = 500`, `90` hari.
