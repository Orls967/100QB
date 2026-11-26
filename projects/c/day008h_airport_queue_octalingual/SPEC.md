# SPEC — Hari 8h (`day008h`): Jaringan Antrean Terminal Bandara Berbasis Agen dengan Validasi DTMC Eksak (48 Tahap, 8 Bahasa)

> **Catatan Duplikasi (Langkah 0)**: Topik jaringan antrean terminal bandara belum ada di `docs/project-explanations/` (`day001`–`day008g`); hanya `docs/ROADMAP.md` item `026` (antrean kasir M/M/1, belum dicentang) yang tidak disentuh.

## 1. Waktu Diskret, Inti Bilangan Bulat & RNG
- **Waktu diskret**: `1` slot = `10` detik; `1` hari = `8 640` slot; `1` jam = `360` slot.
- **Inti murni `int64`**: Seluruh peluang dalam **basis points** (`BP = 10000`), waktu dalam slot, hitungan dalam `int64`. Tanpa `float`/`exp`/`log` di inti simulasi (`S05`, `S06`, `S10`, `S11`, `S14`, `S25`–`S36`, `S39`, `S40`).
- **RNG**: `SplitMix64` (`0x9e3779b97f4a7c15`, `0xbf58476d1ce4e5b9`, `0x94d049bb133111eb`) menginisialisasi `xoshiro256++`:
  - Aliran per tugas: `seed = splitmix64(master_seed ^ fnv1a_u32pair(200, stream_id))` (`0..7`: `arr, svc, route, break, bag, val, sweep, meta`).
  - `rand_int(lo, hi) = lo + ((next_u64() >> 11) % (hi - lo + 1))`
  - `bernoulli(p_bp) = ((next_u64() >> 11) % 10000) < p_bp`
  - `master_seed = 20261126`.

## 2. Fase Slot Tetap (Urutan Iterasi Terkunci)
Dalam setiap slot `t` (iterasi stasiun `0..5` naik, server `0..c_max-1` naik, pelanggan `seq` naik):
1. **Kedatangan**: Profil kedatangan per jam (`0..23`) menghasilkan penumpang (`EKONOMI=0, BISNIS=1, TRANSIT=2`, `DOMESTIK=0, INTERNASIONAL=1`, `0..2` bagasi) + penugasan penerbangan (`dep_slot`).
2. **Penyelesaian layanan**: Untuk server `BUSY` di awal slot, uji `bernoulli(mu_eff)` (`extra_mu_bp` jika kena pemeriksaan tambahan). Pelanggan selesai dicatat ke daftar transisi.
3. **Rutinisasi & transisi**:
   - `CHECKIN -> SECURITY -> (IMIGRASI jika INTERNASIONAL) -> GATE -> selesai` (`TRANSIT` melewati `CHECKIN`).
   - `GATE`: jika `t >= dep_slot - gate_close_lead`, penumpang ketinggalan (`MISSED_FLIGHT`) → `REBOOKING` atau hilang.
   - `Bag` melewati `BAGGAGE_BELT` (kapasitas buffer terbatas `K_belt`): jika penuh, terjadi `BELT_BLOCK`; dapat mengalami `BAG_MISSED`.
4. **Penugasan server ← antrean**: Server `IDLE` & aktif mengambil kepala antrean (`BISNIS` diprioritaskan non-preemptive di `CHECKIN` dan `SECURITY`; `FIFO` di stasiun lain).
5. **Kesabaran / abandonment**: Penumpang yang masih menunggu dengan `t - queue_enter_slot >= patience` keluar (`ABANDON`).
6. **Gangguan / perbaikan**: Server aktif `IDLE`/`BUSY` gagal dengan peluang `b_bp` (`BREAKDOWN`); jika `BUSY`, pelanggan dikembalikan ke **kepala** antrean (`REQUEUE_HEAD`). Server `DOWN` pulih dengan peluang `r_bp` (`REPAIR`).
7. **Pergantian staf (roster per jam)**: Jumlah server aktif stasiun mengikuti tabel roster jam `hr`; server berlebih yang masih `BUSY` menyelesaikan layanan saat ini sebelum menjadi non-aktif.
8. **Pencatatan & invarian**: Catat panjang antrean per stasiun, status server, akumulasi slot-pelanggan (`Σ L(t)`) untuk identitas hukum Little, serta cek konservasi aliran & kapasitas.

## 3. Kebenaran Dasar Eksak DTMC (`S07`, `S08`, `S09`)
- **Stasiun tunggal `Geo/Geo/c/K` (`S07`)**: `c ∈ {1, 2, 3}`, `K ∈ {5, 10, 60}` (potongan `∞`), peluang datang `λ` dan selesai `μ` per slot (mengikuti urutan awal-slot Phase 1→2). Iterasi Gauss–Seidel hingga `max |π_new - π| < 1e-14` menghasilkan `P_block`, `L`, `L_q`, `U`, `W = L / λ_eff`, `W_q = L_q / λ_eff`.
- **Tertutup `c=1, K=∞` (`S09`)**: Rasio geometrik $\rho_g = \frac{\lambda(1-\mu)}{\mu(1-\lambda)}$, $\pi_0 = 1 - \lambda/\mu$, $L = \frac{\lambda(1-\lambda)}{\mu - \lambda}$, $W = \frac{1-\lambda}{\mu - \lambda}$. Cocok dengan `S07(c=1, K=60)` hingga `< 1e-9`.
- **Seri 2-stasiun `Geo/Geo/1/K1 -> Geo/Geo/1/K2` (`S08`)**: `K1 = 6, K2 = 6` (`49` keadaan), iterasi Gauss–Seidel eksak atas `(A, D1, D2) ∈ {0,1}^3` dengan `blocking` saat `n2_after_D2 == K2`.

## 4. Alur 48 Tahap (`s01`–`s48`, `all`) & ≥ 100 Berkas Keluaran
Setiap bahasa menulis berkas identik ke `out/<lang>/`:
- `s01`–`s06`: `s01_selftest.txt`, `s02_config.csv`, `s03_flights.csv`, `s03_arrival_profile.csv`, `s04_stations.csv`, `s04_roster.csv`, `s05_queue_script.csv`, `s06_bag_script.csv`
- `s07`–`s13`: `s07_dtmc_single.csv`, `s08_dtmc_tandem.csv`, `s09_geo_closed_form.csv`, `s10_single_validation.csv`, `s11_tandem_validation.csv`, `s12_warmup.csv`, `s13_calibration.csv`
- `s14`: `passengers_log.csv`, `bags_log.csv`, `station_events_CHECKIN.csv`, `station_events_SECURITY.csv`, `station_events_IMIGRASI.csv`, `station_events_GATE.csv`, `station_events_BELT.csv`, `breakdowns.csv`, `missed_flights.csv`, `abandonments.csv`, `rebookings.csv`, `queue_series.csv`, `server_status_series.csv`, `agents_final.csv` (+ `daily_summary.csv`, `hourly_summary.csv`, `station_0..5.csv`, `class_0..2.csv`, `wait_hist.csv`, `sojourn_hist.csv`, `bag_claims.csv`, `flight_summary.csv`, `gate_events.csv`, `priority_audit.csv`, `loss_breakdown.csv`, `customs_sweep_sample.csv`)
- `s15`–`s28`: `s15_invariants.csv`, `s16_utilization.csv`, `s17_waits.csv`, `s18_sojourn.csv`, `s19_little.csv`, `s20_missed.csv`, `s21_baggage.csv`, `s22_breakdown_effect.csv`, `s23_roster_effect.csv`, `s24_bottleneck.csv`, `s25_negative_controls.csv`, `s26_meta_permute.csv`, `s27_meta_timescale.csv`, `s28_meta_monotone.csv`
- `s29`–`s48`: `s29_sweep_servers.csv`, `s30_sweep_priority.csv`, `s31_sweep_breakdown.csv`, `s32_sweep_belt.csv`, `s33_sweep_gate.csv`, `s34_ensemble.csv`, `s35_ci_methods.csv`, `s36_crn.csv`, `s37_branches.csv`, `s38_vacuous_check.csv`, `s39_replay.csv`, `s40_parity.csv`, `s41_replay_queue.csv`, `s42_bench.csv`, `s43_sizes.csv`, `s44_compress.csv`, `s45_tornado.csv`, `s47_report_tables.csv`, `s48_summary.json`, `files_manifest.csv`
- `s46`: `12` SVG tulisan tangan di `out/<lang>/figures/{queue_stations,util_hourly,wait_hist,wait_quantiles,dtmc_validation,server_heatmap,priority_effect,breakdown_effect,bottleneck,ensemble_band,crn_ratio,tornado}.svg`

## 5. Beban Kerja & Rentang Kalibrasi A Priori
- **Kalibrasi tak jenuh**: utilisasi stasiun terpadat $\in [0.55, 0.92]$, probabilitas menunggu $> 0$ di semua stasiun aktif, kehilangan penumpang $< 15\%$.
- **Paritas (`S39`/`S40`)**: `3` hari × `8 640` slot (`25 920` slot) di semua 8 bahasa.
- **Utama (`S14`)**: C, C++, Rust, Go = `365` hari; TypeScript = `60` hari; Ruby, Python = `7` hari; Notebook = `3` hari. Log rinci slot (`S14`/`S41`) dicatat untuk `min(main_days, 3)` hari pertama agar ukuran berkas terkendali.
