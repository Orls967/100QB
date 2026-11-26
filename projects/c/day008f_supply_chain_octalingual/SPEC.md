# SPEC — Hari 8f (`day008f`): Simulasi Jaringan Rantai Pasok Oktalingual (40 Tahap)

> **Catatan Duplikasi (Langkah 0)**: Pemeriksaan duplikasi menemukan bahwa topik bursa/order-book/market-making sudah ada pada `projects/cpp/day002h_hawkes_orderbook/` dan `projects/rust/day004g_market_making/`. Oleh karena itu, sesi `day008f` menggunakan topik cadangan resmi: **jaringan rantai pasok 3 pabrik / 4 gudang / 10 toko dengan kebijakan `(s,S)`, peramalan, pengiriman, keterlambatan, dan biaya** dalam alur **40 tahap** di **8 bahasa**.

## 1. Diagram Entitas & Aliran
```
+---------------+     produksi (SKU 0,1,2)     +------------------+
| RawSupplier   | ---------------------------> |  3 Pabrik F0..F2 |
| (Akun 21)     | <---- biaya produksi ------- |  PlantStateMach  |
+---------------+                              +--------+---------+
                                                        | kirim via 3 Kurir (C0..C2)
                                                        v
+---------------+   transship & kit D=A+B+C    +------------------+
| FeePool (22)  | <---- biaya kit / penalti -- | 4 Gudang W0..W3  |
+---------------+                              | W3: Hub Kit + TS |
                                               +--------+---------+
                                                        | alokasi lelang / kirim
                                                        v
+---------------+   permintaan stokastik (s,S) +------------------+
| ExtCustomer   | <---- penjualan ritel ------ | 10 Toko S0..S9   |
| (Akun 20)     | ----- pembayaran ritel ----> | EWMA Forecaster  |
+---------------+                              +------------------+
        ^                                               |
        +--- Clearinghouse & CreditRiskManager ---------+
             (EOD settlement, credit call, likuidasi)
```

## 2. Prinsip Paritas Bilangan Bulat & RNG
- **Inti simulasi murni `int64`**: Tanpa `float`/`double` dan tanpa `exp`/`log`/`sin` pada inti simulasi (`S04`–`S10`, `S12`–`S17`, `S29`–`S34`, `S38`). Seluruh nilai uang dalam **sen** (`int64`), seluruh nilai berada di `[-2^48, 2^48]` sehingga aman bagi `Number` 53-bit di TypeScript/JS.
- **RNG**: `SplitMix64` menginisialisasi 4 state `uint64` untuk `xoshiro256++`. Setiap entitas memiliki aliran sendiri:
  `seed_stream = splitmix64(master_seed ^ fnv1a_u32pair(etype, eid))`
  di mana `fnv1a_u32pair(t, id)` melakukan FNV-1a 64-bit atas 8 byte little-endian `(uint32(t), uint32(id))`.
  `rand_int(lo, hi) = lo + int64((next_u64() >> 11) % uint64(hi - lo + 1))`.
- **Tie-break & Urutan**: Seluruh kejadian dijadwalkan dalam **binary min-heap** dengan kunci `(time, seq)`. Seluruh iterasi entitas berurutan indeks `0..N-1`.
- **Statistik (`S18`–`S28`, `S35`)**: Boleh memakai `float64` dengan penjumlahan sekuensial indeks naik.

## 3. Entitas, Kebijakan, & Sesi Harian
- **4 SKU**: `0=A` (harga ritel `10000`, grosir `8000`, produksi `6000`), `1=B` (`10200`, `8200`, `6200`), `2=C` (`9800`, `7800`, `5800`), `3=D` (Kit `1A+1B+1C`, harga ritel `31500`, grosir `26000`, biaya rakit `1500`).
- **3 Pabrik (`F0..F2`)**: `F_k` memproduksi SKU `k`. Memiliki mesin keadaan `OPEN (0) -> DISRUPTED (1) -> RECOVERY (2) -> OPEN (0)`, kapasitas normal `2400`, lembur `3200`.
- **4 Gudang (`W0..W3`)**: `W0` melayani `S0..S2`, `W1` melayani `S3..S5`, `W2` melayani `S6..S9`, `W3` hub pusat (perakitan Kit `D` + *transshipment* antar gudang). Kapasitas simpan `50000` unit per gudang; penolakan jika penuh.
- **10 Toko (`S0..S9`)**: Peramalan EWMA bilangan bulat (`alpha = 250/1000` milli-unit), kebijakan `(s, S)` dengan safety stock dinamis (`S = 450 * cap_scale`, `s = 0.4 * S`, terkalibrasi non-jenuh), *backorder* atau *lost sale* (jika permintaan melebihi batas atau toko `SUSPENDED`).
- **3 Kurir (`C0..C2`)**: Biaya bertingkat volume harian (`<500`: `200`, `500..1999`: `160`, `>=2000`: `120` sen/unit), *lead time* stokastik `{150, 250, 200} + rand(0,25)` tick, kapasitas harian kurir `12000 * cap_scale` unit/kurir.
- **Definisi Tunggal Fill Rate & Kalibrasi A Priori**:
  - `fill_rate = dem_units_filled / dem_units_tot` (dipakai seragam di `s09_calibration.csv`, `s22_service.csv`, dan `s29_sweep_safety_stock.csv`).
  - Rentang kalibrasi baseline (`s09_calibration.csv`, 1 hari $\times$ 2 000 event): `fill_rate in [0.85, 0.99]`, `bo_event_frac < 0.20`, `carrier_util in [0.30, 0.90]`.
- **Statistik Anti-Hampa (`S21`, `S25`, `S31`)**:
  - `S21` (`s21_bullwhip.csv`): kurtosis sampel, ACF(1), dan *bullwhip ratio* dihitung langsung dari 100 *time bin* permintaan & pesanan per SKU.
  - `S25` (`s25_shock_event_study.csv`): `mean_rel_dev` dan `se_rel_dev` dihitung dari jendela 2-bin pra/pasca-guncangan per kelas (`0, 1, 2`).
  - `S31` (`s31_sweep_transship.csv`): 5 replika per kondisi (`0` vs `1`), mencatat `mean_imbalance` dan `n_transships`.
- **Sesi Harian**: `PRE_DAY -> ALLOCATION_AUCTION -> CONTINUOUS_OPS -> END_OF_DAY_CLEARING -> CLOSED`.
  - **Lelang Alokasi Stok Langka (`ALLOCATION_AUCTION`)**: Ketika total permintaan toko atas `(gudang, sku)` melebihi stok tersedia `avail`, alokasi memaksimalkan total unit terisi; tie-break: (1) minimalkan kekurangan maksimum toko (proporsional/water-filling bulat), (2) prioritas tingkat toko lebih tinggi, (3) `seq` terkecil; sisa permintaan ditandai *backorder* atau ditolak sesuai tipe pesanan (`URGENT_DROP` vs `STANDARD_REPLENISH`).
- **Invarian Eksak** (diperiksa di `S11`):
  1. `sum(cash[0..22]) == INITIAL_TOTAL_CASH`
  2. Untuk `k in {0,1,2}`: `sum(on_hand[k]) + sum(in_transit[k]) + total_sold[k] + total_kits == init_stock[k] + total_produced[k]`
  3. Untuk `k = 3`: `sum(on_hand[3]) + sum(in_transit[3]) + total_sold[3] == init_stock[3] + total_kits`
  4. Seluruh `on_hand`, `in_transit`, `backorder`, `delay >= 0`.

## 4. Alur 40 Tahap (`s01`–`s40`, `all`) & 96 Berkas Keluaran
Setiap bahasa menulis berkas identik ke `out/<lang>/`:
- `s01`: `s01_selftest.txt`
- `s02`: `s02_config.csv`
- `s03`: `s03_entities.csv`, `s03_accounts_init.csv`, `s03_carriers.csv`, `s03_network_edges.csv`
- `s04`: `s04_policy_script.csv` (60 skenario kebijakan `(s,S)`, kapasitas, credit hold, disrupsi, kit leg-risk)
- `s05`: `s05_alloc_script.csv` (10 kasus lelang alokasi termasuk tie-break)
- `s06`: `s06_schedule_head.csv`, `s06_sessions.csv`
- `s07`: `s07_shock_schedule.csv`
- `s08`: `s08_warmup_summary.csv`, `s08_warmup_sku_{A,B,C,D}.csv`
- `s09`: `s09_day1_summary.csv`, `s09_calibration.csv`
- `s10`: `orders_log.csv`, `sales_{A,B,C,D}.csv`, `shipments_{A,B,C,D}.csv`, `inventory_snap_{A,B,C,D}.csv`, `store_daily_{0..9}.csv`, `warehouse_daily_{0..3}.csv`, `factory_daily_{0..2}.csv`, `allocations.csv`, `disruptions.csv`, `credit_calls.csv`, `risk_actions.csv`, `production.csv`, `kit_assembly.csv`, `transshipments.csv`, `entities_final.csv`, `accounts_series.csv`, `shocks_log.csv`, `fees_ledger.csv` (41 berkas)
- `s11`: `s11_invariants.csv`
- `s12`: `s12_clearing.csv`, `s12_credit.csv`
- `s13`–`s28`: `s13_risk.csv`, `s14_factories.csv`, `s15_carrier_fees.csv`, `s16_kits.csv`, `s16_kit_imbalance_series.csv`, `s17_allocations.csv`, `s18_cost_entities.csv`, `s19_cost_types.csv`, `s19_cost_daily.csv`, `s20_demand_{A,B,C,D}.csv`, `s21_bullwhip.csv`, `s22_service.csv`, `s23_leadtime_impact.csv`, `s24_kit_coupling.csv`, `s25_shock_event_study.csv`, `s26_disruptions.csv`, `s27_fees.csv`, `s28_capacity_rejects.csv`
- `s29`–`s31`: `s29_sweep_safety_stock.csv`, `s30_sweep_shock_disrupt.csv`, `s31_sweep_transship.csv`
- `s32`–`s38`: `s32_branches.csv`, `s33_replay.csv`, `s34_parity.csv`, `s35_bench.csv`, `s36_sizes.csv`, `s37_compress.csv`, `s38_ledger_replay.csv`
- `s39`: 10 SVG tulisan tangan di `out/<lang>/figures/{inventory_4sku,fill_rate_stores,cost_by_echelon,leadtime_impact,bullwhip_acf,kit_imbalance,shock_study,alloc_vs_cont,carrier_tiers,sweep_safety}.svg`
- `s40`: `s40_summary.json`, `files_manifest.csv`

## 5. Beban Kerja
- **Beban Paritas (`S33`/`S34`)**: `3` hari $\times$ `8 000` event kontinu (`24 000` event, identik 8 bahasa).
- **Beban Utama (`S10`)**: C/C++/Rust/Go = `5` hari $\times$ `400 000` event; TypeScript = `5` hari $\times$ `80 000`; Ruby/Python = `5` hari $\times$ `8 000`; Notebook = `5` hari $\times$ `4 000`.
- **Beban Sapuan (`S29`–`S31`)**: `2` hari $\times$ `2 000` event per replika (`5` replika per titik; `S_lvl in {80, 220, 350}`).
