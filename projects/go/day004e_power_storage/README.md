# day004e — Harga Listrik Negatif dan Arbitrase Baterai (DP Stokastik, Go)

Modul Go 1.22+ mandiri (*zero external dependencies*) untuk memodelkan pasar listrik sintetis 3-rezim dengan harga negatif (*oversupply*) dan lonjakan kelangkaan (*compound Poisson*), serta mengoptimalkan *dispatch* penyimpanan energi baterai melalui *stochastic dynamic programming* (DP) Bellman 24 jam secara *nonanticipative*.

> **Catatan Ruang Lingkup**: Seluruh parameter harga dalam proyek ini adalah **model sintetis bergaya pasar listrik**, bukan kalibrasi terhadap data pasar nyata. Hasil evaluasi strategi hanya berlaku bagi spesifikasi model, grid harga, aksi diskret baterai, dan fungsi nilai terminal yang didefinisikan di sini.

## 1. Intuisi Pasar Listrik & Penyimpanan Energi

Berbeda dari aset keuangan biasa, listrik hampir tidak dapat disimpan di jaringan tanpa perangkat *storage*, sementara permintaan berfluktuasi secara musiman sepanjang hari dan pasokan terbarukan (*solar/wind*) maupun pembangkit termal lambat (*must-run*) dapat memicu dua ekstrem harga:
1. **Harga negatif (`oversupply`)**: saat kelebihan pasokan melampaui kapasitas serap atau fleksibilitas penurunan pembangkitan, harga spot dapat turun di bawah `0 EUR/MWh`. Baterai justru **dibayar untuk mengisi daya** (*charging*) selama imbalan harga negatif melebihi biaya degradasi.
2. **Lonjakan kelangkaan (`scarcity`)**: saat cadangan pembangkitan menipis atau terjadi gangguan pasokan, harga dapat melonjak ekstrem melalui *compound Poisson jumps*. Menyimpan energi di jam murah memberi nilai opsi untuk menjual saat kelangkaan (*arbitrase antar-waktu*).

Kebijakan *threshold* statis sederhana (`charge` jika $P \le 0$, `discharge` jika $P \ge 200$) melewatkan arbitrase diurnal pada hari normal maupun kesiapan *State of Charge* (SoC) sebelum peralihan rezim. Sebaliknya, **Stochastic DP** memodelkan nilai kelanjutan $V_{h+1}(r', b', s')$ berdasarkan jam, rezim Markov, dan level harga saat ini tanpa mengintip realisasi masa depan.

## 2. Spesifikasi Model Sintetis & Baterai

### Model Harga Listrik 3-Rezim (1 Jam, `EUR/MWh`)
- **Rezim Markov**: `normal` (0), `scarcity` (1), `oversupply` (2) dengan matriks transisi baris-ke-baris:
  - `normal -> [0.985, 0.010, 0.005]`
  - `scarcity -> [0.250, 0.740, 0.010]`
  - `oversupply -> [0.150, 0.050, 0.800]`
  - Distribusi stasioner eksak: $\pi = [0{,}933817,\; 0{,}040798,\; 0{,}025385]$ (residual `1,04e-17`).
- **Mean musiman**: $m(r, h) = \mu_r + 15 \cos\bigl(2\pi((h \bmod 24) - 19)/24\bigr)$ (puncak sore pukul `19:00`).
- **Parameter per rezim**:
  - `normal`: $\mu=45$, $\phi=0{,}82$, $\sigma=18\text{ EUR/MWh}/\sqrt{\text{jam}}$, $\lambda=0{,}002$, $\text{mean\_jump}=0$ ($J = 0$).
  - `scarcity`: $\mu=90$, $\phi=0{,}60$, $\sigma=100$, $\lambda=0{,}08$, $\text{mean\_jump}=120$, $J = +\sum_{k=1}^N \text{Exp}(120)$.
  - `oversupply`: $\mu=-25$, $\phi=0{,}75$, $\sigma=40$, $\lambda=0{,}05$, $\text{mean\_jump}=60$, $J = -\sum_{k=1}^N \text{Exp}(60)$.
- **Transisi harga tak-terpotong (*unbounded*)**: $P_{t+1} = m(r', h+1) + \phi_{r'}(P_t - m(r', h)) + \sigma_{r'}Z + J$, diawali dari distribusi stasioner $\pi$ dengan *burn-in* `60` hari (`1.440` jam) sebelum jendela evaluasi `24` jam.

### Model Baterai & Rekursi Bellman
- **Baterai**: $E = 4\text{ MWh}$, $\text{SoC}_0 = 2\text{ MWh}$, grid SoC `0,1 MWh` (`41` state `0..40`), $\eta_c = \eta_d = 0{,}90$, batas charge `2 MW`, batas discharge `1,8 MW`, biaya degradasi `2 EUR/MWh` throughput sisi-grid.
- **Lima aksi diskret tanpa charge/discharge bersamaan**:
  1. `charge 2 MW`: $\Delta\text{SoC} = +1{,}8\text{ MWh}$, arus kas $-2P_t - 4$.
  2. `charge 1 MW`: $\Delta\text{SoC} = +0{,}9\text{ MWh}$, arus kas $-P_t - 2$.
  3. `idle`: $\Delta\text{SoC} = 0$, arus kas `0`.
  4. `discharge 0.9 MW`: $\Delta\text{SoC} = -1{,}0\text{ MWh}$, arus kas $+0{,}9 P_t - 1{,}8$.
  5. `discharge 1.8 MW`: $\Delta\text{SoC} = -2{,}0\text{ MWh}$, arus kas $+1{,}8 P_t - 3{,}6$.
- **Nilai terminal & definisi profit bersih**: $\text{terminalValue}(\text{soc}) = 45 \cdot \eta_d \cdot \text{soc}$ (`81 EUR` saat $\text{SoC}_0 = 2\text{ MWh}$) sebagai nilai konservatif energi tersisa. Profit arbitrase bersih harian didefinisikan sebagai $\sum_{h=0}^{23} \text{cashflow}_h + \text{terminalValue}(\text{soc}_{24}) - \text{terminalValue}(\text{soc}_0)$, sehingga strategi `no_battery` dan baterai selalu *idle* tepat bernilai `0 EUR`.
- **Equivalent Full Cycles (EFC)**: didefinisikan secara eksplisit sebagai total perubahan energi tersimpan absolut dibagi dua kali kapasitas: $\text{EFC} = \sum_{h=0}^{23}|\Delta\text{SoC}_h| / (2E)$ (`2E = 8 MWh`).
- **Grid harga DP & kernel**: grid kasar `17` bin (`1.224` state) dan grid halus `33` bin (`2.376` state), keduanya menyertakan *tail bin* terbuka `< -500` (representatif `-600`) dan `>= 2000` (representatif `2500`), dengan `10.000` *draw* kondisional per state (`KernelSeed = 2026112215`, terpisah dari `EvalSeed = 2026112225`).

## 3. Ringkasan Hasil Utama (`outputs/day004e/results.json`, `0,70s` pada `8` worker)

### Verifikasi Teknis (`c1..c8 = tercapai`, `all_technical_passed = true`)
- **Momen Jump Compound Poisson (`200k` draw)**: `scarcity` empiris `mean = 9,8541 ± 0,1073` (teori `9,6`), `var = 2420,32 ± 50,45` (teori `2304`); `oversupply` empiris `-2,9944 ± 0,0424` (teori `-3,0`), `var = 355,24 ± 8,86` (teori `360`); `normal` tepat `(0, 0)`. Seluruh kasus dalam `4 MCSE`.
- **Momen Transisi Satu Langkah (`200k` draw)**: pada 3 titik uji (`normal h=18 P=55`, `scarcity h=23 P=210`, `oversupply h=4 P=-85`), mean & varians empiris cocok terhadap rumus total ekspektasi & varians dalam `4 MCSE`.
- **Bellman vs Brute Force (`4` jam, `625` lintasan)** & Residual Bellman: cocok hingga `1,14e-13` (`626,20 EUR` gross / `545,20 EUR` net), `max_bellman_residual = 0` pada kedua grid (`17` & `33` bin).
- **Determinisme worker**: kernel dan lintasan evaluasi identik bit-demi-bit pada `1` vs `8` worker.
- **Realisasi lintasan evaluasi (`N = 2.000` hari = `48.000` jam)**: `5.608` jam harga negatif (`11,68%`, di `1.240/2.000` hari), `2.061` jam `scarcity`, `1.029` jam `oversupply`, `204` jump non-nol (`151` *scarcity*, `53` *oversupply*, `8` jam multi-jump $N>1$), rentang harga `[-354,69, 1041,92] EUR/MWh` (`q01 = -68,73`, `q50 = 45,51`, `q99 = 215,34`).
- **Massa probabilitas *tail bin***: pada grid `17` bin (`12,24M` draw), tail bawah `< -500` sebesar `0,829%` (mean `-517,88`) dan tail atas `>= 2000` sebesar `2,720%` (mean `2058,27`); pada grid `33` bin (`23,76M` draw), tail bawah `0,431%` (`-518,05`) dan tail atas `1,400%` (`2058,18`).
- **Audit tail diperluas (`expanded_tail_38`, domain `[-1000, 5000]` EUR/MWh, CRN pada bin interior `fine_33`)**: massa tail `< -1000` turun menjadi `0,00616%` (`1.685/27,36M`, mean `-1053,15`) dan `>= 5000` menjadi `0,000%` (`0/27,36M`). Kebijakan `polExp` identik pada `91.511/91.512` state interior (`99,9989%`) dan `48.000/48.000` langkah evaluasi (`100%`), sehingga selisih profit dan regret terhadap `DP-33` tepat `0,00 EUR/hari` (selisih `DP-33` vs `DP-17` murni berasal dari resolusi interior, bukan representasi tail).

### Perbandingan Ekonomi Antar-Strategi (`N = 2.000` hari *Common Random Numbers*)

| Strategi | Mean ± MCSE (EUR/hari) | Median | SD | VaR 5% | ES 5% | Mean EFC | Throughput (MWh) | Mean Regret vs PF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `stochastic_dp_33bin` | `236,09 ± 4,73` | `192,98` | `211,46` | `27,59` | `+4,87` | `0,932` | `7,45` | `252,30 ± 4,56` |
| `stochastic_dp_17bin` | `231,90 ± 4,75` | `189,54` | `212,48` | `22,55` | `-1,85` | `1,020` | `8,16` | `256,49 ± 4,60` |
| `threshold_policy` | `116,01 ± 4,41` | `77,55` | `197,33` | `0,00` | `0,00` | `0,175` | `1,50` | `372,38 ± 4,40` |
| `perfect_foresight_dp` (batas atas) | `488,39 ± 7,16` | `394,22` | `320,34` | `219,18` | `181,26` | `2,505` | `20,08` | `0,00` |
| `no_battery` (benchmark nol) | `0,00` | `0,00` | `0,00` | `0,00` | `0,00` | `0,000` | `0,00` | `488,39 ± 7,16` |

- **Selisih berpasangan & interval 95% MC (`z = 1,96`)**:
  - `DP-33 - Threshold`: **`+120,08 ± 3,47 EUR/hari`** (`95% CI [113,29, 126,88]`, `z = 34,63`, lebih unggul pada `87,85%` hari).
  - `DP-17 - Threshold`: **`+115,89 ± 3,49 EUR/hari`** (`95% CI [109,06, 122,73]`, `z = 33,24`, lebih unggul pada `85,45%` hari).
  - `DP-33 - DP-17`: **`+4,19 ± 1,26 EUR/hari`** (`95% CI [1,72, 6,66]`, `z = 3,32`, lebih unggul pada `37,90%` hari) dengan *turnover* lebih hemat (`EFC 0,932` vs `1,020`) dan ekor bawah `ES 5%` membaik (`+4,87` vs `-1,85 EUR`).
- **Validitas batas**: `0` pelanggaran SoC atau daya dan `0` pelanggaran dominasi *perfect foresight* di seluruh strategi.

## 4. Keterbatasan
- Model harga sintetis mengasumsikan parameter stasioner yang diketahui oleh DP (tanpa risiko estimasi parameter maupun *bidding* pasar *day-ahead* vs *intraday*).
- DP mengasumsikan rezim Markov saat ini (`normal`/`scarcity`/`oversupply`) dapat diamati setiap jam dan menggunakan diskretisasi harga 17/33 bin dengan langkah satu jam tanpa dampak pasar (*price-taker*).
- Keuntungan di atas `threshold_policy` berlaku hanya dalam model sintetis ini dan **tidak** menyiratkan profitabilitas pasar listrik nyata.

## 5. Cara Build, Run, dan Test

```bash
cd projects/go/day004e_power_storage
gofmt -l .
go vet ./...
go test -race -count=1 -v ./...  # 11 paket, 15 tes utama
go run ./cmd/run               # Runner penuh N=2000, draws=10000 (~0.7s pada 8 worker)
go run ./cmd/run -quick        # Mode cepat pengembangan ke temp dir (N=50, draws=500)
```

## 6. Struktur Kode

- `internal/rng/` — Pembangkit deterministik `xoshiro256++` + `SplitMix64` (`MainSeed`, `KernelSeed`, `EvalSeed`, `VerifySeed`), sampler Normal, Eksponensial, dan Poisson.
- `internal/parallel/` — `DefaultWorkers` dan `RunChunks` deterministik.
- `internal/regime/` — Matriks transisi Markov 3-rezim dan penyelesai linear eksak distribusi stasioner $\pi P = \pi$.
- `internal/jumps/` — Sampler dan momen kondisional teoretis *compound Poisson* per rezim.
- `internal/prices/` — Dinamika harga musiman AR(1)+rezim+jump tak-terpotong, momen satu langkah, dan simulator *burn-in* `60` hari.
- `internal/battery/` — Spesifikasi baterai `4 MWh`, 5 aksi diskret tanpa *rounding*, arus kas degradasi, nilai terminal, dan EFC.
- `internal/dp/` — Grid harga 17 & 33 bin (dengan *tail bin*), pembangun kernel transisi paralel (`10.000` draw/state), solver Bellman 24 jam + cek residual, kebijakan *threshold*, dan DP *perfect-foresight*.
- `internal/stats/`, `internal/jsonout/`, `internal/svg/`, `cmd/run/` — Statistik (`VaRAndES`, `PairedDiffCI`), JSON 17-digit, pembuat 5 grafik SVG, dan CLI runner.

## 7. Grafik & Artefak Terkait

- JSON keluaran: [`outputs/day004e/results.json`](../../../outputs/day004e/results.json)
- Grafik SVG:
  - [`day004e_sample-price-regimes-jumps.svg`](../../../outputs/figures/day004e_sample-price-regimes-jumps.svg)
  - [`day004e_example-day-dispatch-soc.svg`](../../../outputs/figures/day004e_example-day-dispatch-soc.svg)
  - [`day004e_daily-profit-distribution.svg`](../../../outputs/figures/day004e_daily-profit-distribution.svg)
  - [`day004e_dp-value-vs-soc-regime.svg`](../../../outputs/figures/day004e_dp-value-vs-soc-regime.svg)
  - [`day004e_grid-resolution-sensitivity.svg`](../../../outputs/figures/day004e_grid-resolution-sensitivity.svg)
- Penjelasan proyek 11 bagian: [`day004e_quant_power-storage.md`](../../../docs/project-explanations/day004e_quant_power-storage.md)
