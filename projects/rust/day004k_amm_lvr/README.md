# Hari 004k — AMM, Impermanent Loss, LVR, Biaya, dan Likuiditas Terkonsentrasi (Rust)

Proyek Rust mandiri (*pure standard library*, `edition = "2021"`, nol *crate* eksternal) yang mempelajari mikrostruktur **Automated Market Maker (CPMM $x\cdot y = L^2$)** dan **Loss-Versus-Rebalancing (LVR)** pada data sintetis: bentuk tertutup nilai pool & *impermanent loss*, *anchor* martingal GBM ($R = 10^6$), verifikasi *pathwise* Itô, arbitrase per blok diskret (kuadratur Gauss–Legendre vs MC), limit kontinu gerak Brownian terefleksi ($\Delta \to 0$), titik impas aliran *noise trader* Poisson, likuiditas terkonsentrasi ala Uniswap v3, perbandingan 5 strategi LP, serta uji salah spesifikasi (lompatan Merton, volatilitas stokastik CIR, *drift* $\mu \ne 0$, dan ukuran blok).

> **Catatan**: Seluruh data pada proyek ini bersifat **sintetis** (PRNG deterministik `xoshiro256++`) dan tidak memuat klaim empiris tentang pasar atau protokol nyata, serta bukan nasihat investasi.

Dokumentasi penjelasan proyek lengkap tersedia di [day004k_quant_amm-lvr.md](../../../docs/project-explanations/day004k_quant_amm-lvr.md).

## Intuisi Ekonomi & Matematis

- **LVR sebagai biaya konveksitas negatif**: Penyedia likuiditas (LP) CPMM memegang nilai cekung $V(P) = 2L\sqrt{P}$ dengan gamma negatif $V_{PP} = -\frac{1}{2}LP^{-3/2} < 0$. Karena harga pool tertinggal di antara blok dan arbitraseur mengeksekusi terhadap kurva cekung pada harga usang, portofolio replikasi delta tanpa biaya mengungguli pool pada laju sesaat
  $$\mathrm{LVR}_{\mathrm{rate}} = -\frac{1}{2}\sigma^2 P^2 V_{PP} = \frac{\sigma^2}{8}V.$$
- **Pita biaya & limit refleksi**: Dengan biaya input $\varphi > 0$, arbitrase hanya menguntungkan bila $|\ln(P'/p)| > -\ln(1-\varphi) \approx \varphi$. Pada limit $\Delta \to 0$, salah-harga terefleksi di $[-\delta, +\delta]$ dengan laju waktu lokal total $\frac{\sigma^2}{2\delta}$, menghasilkan laju volume arbitrase $\frac{\sigma^2}{8\delta}V$ dan laju pendapatan biaya LP $\varphi\frac{\sigma^2}{8\delta}V \approx \frac{\sigma^2}{8}V$.
- **Likuiditas terkonsentrasi ($[P_a, P_b]$)**: Menyempitkan rentang meningkatkan efisiensi modal $\mathrm{CE}$ untuk $V_0 = 1$ yang sama; karena $V_{PP} = -\frac{1}{2}LP^{-3/2}$ sebanding dengan $L$, laju $\mathrm{LVR}/V_0$ meningkat sebesar faktor $\mathrm{CE}$ selama harga berada di dalam rentang.

## Rumus Kunci

- **CPMM & *Impermanent Loss***: $V(P) = 2L\sqrt{P}$, $\Delta = L/\sqrt{P}$, $\Gamma = -\frac{1}{2}LP^{-3/2}$, $\mathrm{IL}(r) = \frac{2\sqrt{r}}{1+r} - 1$.
- **Momen GBM ($\mu = 0$)**: $\mathbb{E}[V(P_T)] = V_0 e^{-\sigma^2 T/8}$, $\operatorname{Var}[V(P_T)] = V_0^2(1 - e^{-\sigma^2 T/4})$.
- **Ambang blok & laba arbitrase**: Untuk $z > -\ln(1-\varphi)$, $u = \sqrt{(1-\varphi)e^z}$, laba arbitrase $\frac{L\sqrt{p}}{1-\varphi}(u-1)^2$ dan biaya LP $\frac{\varphi}{1-\varphi}L\sqrt{p}(u-1)$.
- **Syarat impas volume/TVL**: $\varphi \cdot \frac{\text{Volume Total}}{V} \ge \frac{\sigma^2}{8}$ (`0,0731/hari` pada $\sigma = 0{,}8, \varphi = 0{,}003$).

## Struktur Kode

```
projects/rust/day004k_amm_lvr/
├── Cargo.toml
├── meta.json
├── README.md
├── src/
│   ├── lib.rs         # Re-export modul
│   ├── rng.rs         # xoshiro256++ + SplitMix64 + Box-Muller + Poisson
│   ├── parallel.rs    # Worker pool std::thread::scope deterministik (blok tetap 64)
│   ├── stats.rs       # Mean, SD, SE, Sharpe, korelasi
│   ├── special.rs     # PDF/CDF normal (~1e-15) & kuadratur Gauss-Legendre
│   ├── cpmm.rs        # Invarian, V, delta, gamma, IL, momen GBM (R=10^6), LVR Ito
│   ├── blocks.rs      # Mekanik arbitrase per blok (biaya input) & integral GL vs MC
│   ├── reflect.rs     # Limit Brownian terefleksi & titik impas noise trader Poisson
│   ├── conc.rs        # Likuiditas terkonsentrasi [Pa, Pb], gamma, waktu keluar GBM
│   ├── strategies.rs  # Perbandingan 5 strategi LP 30 hari (CRN)
│   ├── misspec.rs     # Lompatan Merton, vol stokastik CIR, drift mu!=0
│   ├── jsonout.rs     # Penulis JSON 17 digit
│   ├── svg.rs         # Generator 8 grafik SVG
│   └── bin/run.rs     # Runner utama Bagian 1-4 -> results.json & SVG
└── tests/
    └── integration_tests.rs
```

## Cara Build, Test, dan Run

```bash
cd projects/rust/day004k_amm_lvr
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run
```

Keluaran utama ditulis ke [../../../outputs/day004k/results.json](../../../outputs/day004k/results.json) dan 8 grafik SVG ke [../../../outputs/figures/](../../../outputs/figures/).

## Ringkasan Hasil Eksperimen (`main_seed = 20261122`, `8` worker, `86,58s`)

| Kriteria | Nilai Terukur | Acuan / Ambang | Status |
| --- | --- | --- | --- |
| `c1_gbm_pool_moments_3se` | `max_z = 2,31` ($R=10^6$, 4 horizon) | `<= 3 SE` | `tercapai` |
| `c2_pathwise_lvr_ito_matches` | `0,0065501` vs `0,0065516` (`diff = 1,5e-6`) | `< 2e-4` | `tercapai` |
| `c3_block_gl_matches_mc_3se` | `max_z = 1,65` pada `12` pasangan `(Delta, phi)` | `<= 3 SE` | `tercapai` |
| `c4_reflected_limit_matches_1s_block` | Limit `0,07988` vs `1s` MC `0,07761` | Selisih `< 25%` | `tercapai` |
| `c5_concentrated_lvr_scales_with_ce` | `+-5%: 3,238`, `+-10%: 1,635`, `Full: 0,080` (`max_z = 1,13`) | `<= 3 SE` & monoton | `tercapai` |
| `c6_merton_jumps_increase_lvr` | Lompatan `0,08678` > difusi `0,07973` | `LVR_jump > LVR_diff` | `tercapai` |
| `c7_stoch_vol_lvr_matches_ev_over_8` | `0,07968` vs `0,07968` (`z = 0,07`, `corr = 0,923`) | `<= 3 SE` | `tercapai` |
| `c8_drift_invariance_lvr_per_v` | `LVR/int V dt = 0,08001 (mu=-0,5)` & `0,07999 (+0,5)` | `0,0800 ± 0,005` | `tercapai` |
| `c9_runtime_under_480s` | `86,58s` (`8` *worker*) | `< 480s` | `tercapai` |

## Daftar Grafik

- [day004k_vp-and-il.svg](../../../outputs/figures/day004k_vp-and-il.svg)
- [day004k_ev-vs-mc.svg](../../../outputs/figures/day004k_ev-vs-mc.svg)
- [day004k_lvr-vs-delta.svg](../../../outputs/figures/day004k_lvr-vs-delta.svg)
- [day004k_block-profit-vs-z.svg](../../../outputs/figures/day004k_block-profit-vs-z.svg)
- [day004k_reflect-vs-block.svg](../../../outputs/figures/day004k_reflect-vs-block.svg)
- [day004k_conc-range-lvr.svg](../../../outputs/figures/day004k_conc-range-lvr.svg)
- [day004k_lp-strategies.svg](../../../outputs/figures/day004k_lp-strategies.svg)
- [day004k_misspec-lvr.svg](../../../outputs/figures/day004k_misspec-lvr.svg)
