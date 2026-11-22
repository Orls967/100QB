# Hari 004j — Deep Hedging dari Nol: Jaringan Saraf, Risiko Entropik, Biaya Transaksi, dan Pasar Tak Lengkap (Rust)

Proyek Rust mandiri (*pure standard library*, `edition = "2021"`, nol *crate* eksternal) yang mengimplementasikan **Deep Hedging** (Buehler et al., 2019) dari nol untuk *short call* Eropa ($S_0 = 1, K = 1, T = 0{,}25, r = 0, n = 30$): arsitektur MLP 2 lapis tersembunyi lebar 32 dengan *backpropagation through time* (BPTT) manual, optimizer Adam, paralelisme deterministik blok tetap (`1` vs `8` *worker* identik bit-demi-bit), verifikasi di dunia Black–Scholes tanpa biaya, zona tanpa-perdagangan (*no-trade band*) pada biaya proporsional vs asimtotik Whalley–Wilmott dan *anchor* DP eksak ($n=10$), pasar tak lengkap Heston QE & harga/delta Fourier Lewis, serta uji salah spesifikasi dan ablasi.

> **Catatan**: Seluruh data pada proyek ini bersifat **sintetis** (PRNG deterministik `xoshiro256++`) dan tidak memuat klaim empiris tentang pasar nyata maupun nasihat investasi.

Dokumentasi penjelasan proyek lengkap tersedia di [day004j_quant_deep-hedging.md](../../../docs/project-explanations/day004j_quant_deep-hedging.md).

## Intuisi Ekonomi & Matematis

- **Risiko entropik & *indifference price***: Untuk P&L lindung nilai $X = p_0 - (S_T - K)_+ + \sum_t \delta_t \Delta S_t - \sum_t c S_t |\delta_t - \delta_{t-1}|$, ukuran risiko entropik adalah
  $$\rho_\lambda(X) = \frac{1}{\lambda}\ln\mathbb{E}\!\left[e^{-\lambda X}\right].$$
  Karena sifat invarian kas $\rho_\lambda(X) = \rho_\lambda(-Z) - p_0$ (dengan $Z = (S_T - K)_+ - \sum \delta_t\Delta S_t + \text{biaya}$), *indifference price* yang membuat $\rho_\lambda(\mathrm{P\&L}(p_0)) = 0$ adalah $p_0^* = \rho_\lambda(-Z) = p_0 + \rho_\lambda(X)$.
- **Zona tanpa-perdagangan (Whalley–Wilmott, 1997)**: Dengan biaya transaksi proporsional $c > 0$, melakukan *rebalance* terus-menerus menuju delta BS menjadi mahal; strategi optimal menahan posisi $\delta_{t-1}$ di dalam pita $[\delta_{\mathrm{BS}} - H, \delta_{\mathrm{BS}} + H]$ dengan setengah-lebar asimtotik
  $$H_{\mathrm{WW}} = \left(\frac{3 e^{-r\tau} c S \Gamma^2}{2\lambda}\right)^{1/3}.$$
- **Pasar tak lengkap Heston**: Dengan hanya *underlying* $S_t$ sebagai instrumen lindung nilai, risiko varians stokastik $v_t$ tidak dapat dihilangkan sepenuhnya. Menambahkan $v_t$ ke input jaringan memungkinkan NN menangkap komponen lindung nilai *minimum-variance* akibat korelasi $\rho < 0$ ($\delta_{\mathrm{MV}} = \Delta_{\mathrm{Heston}} + \frac{\rho\xi}{S}\frac{\partial C}{\partial v}$).

## Rumus Kunci

- **Gradien risiko entropik**: $\dfrac{\partial \rho_\lambda}{\partial X_i} = -\dfrac{e^{-\lambda X_i}}{\sum_j e^{-\lambda X_j}}$ (dihitung stabil dengan pergeseran *log-sum-exp*).
- **SD lindung nilai diskret BS**: $\mathrm{SD}(\mathrm{P\&L}) \approx \sqrt{\dfrac{T}{2n}\int_0^T \mathbb{E}[(S_t^2 \Gamma_t \sigma^2)^2]\,dt}$ (`0,006307` vs MC `0,006244`, rasio `0,9900`).
- **Bellman DP entropik ($n=10$)**: $V_t(x, \delta_{t-1}) = \min_\delta \bigl\{ c|\delta - \delta_{t-1}|S + \frac{1}{\lambda}\ln\mathbb{E}[e^{-\lambda\delta\Delta S + \lambda V_{t+1}(x', \delta)}] \bigr\}$.
- **Heston Lewis (2001) & *Little Heston Trap***: $C = S_0 - \dfrac{\sqrt{S_0 K}e^{-rT}}{\pi}\int_0^\infty \operatorname{Re}\!\left[\dfrac{e^{iu\ln(F/K)}\phi(u - i/2)}{u^2 + 1/4}\right]du$.

## Struktur Kode

```
projects/rust/day004j_deep_hedging/
├── Cargo.toml
├── meta.json
├── README.md
├── src/
│   ├── lib.rs       # Re-export modul
│   ├── rng.rs       # xoshiro256++ + SplitMix64 + Box-Muller
│   ├── parallel.rs  # Worker pool std::thread::scope deterministik (blok tetap 64)
│   ├── stats.rs     # Mean, SD, SE, kuantil, CVaR 95%, entropic_risk
│   ├── nn.rs        # MLP 2x32 (tanh/softplus, keluaran linear), BPTT & gradient check
│   ├── optim.rs     # Optimizer Adam (beta1=0.9, beta2=0.999, eps=1e-8)
│   ├── bs.rs        # BS call/delta/gamma/IV, Leland, Whalley-Wilmott, SD diskret
│   ├── dp.rs        # Solver DP entropik eksak n=10 pada grid (log S, delta_{t-1})
│   ├── heston.rs    # Simulator QE Andersen (2008) dengan koreksi martingal K0*
│   ├── fourier.rs   # Fungsi karakteristik Little Heston Trap + integral Lewis + delta FD
│   ├── hedge.rs     # FeatureSpec & evaluator strategi P&L paralel (R=100k)
│   ├── train.rs     # Loop pelatihan Adam + BPTT mini-batch paralel deterministik
│   ├── eval.rs      # Paired CRN diff, OLS slope, pengukur pita no-trade, tabel delta Heston
│   ├── jsonout.rs   # Penulis JSON 17 digit
│   ├── svg.rs       # Generator 8 grafik SVG
│   └── bin/run.rs   # Runner utama Bagian 1-5 -> results.json & SVG
└── tests/
    └── integration_tests.rs
```

## Cara Build, Test, dan Run

Jalankan dari direktori proyek `projects/rust/day004j_deep_hedging`:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run
```

Keluaran utama ditulis ke [../../../outputs/day004j/results.json](../../../outputs/day004j/results.json) dan 8 grafik SVG ke [../../../outputs/figures/](../../../outputs/figures/).

## Ringkasan Hasil Eksperimen (`main_seed = 20261122`, `8` worker, `1305,62s`)

| Bagian / Kriteria | Nilai Terukur | Acuan / Ambang | Status |
| --- | --- | --- | --- |
| `c1_grad_check_and_adam_under_1e6` | BPTT `<= 3,04e-9`; Adam `<= 1,78e-15` | `< 1e-6` | `tercapai` |
| `c2_worker_invariance_1_vs_8` | `n_params = 353`, `bit_diffs = 0` | Bit-identik | `tercapai` |
| `c3_bs_nn_delta_l2_and_price` | `p0_diff = 1,69e-6`; `L2_mse = 0,00611`; `SD_nn = 0,00627`; `L2_ent = 0,0466` | `L2 <= 0,03`, `SD <= 0,00629` | `belum tercapai` |
| `c4_discrete_hedging_n_scaling` | `slope_BS = -0,4817`; `slope_NN = -0,4532` | `slope_NN in [-0,60; -0,40]` | `tercapai` |
| `c5_cost_notrade_band_and_ww` | `hw_nn = 0,0129 ± 0,0011`, `ratio = 0,077 ± 0,007`, `exp_c = 0,585` | `exp in [0,23; 0,43]`, `ratio ~ 1` | `belum tercapai` |
| `c6_dp_anchor_n10_agrees` | `DP = 0,04191` vs `NN = 0,04171`, `agree = 14,63%` (visit `20,87%`) | Utilitas & kebijakan `>= 90%` | `sebagian` |
| `c7_nn_beats_bs_under_cost` | `NN 0,001834` < `BS 0,002393` & `Leland 0,002361` (> `WW 0,001645`) | Ungguli semua benchmark | `sebagian` |
| `c8_heston_martingale_and_fourier_3se` | `z(E[S_T]=1) = 0,51`; `z(MC vs Fourier) = 0,31` | `<= 3 SE` | `tercapai` |
| `c9_heston_vt_reduces_residual_risk` | `SD: 0,01167 (no v) -> 0,01139 (v_t)`; `rho10: 0,000682 -> 0,000646` | `v_t` menurunkan risiko | `tercapai` |
| `c10_runtime_under_600s` | `1305,62s` (`8` worker) | `< 480s` | `belum tercapai` |

## Daftar Grafik

- [day004j_learning-curve.svg](../../../outputs/figures/day004j_learning-curve.svg)
- [day004j_nn-vs-bs-delta.svg](../../../outputs/figures/day004j_nn-vs-bs-delta.svg)
- [day004j_pnl-distribution.svg](../../../outputs/figures/day004j_pnl-distribution.svg)
- [day004j_notrade-band-vs-ww.svg](../../../outputs/figures/day004j_notrade-band-vs-ww.svg)
- [day004j_dp-vs-nn-n10.svg](../../../outputs/figures/day004j_dp-vs-nn-n10.svg)
- [day004j_heston-strategies.svg](../../../outputs/figures/day004j_heston-strategies.svg)
- [day004j_ablation-width4-noprev.svg](../../../outputs/figures/day004j_ablation-width4-noprev.svg)
- [day004j_misspecification.svg](../../../outputs/figures/day004j_misspecification.svg)
