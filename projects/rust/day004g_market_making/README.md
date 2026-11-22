# Hari 004g — Market Making Optimal: Guéant–Lehalle–Tapia, HJB, dan Salah Spesifikasi (Rust)

Proyek Rust mandiri (*pure standard library*, `edition = "2021"`, nol *crate* eksternal) yang memodelkan **market making optimal** pada batas inventori $q \in [-Q, Q]$ dalam kerangka **Guéant–Lehalle–Tapia (GLT, 2013)** dengan dinamika harga mid $dS_t = \sigma\,dW_t$, intensitas eksekusi eksponensial $\lambda(\delta) = A e^{-k\delta}$, dan fungsi utilitas CARA $U(W) = -e^{-\gamma W}$.

Dokumentasi penjelasan proyek lengkap tersedia di [day004g_quant_market-making.md](../../../docs/project-explanations/day004g_quant_market-making.md).

## Intuisi Ekonomi & Matematis

Seorang *market maker* memasang kuotasi *bid* $S_t - \delta_t^b$ dan *ask* $S_t + \delta_t^a$ secara kontinu:
- **Spread kecil ($\delta$ kecil)** meningkatkan probabilitas eksekusi $\lambda(\delta) = A e^{-k\delta}$, tetapi memberi laba per transaksi yang tipis dan memperbesar risiko inventori.
- **Risiko inventori ($q_t \neq 0$)** mengekspos kekayaan *mark-to-market* terhadap fluktuasi harga $S_t$. Dengan utilitas CARA, agen yang menghindari risiko menggeser kuotasi (*inventory skew*): saat $q > 0$ (*long*), ia menurunkan $\delta^a$ (mempermudah jual) dan menaikkan $\delta^b$ (memperlambat beli), sehingga inventori terdorong kembali menuju $0$.
- **Transformasi eksponensial GLT**: Melalui $v_q(t) = \exp(k\,\theta(t,q))$ (dengan $\alpha = \frac{1}{2}k\gamma\sigma^2$ dan $\eta = A(1+\gamma/k)^{-(1+k/\gamma)}$), persamaan Hamilton–Jacobi–Bellman (HJB) nonlinier tereduksi menjadi sistem ODE linier berkoefisien konstan $\frac{dv}{dt} = M v$ dengan matriks tridiagonal simetris $M$.

## Rumus Kunci

- **Fungsi nilai**: $u(t, s, q, x) = -\exp\bigl(-\gamma(x + q s + \theta(t, q))\bigr)$, di mana $\theta(t, q) = \frac{1}{k}\ln v_q(t)$.
- **Kuotasi optimal**:
  $$\delta^{b*}(t, q) = \frac{1}{k}\ln\!\left(\frac{v_q(t)}{v_{q+1}(t)}\right) + \frac{1}{\gamma}\ln\!\left(1 + \frac{\gamma}{k}\right), \qquad \delta^{a*}(t, q) = \frac{1}{k}\ln\!\left(\frac{v_q(t)}{v_{q-1}(t)}\right) + \frac{1}{\gamma}\ln\!\left(1 + \frac{\gamma}{k}\right).$$
- **Sistem ODE linear GLT**:
  $$\frac{dv_q}{dt} = \alpha q^2 v_q - \eta\bigl(\mathbf{1}_{q < Q}\,v_{q+1} + \mathbf{1}_{q > -Q}\,v_{q-1}\bigr), \qquad v_q(T) = 1.$$
- **Ekivalen Kepastian (*Certainty Equivalent*)**:
  $$\mathrm{CE} = -\frac{1}{\gamma}\ln\!\bigl(\mathbb{E}[e^{-\gamma W_T}]\bigr) \approx \theta(0, 0).$$

## Struktur Kode

```
projects/rust/day004g_market_making/
├── Cargo.toml
├── meta.json
├── README.md
├── src/
│   ├── lib.rs          # Konstanta skenario, tipe hasil, dan re-export modul
│   ├── rng.rs          # xoshiro256++ + Box-Muller + SplitMix64 (deterministik per indeks path)
│   ├── parallel.rs     # Worker pool std::thread::scope deterministik
│   ├── stats.rs        # Rata-rata, SD, MCSE, korelasi, dan paired z-score
│   ├── linalg.rs       # Eigensolver Jacobi siklik & eksponensial matriks scaling-and-squaring
│   ├── ode.rs          # RK4 + ekstrapolasi Richardson
│   ├── glt.rs          # Solusi eksak v-system, celah spektral, state dasar stasioner, hampiran closed-form
│   ├── hjb_fd.rs       # Solver FD eksplisit mundur langsung pada bentuk nonlinier theta(t,q)
│   ├── kolmogorov.rs   # ODE Forward Kolmogorov P_q(t) untuk E[exec] dan Var(q_T) teoretis
│   ├── strategies.rs   # Strategi GLT, Avellaneda-Stoikov, Symmetric Fixed, Tight Fixed, dan AS-HJB
│   ├── sim.rs          # Simulator Monte Carlo R=200.000 lintasan CRN & estimasi CE
│   ├── misspec.rs      # Studi salah spesifikasi: adverse selection, sigma, intensitas (A, k), drift mu
│   ├── jsonout.rs      # Penulis JSON 17 digit presisi tanpa dependensi eksternal
│   ├── svg.rs          # Generator 7 grafik SVG
│   └── bin/run.rs      # Runner utama yang menghasilkan outputs/day004g/results.json dan SVG
└── tests/
    └── integration_tests.rs
```

## Cara Build, Test, dan Run

Jalankan dari direktori proyek `projects/rust/day004g_market_making`:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run
```

Keluaran utama ditulis ke [../../../outputs/day004g/results.json](../../../outputs/day004g/results.json) dan 7 grafik ke [../../../outputs/figures/](../../../outputs/figures/).

## Ringkasan Hasil Eksperimen (`main_seed = 20261122`, `R = 200.000`, `8` worker, `47,05s`)

Parameter acuan: $T = 1$, $\sigma = 0{,}3$, $\gamma = 0{,}1$, $k = 1{,}5$, $A = 140$, $Q = 20$, $\alpha = 0{,}00675$, $\eta = 49{,}850378$.

| Panel / Besaran | Nilai Terukur | Acuan / Pembanding | Catatan |
| --- | --- | --- | --- |
| $\theta(0, 0)$ eksak (GLT Jacobi) | `66,261588` | `max rel err` vs `expm`: `4,97e-14`; vs `RK4-Rich`: `2,80e-11` | `c1 = tercapai` |
| Konvergensi HJB FD ($\Delta t$) | Error `1,48e-4` → `1,84e-5` | Orde empiris `1,0032` → `1,0008` | `c2 = tercapai` |
| Simetri & monotonisitas kuotasi | `max asym = 9,38e-15` | $\delta^b(0) = \delta^a(0) = 0{,}648923$, spread `1,297846` | `c3 = tercapai` |
| Celah spektral & batas stasioner | $\Lambda_1 - \Lambda_0 = 1{,}277564$ | Beda kuotasi $T=50$ vs *ground state* $= 0$ | — |
| Ekstrapolasi $\mathrm{CE}_{\mathrm{MC}}(\Delta t)$ | `66,3901` → `66,2981` (`Rich = 66,2684`, $p=1{,}017$) | Eksak `66,2616` ($z_{\mathrm{Rich}} = 0{,}23$, $z_{\mathrm{cont}} = 0{,}88$) | `c4 = tercapai` |
| Kolmogorov vs MC kontinu | $\mathbb{E}[\mathrm{exec}]$: `106,0105` vs `106,0240` ($z=0{,}59$) | $\operatorname{Var}(q_T)$: `62,0174` vs `62,1197` ($z=0{,}60$) | `c5 = tercapai` |
| Perbandingan Strategi (`CE`, CRN) | GLT `66,3285 ± 0,0164` | AS `66,2892`, Sym `66,0818`, Tight `23,3540` | `c6 = tercapai` |
| Salah spesifikasi *Adverse Selection* ($\varepsilon = 0{,}1$) | Unadj `57,5183` vs Shift `57,8590` vs AS-HJB `60,7627` | Reduksi $\varepsilon=0$: `66,3191` vs `66,3191` | `c7 = tercapai` |
| Salah spesifikasi $\sigma$ ($\{0{,}5\sigma, \sigma, 1{,}5\sigma\}$) | CE `66,3074` / `66,3148` / `66,2944` | `max_deg = 0,0203 < 0,1`; $\mathrm{SD}(q_T)$: `8,339` / `7,807` / `7,106` | `c8 = tercapai` |
| Anggaran waktu runner | `53,53s` (`8` worker) | Batas `300s` | `c9 = tercapai` |
| Salah spesifikasi intensitas ($A, k \pm 30\%$) | GLT `46,4034..88,7541` | AS `46,3926..88,6100` (GLT unggul konsisten) | — |
| Salah spesifikasi *drift* ($\mu \in [-0{,}2, 0{,}2]$) | CE `66,2553` ($\pm 0{,}2$) vs `66,2920` ($\mu=0$) | Penalti kecil (~`0,037`), simetris | — |

## Catatan Diagnostik (`c4`, `c5`, `c7`)

- **Orde Richardson `c4`**: Ekstrapolasi memakai orde empiris $p = 1{,}017$ (bukan orde 2) dengan propagasi SE CRN (`0,0295`).
- **Simulator Kontinu `c5`**: Selisih $\operatorname{Var}(q_T)$ diskret terhadap Kolmogorov menyusut seorde $O(\Delta t)$ (`2,617` → `1,359` → `0,468` → `0,313`), dan simulator kontinu *Lewis–Shedler thinning* cocok dengan Kolmogorov dalam `< 1 SE`.
- **Reduksi $\varepsilon = 0$ `c7`**: `exact_as_hjb` kini memakai tabel waktu penuh $(t, q)$ sehingga pada $\varepsilon = 0$ identik dengan GLT (`66,3191`, $\operatorname{Var}(q_T) = 60{,}740$).

## Daftar Grafik

- [day004g_quotes-and-spread.svg](../../../outputs/figures/day004g_quotes-and-spread.svg)
- [day004g_hjb-dt-convergence.svg](../../../outputs/figures/day004g_hjb-dt-convergence.svg)
- [day004g_mc-utility-dt-extrapolation.svg](../../../outputs/figures/day004g_mc-utility-dt-extrapolation.svg)
- [day004g_strategy-comparison.svg](../../../outputs/figures/day004g_strategy-comparison.svg)
- [day004g_adverse-selection-misspec.svg](../../../outputs/figures/day004g_adverse-selection-misspec.svg)
- [day004g_sigma-and-drift-misspec.svg](../../../outputs/figures/day004g_sigma-and-drift-misspec.svg)
- [day004g_example-path-s-q-wealth.svg](../../../outputs/figures/day004g_example-path-s-q-wealth.svg)
