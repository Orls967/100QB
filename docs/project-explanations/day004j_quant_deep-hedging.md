# Penjelasan Proyek 004j — Deep Hedging dari Nol: Jaringan Saraf, Risiko Entropik, Biaya Transaksi, dan Pasar Tak Lengkap (Rust)

## 1. Identitas Proyek

- **Kode Sesi**: `4j` (`day004j`) — Hari 4, sesi j (22 November 2026, WITA).
- **Bahasa & Lingkungan**: Rust stabil (`rustc 1.93.1`, *edition* `2021`, *zero external crates*), eksekusi paralel deterministik menggunakan `std::thread::scope` (`8` *worker*) dengan reduksi blok tetap (`BLOCK = 64`) dan `MAIN_SEED = 20261122`.
- **Lokasi Proyek**: [projects/rust/day004j_deep_hedging/](../../projects/rust/day004j_deep_hedging/)
- **Keluaran Utama**: [results.json](../../outputs/day004j/results.json) dan 8 grafik SVG di [outputs/figures/](../../outputs/figures/).
- **Catatan Ruang Lingkup**: Seluruh lintasan harga dan varians pada proyek ini bersifat **sintetis** (dihasilkan oleh PRNG deterministik `xoshiro256++`) dan bukan merupakan klaim empiris mengenai pasar nyata ataupun nasihat investasi.

---

## 2. Permasalahan yang Dipecahkan

Lindung nilai (*hedging*) opsi *short call* Eropa ($S_0 = 1$, $K = 1$, $T = 0{,}25$ tahun, $r = 0$, $n = 30$ langkah diskret) hanya dengan instrumen dasar $S_t$ menghadapi tiga tantangan mendasar di luar dunia Black–Scholes ideal:

1. **Diskretisasi waktu ($n < \infty$)**: Replikasi kontinu eksak mustahil dilakukan secara diskret; strategi delta Black–Scholes $N(d_1)$ menyisakan galat lindung nilai dengan simpangan baku $\mathrm{SD} \propto n^{-1/2}$.
2. **Biaya transaksi proporsional ($c > 0$)**: Setiap perubahan posisi membebankan biaya $c S_t |\delta_t - \delta_{t-1}|$. Mengikuti delta BS tanpa penyangga (*buffer*) menimbulkan *turnover* besar dan menggerus kekayaan, sedangkan strategi optimal harus membentuk zona tanpa-perdagangan (*no-trade band*) di sekitar delta target agar biaya dan risiko ekor seimbang.
3. **Pasar tak lengkap (volatilitas stokastik Heston)**: Ketika varians $v_t$ berfluktuasi secara stokastik dengan korelasi $\rho = -0{,}7$ dan *vol-of-vol* $\xi = 0{,}5$, opsi tidak dapat direplikasi sempurna hanya dengan $S_t$. Harga *indifference* bergantung pada preferensi risiko ($\rho_\lambda$) serta ketersediaan informasi $v_t$ (teramati langsung vs *proxy* varians realisasi tertunda).

Proyek ini membangun kerangka **Deep Hedging** (Buehler et al., 2019) sepenuhnya dari nol di Rust tanpa pustaka eksternal—mencakup MLP, *backpropagation through time* (BPTT) manual, Adam, solver *Dynamic Programming* (DP) eksak, serta simulator QE Heston dan *pricer* Fourier Lewis—lalu mengujinya secara ketat terhadap *anchor* analitik dan numerik pada seluruh lima bagian eksperimen.

---

## 3. Daftar Fungsi & Modul Utama

1. `rng.rs` — `splitmix64`, `derive_seed`, `Rng::normal`, `Rng::uniform`: generator deterministik `xoshiro256++` + Box–Muller per indeks lintasan.
2. `parallel.rs` — `run_chunks`, `reduce_vec_fixed_blocks`: eksekusi paralel `std::thread::scope` dengan penjumlahan blok tetap (`64` item per blok) agar hasil pelatihan bit-identik untuk `1` vs `8` *worker*.
3. `stats.rs` — `mean`, `sd`, `se`, `cvar_loss`, `entropic_risk`: statistik sampel dan risiko entropik stabil secara numerik dengan pergeseran *log-sum-exp*.
4. `nn.rs` — `MlpPolicy` (2 lapis tersembunyi lebar 32, aktivasi `tanh` atau `softplus`, keluaran linear $\delta_t \in \mathbb{R}$), `forward_step`, `backward_step`, `loss_and_d_pnl` untuk kerugian `Mse`, `Entropic`, dan `Cvar`.
5. `optim.rs` — `Adam::step` ($\beta_1 = 0{,}9, \beta_2 = 0{,}999, \varepsilon = 10^{-8}$) buatan sendiri beserta uji konvergensi kuadratik dan OLS.
6. `bs.rs` — `norm_cdf` (deret Taylor + pecahan kontinu $\approx 10^{-15}$), `bs_call`, `bs_delta`, `bs_gamma`, `implied_vol_call`, `leland_sigma`, `whalley_wilmott_half_width`, `discrete_hedging_sd_approx`.
7. `dp.rs` — `solve_entropic_dp_gbm`, `DpResult::act`: solver Bellman eksak bertransformasi eksponensial untuk kriteria entropik pada $n = 10$ dengan grid $(\ln(S/K), \delta_{t-1})$ ukuran $81 \times 81$ dan kuadratur normal 25 titik.
8. `heston.rs` — `generate_heston_qe_paths`: simulator skema *Quadratic Exponential* (Andersen, 2008) dengan koreksi *drift* martingal eksak $K_0^*$.
9. `fourier.rs` — `heston_call_fourier`, `heston_delta_fourier`: fungsi karakteristik *Little Heston Trap* (Albrecher et al., 2007) dan integral tunggal Lewis (2001) dengan aturan Simpson ($n = 1000$).
10. `hedge.rs`, `train.rs`, `eval.rs` — `FeatureSpec`, `eval_strategy`, `train_policy`, `delta_l2_vs_bs`, `measure_nn_notrade_half_width`, `HestonDeltaTable`, `paired_diff`, `ols_slope`.
11. `jsonout.rs`, `svg.rs`, `bin/run.rs` — penulis JSON 17 digit, generator 8 grafik SVG, dan orkestrasi seluruh eksperimen Bagian 1–5.

---

## 4. Uraian Kode & Desain Arsitektur

### 4.1 Jaringan Saraf, Aktivasi, dan BPTT Manual (`nn.rs`, `train.rs`)

Setiap langkah waktu $t \in \{0, \dots, n-1\}$ memetakan vektor fitur $x_t$ (seperti $\ln(S_t/K)$, waktu sisa $\tau_t = T - t\Delta t$, posisi sebelumnya $\delta_{t-1}$, serta opsional $\sqrt{v_t}$ atau $\sqrt{\widehat{\mathrm{RV}}_5}$) ke posisi $\delta_t \in \mathbb{R}$ melalui MLP $d_{\mathrm{in}} \to 32 \to 32 \to 1$.

- **Pemilihan aktivasi (`tanh`) dan keluaran linear ($\delta_t \in \mathbb{R}$)**: `tanh` dipilih karena bersifat $C^\infty$, simetris, dan memiliki turunan analitik cepat $1 - a^2$ tanpa ledakan gradien; keluaran dibiarkan tak terbatas ($\delta_t \in \mathbb{R}$) alih-alih dipaksa ke $[0, 1]$ lewat sigmoid agar jaringan tidak mengalami kejenuhan gradien (*vanishing gradient*) di daerah ITM/OTM ekstrem dan bebas menyesuaikan posisi pada skenario pasar tak lengkap atau biaya tinggi (di mana *over-hedging* atau posisi negatif kecil dapat optimal).
- **Penurunan BPTT & titik non-diferensiabel**: Untuk P&L terminal
  $$X_i = p_0 - (S_{i,n} - K)_+ + \sum_{t=0}^{n-1} \left[\delta_{i,t}(S_{i,t+1} - S_{i,t}) - c\,S_{i,t}\,|\delta_{i,t} - \delta_{i,t-1}|\right],$$
  turunan berantai mundur (*reverse-mode* sepanjang $t = n-1, \dots, 0$) menggabungkan kontribusi langsung terhadap P&L langkah $t$ dan kontribusi rekuren dari langkah $t+1$ karena $\delta_t$ menjadi input $\delta_{t-1}$ pada langkah berikutnya:
  $$\frac{\partial \mathcal{L}}{\partial \delta_t} = \frac{\partial \mathcal{L}}{\partial X_i}\bigl[(S_{t+1} - S_t) - c\,S_t\,\operatorname{sgn}(\delta_t - \delta_{t-1})\bigr] + \frac{\partial \mathcal{L}}{\partial \delta_{t-1}^{(t+1)}} + \frac{\partial \mathcal{L}}{\partial X_i} c\,S_{t+1}\,\operatorname{sgn}(\delta_{t+1} - \delta_t).$$
  Pada titik non-diferensiabel biaya $|\delta_t - \delta_{t-1}|$ di $\delta_t = \delta_{t-1}$, digunakan subgradien $\operatorname{sgn}(0) = +1$, sedangkan untuk CVaR 95% digunakan subgradien rata-rata ekor kerugian terburuk $k = \lceil (1-\alpha)R \rceil$.

### 4.2 Paralelisme Deterministik (`parallel.rs`)

Pelatihan mini-batch dijalankan dalam dua fase deterministik per iterasi Adam: (1) fase *forward* paralel mengumpulkan vektor $X$ untuk seluruh $B = 1024$ lintasan guna menghitung kerugian dan vektor $\partial \mathcal{L}/\partial X_i$ (yang saling bergantung lintas lintasan pada `Entropic` dan `Cvar`); (2) fase *backward* menggunakan `reduce_vec_fixed_blocks` dengan ukuran blok tetap `64` lintasan dan penjumlahan blok berurutan, menjamin hasil pelatihan `20` iterasi pada `1` vs `8` *worker* identik hingga bit IEEE-754 terakhir.

---

## 5. Metode Eksperimen & Parameter

1. **Bagian 1 (Pustaka & Cek Gradien)**: Verifikasi gradien BPTT vs beda hingga pusat + ekstrapolasi Richardson ($h = 10^{-5}$) pada bobot acak untuk MSE, Entropic ($\lambda = 10$), dan CVaR 95%; verifikasi Adam pada kuadratik 2D dan OLS 6 titik; uji invariansi `1` vs `8` *worker* selama 20 iterasi.
2. **Bagian 2 (Black–Scholes Tanpa Biaya)**: $S_0 = 1, K = 1, r = 0, \sigma = 0{,}2, \mu = 0, T = 0{,}25, n = 30$. Pelatihan pada $N_{\mathrm{train}} = 16.384$ lintasan (`360` iterasi Adam, `batch_size = 1024`, `lr = 3e-3` turun $0{,}5\times$ tiap `120` iterasi) dan evaluasi *out-of-sample* pada $R = 100.000$ lintasan independen. Termasuk sensitivitas `5` *seed* inisialisasi dan sapuan $n \in \{10, 20, 30, 60\}$.
3. **Bagian 3 (Biaya Transaksi & Anchor DP)**: Sapuan $c \in \{0{,}0005; 0{,}001; 0{,}002\}$ dan $\lambda \in \{10, 50\}$; pengukuran setengah-lebar pita *no-trade* NN pada $(t=15, S=1)$ vs Whalley–Wilmott $H_{\mathrm{WW}}$; *anchor* DP eksak pada $n=10, c=0{,}001, \lambda=10$; serta perbandingan strategi CRN (BS delta, WW, Leland, NN) pada $R = 100.000$.
4. **Bagian 4 (Heston QE & Fourier)**: $(\kappa, \theta, \xi, \rho, v_0) = (2; 0{,}04; 0{,}5; -0{,}7; 0{,}04)$, $R = 100.000$. Cek martingal $\mathbb{E}[S_T] = 1$ dan harga *call* MC vs Lewis Fourier ($3\ \mathrm{SE}$); perbandingan 5 strategi (*BS ATM IV*, *Heston Fourier delta* via tabel $(t, x, v)$ terkalibrasi, *NN tanpa $v_t$*, *NN dengan $v_t$*, *NN dengan RV 5-langkah*).
5. **Bagian 5 (Salah Spesifikasi & Ablasi)**: Uji silang `BS -> Heston`, `BS -> BS(σ=0,3)`, `Heston -> BS(σ=0,2)`, biaya salah $c_{\mathrm{latih}} \neq c_{\mathrm{uji}}$, serta ablasi lebar jaringan `4` dan penghapusan input $\delta_{t-1}$.

---

## 6. Target & Kriteria A Priori

Seluruh 10 kriteria berikut ditetapkan sebelum eksekusi dan dievaluasi secara otomatis di [results.json](../../outputs/day004j/results.json):

1. `c1_grad_check_and_adam_under_1e6`: Galat relatif gradien BPTT vs Richardson `< 1e-6` (MSE, Entropic, CVaR) dan galat Adam vs solusi eksak `< 1e-6`.
2. `c2_worker_invariance_1_vs_8`: Bobot setelah 20 iterasi pelatihan identik bit-demi-bit antara `1` dan `8` *worker*.
3. `c3_bs_nn_delta_l2_and_price`: Pada BS tanpa biaya ($R=100\mathrm{k}$), `L2_mse <= 0,03`, `L2_ent <= 0,03`, $|p_{0,\mathrm{MSE}} - C_{\mathrm{BS}}| \le 0{,}002$, dan $\mathrm{SD}_{\mathrm{NN,MSE}} \le \mathrm{SD}_{\mathrm{BS}} + 3\,\mathrm{SE}$.
4. `c4_discrete_hedging_n_scaling`: Eksponen skala $\log \mathrm{SD}_{\mathrm{NN}}$ terhadap $\log n$ berada di `[-0,65; -0,35]` (teori $-0{,}5$).
5. `c5_cost_notrade_band_and_ww`: Pada $c = 0{,}001, \lambda = 10$, muncul zona tanpa-perdagangan terukur dengan setengah-lebar `hw_nn > 0,01` dan eksponen skala terhadap $c$ dilaporkan terhadap $1/3$.
6. `c6_dp_anchor_n10_agrees`: Harga *indifference* NN di $n=10, c=0{,}001, \lambda=10$ cocok dengan DP eksak dalam $\max(0{,}0015,\, 3\,\mathrm{SE})$.
7. `c7_nn_beats_bs_under_cost`: Pada $c = 0{,}001, \lambda = 10$, $\rho_{10}(\mathrm{NN}) < \rho_{10}(\mathrm{BS})$.
8. `c8_heston_martingale_and_fourier_3se`: $\mathbb{E}[S_T] = 1$ dan harga *call* MC Heston QE cocok dengan Fourier Lewis dalam $3\,\mathrm{SE}$.
9. `c9_heston_vt_reduces_residual_risk`: Menambahkan $v_t$ menurunkan $\mathrm{SD}$ dan $\rho_{10}$ dibanding NN tanpa $v_t$.
10. `c10_runtime_under_600s`: Total waktu eksekusi runner penuh `< 600s`.

---

## 7. Status Pencapaian Berbasis Bukti

Seluruh runner penuh selesai dalam **`132,64s`** pada `8` *worker* (`part2 = 53,87s`, `part3 = 40,43s`, `part4 = 19,44s`, `part5 = 18,90s`), dengan `run_completed_without_errors = true`, **8 kriteria `tercapai`** dan **2 kriteria `belum tercapai`** (`c3` dan `c4`):

| Kriteria | Status | Nilai Terukur | Ambang A Priori |
| --- | --- | --- | --- |
| `c1_grad_check_and_adam_under_1e6` | `tercapai` | Semua cek unit lulus (`< 1e-6`) | `< 1e-6` |
| `c2_worker_invariance_1_vs_8` | `tercapai` | 20 iterasi `1` vs `8` worker identik bit | Bit-identik |
| `c3_bs_nn_delta_l2_and_price` | `belum tercapai` | `p0_diff = 0,000047`; `L2_mse = 0,0880`; `SD_nn = 0,01059` vs `SD_bs = 0,00624` | `L2 <= 0,03`, `SD <= 0,00629` |
| `c4_discrete_hedging_n_scaling` | `belum tercapai` | `slope_BS = -0,4817`; `slope_NN = -0,1412` | `slope_NN in [-0,65; -0,35]` |
| `c5_cost_notrade_band_and_ww` | `tercapai` | `hw_nn = 0,0286` (`c=0,001, lam=10`); eksponen `c = 0,2925` | `hw_nn > 0,01`; WW `1/3` |
| `c6_dp_anchor_n10_agrees` | `tercapai` | `DP_indiff = 0,04191` vs `NN_sim = 0,04202` (`diff = -0,000105`) | `<= max(0,0015; 3 SE)` |
| `c7_nn_beats_bs_under_cost` | `tercapai` | `rho10_NN = 0,002137` vs `rho10_BS = 0,002393` | `rho10_NN < rho10_BS` |
| `c8_heston_martingale_and_fourier_3se` | `tercapai` | `E[S_T] = 1,00015 ± 0,00031` (`z=0,51`); `MC 0,03769` vs `Fourier 0,03773` (`z=0,31`) | `<= 3 SE` |
| `c9_heston_vt_reduces_residual_risk` | `tercapai` | `SD: 0,01563 (tanpa v) -> 0,01356 (v_t)`; `rho10: 0,00113 -> 0,00086` | `v_t` menurunkan risiko |
| `c10_runtime_under_600s` | `tercapai` | `132,64s` (`8` *worker*) | `< 600s` |

### Rincian Hasil Terukur per Bagian

- **Bagian 2 (BS tanpa biaya)**: Harga BS eksak `0,039878`; aproksimasi SD diskret yang diturunkan `0,006307` cocok hingga `1%` terhadap SD MC BS `0,006244 ± 0,000014` (`rasio = 0,9900`). Harga NN MSE `p0 = 0,039831` (`diff = 0,000047`). Sensitivitas 5 *seed* inisialisasi: `L2 in [0,0918; 0,0988]`, `SD = 0,01125 ± 0,00023`.
- **Bagian 3 (Biaya & DP)**: Pada `c = 0,001, lam = 10`: $\rho_{10}(\text{WW}) = 0{,}001645$ (biaya `0,00100`, *turnover* `0,997`) < $\rho_{10}(\text{NN}) = 0{,}002137$ (biaya `0,00138`, *turnover* `1,383`) < $\rho_{10}(\text{Leland}) = 0{,}002361$ < $\rho_{10}(\text{BS}) = 0{,}002393$ (biaya `0,00223`, *turnover* `2,227`). Eksponen skala pita NN terhadap $c$ pada $\lambda = 10$ adalah `0,2925` (mendekati `1/3`), meski rasio terhadap asimtotik Whalley–Wilmott berkisar `0,15..0,17` (karena $n=30$ diskret dan $\lambda=10$ tidak berada di limit kontinu $c \to 0$). Pada $n=10$, selisih utilitas simulasi CRN `NN - DP = -0,000105 ± 0,000037` dan kesepakatan kebijakan dalam ½ langkah grid (`0,006875`) adalah `14,63%` (`MAE = 0,0245`) karena kurva kebijakan `tanh` halus di sudut pita *no-trade* yang tajam pada DP.
- **Bagian 4 (Heston)**: SD tak ter lindung nilai `0,04803`; *BS ATM IV* (`18,92%`) `SD = 0,01193, rho10 = 0,000709, CVaR95 = 0,03049`; *Heston Fourier delta* `SD = 0,01387, rho10 = 0,000936, CVaR95 = 0,03146`; *NN tanpa v* `SD = 0,01563, rho10 = 0,001130`; *NN dengan RV5* `SD = 0,01385, rho10 = 0,000913, CVaR95 = 0,03038`; *NN dengan $v_t$* `SD = 0,01356, rho10 = 0,000855, CVaR95 = 0,02900` dengan *indifference price* `0,03859` (> harga risiko-netral Fourier `0,03773` sebesar `+0,00086` sebagai premi risiko tak lengkap).
- **Bagian 5 (Salah Spesifikasi & Ablasi)**:
  - `BS NN -> Heston`: `rho10 = 0,00182` (degradasi `+61%` dibanding `Heston NN 0,00113`); `BS NN -> BS(0,3)`: `rho10 = 0,00339` vs `BS delta(0,2)`: `0,00086`; `Heston NN -> BS(0,2)`: `rho10 = 0,00104`.
  - Biaya salah: `c_train=0,0005` diuji pada `c=0,002` memberi `rho10 = 0,00407` (terlalu banyak *turnover*); `c_train=0,002` diuji pada `c=0,0005` memberi `rho10 = 0,00162`.
  - Ablasi pada `c = 0,001, lam = 10`: lebar `4` memberi `rho10 = 0,00228` (`hw = 0,0214`), sedangkan penghapusan input $\delta_{t-1}$ menghapus pita *no-trade* (`hw = 0,0071`) dan melonjakkan `rho10` menjadi `0,00519` (`+143%`).

---

## 8. Pemisahan Fakta vs Hipotesis (Diagnosis `c3` dan `c4`)

- **FAKTA**:
  1. Gradien BPTT analitik dan optimizer Adam terverifikasi eksak (`< 1e-6`), dan $p_{0,\mathrm{MSE}} = 0{,}039831$ cocok hingga `4,7e-5` terhadap harga BS `0,039878`.
  2. Namun pada `c3`, `L2_mse = 0,0880 > 0,03` (`sup = 0,426` di ujung rentang `[0,8; 1,2]`) dan `SD_nn = 0,01059 > SD_bs = 0,00624`; akibat lantai galat optimisasi `~0,010` yang muncul merata di $n \in \{10, 20, 30, 60\}$ (`0,0135 -> 0,0117 -> 0,0111 -> 0,0105`), kemiringan skala `c4` NN hanya `-0,141` (sementara BS delta sendiri memberi `-0,482`, tepat sesuai teori `-0,5`).
- **HIPOTESIS**: Dengan `360` iterasi Adam (`batch_size = 1024`, `lr = 3e-3 -> 7,5e-4`), jaringan cukup cepat belajar harga rata-rata dan strategi pita biaya (`c5, c6, c7` tercapai), tetapi belum cukup konvergen pada gradien orde kedua di wilayah ekor $\ln(S/K) \approx \pm 0{,}2$ dan dekat jatuh tempo $\tau \to 0$ (di mana delta BS berbentuk fungsi tangga curam), sehingga lantai galat aproksimasi/optimisasi mendominasi di atas galat diskretisasi $n^{-1/2}$. Sesuai aturan integritas, ambang `c3` dan `c4` tidak dilonggarkan dan dilaporkan `belum tercapai` apa adanya.

---

## 9. Landasan Teori (10 Poin Terkait Kode)

1. **Deep Hedging & Risiko Entropik (Buehler et al., 2019)**: Kebijakan lindung nilai diparameterisasi oleh jaringan saraf berinput fitur keadaan dan posisi sebelumnya untuk meminimalkan $\rho_\lambda(X) = \frac{1}{\lambda}\ln\mathbb{E}[e^{-\lambda X}]$; di [nn.rs](../../projects/rust/day004j_deep_hedging/src/nn.rs) dan [train.rs](../../projects/rust/day004j_deep_hedging/src/train.rs).
2. **Invariansi Kas & *Indifference Price***: $\rho_\lambda(p_0 - Z) = \rho_\lambda(-Z) - p_0$, sehingga harga *indifference* $p_0^* = \rho_\lambda(-Z)$ memisahkan penentuan harga dari optimisasi kebijakan $\delta_t$; di [hedge.rs](../../projects/rust/day004j_deep_hedging/src/hedge.rs).
3. **Backpropagation Through Time (BPTT) pada Kebijakan Rekuren**: Karena biaya transaksi menggandeng $\delta_t$ dengan $\delta_{t-1}$, turunan total terhadap bobot harus merambatkan sensitivitas posisi masa depan mundur dalam waktu; di [nn.rs](../../projects/rust/day004j_deep_hedging/src/nn.rs) dan [train.rs](../../projects/rust/day004j_deep_hedging/src/train.rs).
4. **Optimisasi Adam & Ekstrapolasi Richardson**: Momen pertama/kedua adaptif Adam ($\beta_1=0{,}9, \beta_2=0{,}999$) diverifikasi pada fungsi konveks, sementara gradien BPTT diverifikasi terhadap beda hingga orde-4 Richardson; di [optim.rs](../../projects/rust/day004j_deep_hedging/src/optim.rs) dan [nn.rs](../../projects/rust/day004j_deep_hedging/src/nn.rs).
5. **Galat Lindung Nilai Diskret (Boyle–Emanuel / Bertsimas–Kogan–Lo)**: Varians P&L *delta-hedging* diskret menyusut sebagai $O(1/n)$ dengan $\mathrm{SD} \approx \sqrt{\frac{T}{2n}\int_0^T \mathbb{E}[(S_t^2 \Gamma_t \sigma^2)^2]\,dt}$; di [bs.rs](../../projects/rust/day004j_deep_hedging/src/bs.rs).
6. **Asimtotik Pita Tanpa-Perdagangan Whalley–Wilmott (1997)**: Pertukaran biaya $c$ dan varians risiko lokal menghasilkan setengah-lebar pita $H = (3 e^{-r\tau} c S \Gamma^2 / (2\lambda))^{1/3}$; di [bs.rs](../../projects/rust/day004j_deep_hedging/src/bs.rs) dan [eval.rs](../../projects/rust/day004j_deep_hedging/src/eval.rs).
7. **Volatilitas Tersesuaikan Leland (1985)**: Penyesuaian $\sigma_{\mathrm{Leland}} = \sigma\sqrt{1 + \sqrt{8/\pi}\,c/(\sigma\sqrt{\Delta t})}$ untuk biaya satu arah $c$; di [bs.rs](../../projects/rust/day004j_deep_hedging/src/bs.rs).
8. **Pemrograman Dinamis Eksponensial CARA**: Transformasi $V_t = \frac{1}{\lambda}\ln J_t$ mengubah ekspektasi utilitas multiplikatif menjadi rekursi Bellman aditif di ruang $(\ln S, \delta_{t-1})$; di [dp.rs](../../projects/rust/day004j_deep_hedging/src/dp.rs).
9. **Skema Quadratic Exponential Andersen (2008)**: Simulasi varians CIR tanpa nilai negatif dengan koreksi martingal $K_0^*$ eksak pada langkah log-spot; di [heston.rs](../../projects/rust/day004j_deep_hedging/src/heston.rs).
10. **Integral Lewis (2001) & *Little Heston Trap* (Albrecher et al., 2007)**: Harga opsi Heston bebas diskontinuitas cabang lewat kontur $\operatorname{Im}(z) = -1/2$ dan integrasi Simpson; di [fourier.rs](../../projects/rust/day004j_deep_hedging/src/fourier.rs).
11. **Conditional Value-at-Risk (CVaR 95%)**: Ukuran risiko koheren $\mathrm{CVaR}_\alpha(-X)$ yang mengukur rata-rata kerugian pada ekor terburuk $5\%$ beserta subgradiennya; di [stats.rs](../../projects/rust/day004j_deep_hedging/src/stats.rs) dan [nn.rs](../../projects/rust/day004j_deep_hedging/src/nn.rs).
12. **Common Random Numbers (CRN) & Selisih Berpasangan**: Penggunaan lintasan identik antar strategi memangkas varians pembanding secara tajam sehingga perbedaan kinerja $\pm\,\mathrm{SE}$ dapat terukur akurat; di [eval.rs](../../projects/rust/day004j_deep_hedging/src/eval.rs).
13. **Pasar Tak Lengkap & Lindung Nilai Minimum-Variance**: Di model Heston tanpa perdagangan opsi, korelasi $\rho < 0$ menciptakan komponen lindung nilai silang $\frac{\rho\xi}{S}\frac{\partial C}{\partial v}$ yang dapat dipelajari jaringan bila $v_t$ diamati; di [hedge.rs](../../projects/rust/day004j_deep_hedging/src/hedge.rs) dan [bin/run.rs](../../projects/rust/day004j_deep_hedging/src/bin/run.rs).

---

## 10. Kesimpulan, Batasan, dan Deviasi dari Prompt

1. **Kesimpulan**: Seluruh mesin Deep Hedging dari nol di Rust berjalan dalam `132,64s` (`8` *worker*) dengan determinisme bit-identik, gradien BPTT eksak (`< 1e-6`), kesesuaian DP eksak pada $n=10$ (`diff = -0,000105`), pembentukan pita *no-trade* dengan eksponen skala biaya `0,2925`, validasi Heston QE vs Fourier (`z = 0,31`), serta pembuktian bahwa informasi $v_t$ memangkas risiko residual pasar tak lengkap (`SD: 0,01563 -> 0,01356`).
2. **Batasan**: Seluruh data sintetis; anggaran pelatihan dibatasi `360` iterasi Adam agar total waktu runner jauh di bawah batas `600s`, yang menyebabkan `c3` dan `c4` berstatus `belum tercapai`.
3. **Deviasi dari Prompt**:
   - Pada Bagian 1, `lambda = 1` dicakup oleh formulasi umum $\rho_\lambda$ tetapi sapuan utama Bagian 3 menggunakan $\lambda \in \{10, 50\}$ (karena pada $\lambda = 1$ kurvanya hampir identik dengan netral-risiko pada horizon $0{,}25$ tahun).
   - Pada Bagian 4, *Heston Fourier delta* dihitung melalui tabel interpolasi bilinear $(t, \ln(S/K), v)$ berukuran `30 x 45 x 25` yang dibangun langsung dari `heston_delta_fourier` agar evaluasi $R = 100.000$ lintasan tetap dalam anggaran waktu `< 10` menit.
   - Kriteria `c3` dan `c4` tetap dipertahankan dengan ambang ketatnya dan dilaporkan `belum tercapai` tanpa penyetelan parameter pasca-uji.

---

## 11. Referensi dan Jejak Bukti

- Andersen, L. (2008). *Simple and efficient simulation of the Heston stochastic volatility model*. Journal of Computational Finance, 11(3), 1–42.
- Buehler, H., Gonon, L., Teichmann, J., & Wood, B. (2019). *Deep hedging*. Quantitative Finance, 19(8), 1271–1291.
- Lewis, A. L. (2001). *A simple option formula for general jump-diffusion and other exponential Lévy processes*. Envision Financial Systems.
- Whalley, A. E., & Wilmott, P. (1997). *An asymptotic analysis of an optimal hedging model for option pricing with transaction costs*. Mathematical Finance, 7(3), 307–324.
- Artefak hasil & tes: [outputs/day004j/results.json](../../outputs/day004j/results.json), [README.md](../../projects/rust/day004j_deep_hedging/README.md), [integration_tests.rs](../../projects/rust/day004j_deep_hedging/tests/integration_tests.rs).
- Grafik SVG (`outputs/figures/`): [learning-curve](../../outputs/figures/day004j_learning-curve.svg), [nn-vs-bs-delta](../../outputs/figures/day004j_nn-vs-bs-delta.svg), [pnl-distribution](../../outputs/figures/day004j_pnl-distribution.svg), [notrade-band-vs-ww](../../outputs/figures/day004j_notrade-band-vs-ww.svg), [dp-vs-nn-n10](../../outputs/figures/day004j_dp-vs-nn-n10.svg), [heston-strategies](../../outputs/figures/day004j_heston-strategies.svg), [ablation-width4-noprev](../../outputs/figures/day004j_ablation-width4-noprev.svg), [misspecification](../../outputs/figures/day004j_misspecification.svg).
