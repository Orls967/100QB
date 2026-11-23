# Penjelasan Proyek 100QB

Folder ini menyimpan dokumen penjelasan teknis berbahasa Indonesia untuk setiap sesi/proyek di repositori **100QB**, mengikuti spesifikasi struktur 11 bagian (identitas, masalah, fungsi, uraian kode, metode, target, status pencapaian berbasis bukti, pemisahan fakta vs hipotesis, minimal 10 landasan teori, kesimpulan, dan jejak bukti).

## Daftar Proyek (Urut Hari & Sesi)

| Hari / Sesi | Topik No | Jenis | Bahasa | Dokumen Penjelasan |
|:---:|:---:|:---:|:---:|---|
| `1a` (`day001a`) | `001` | `quant` | `python` | [Return harian & kumulatif](day001a_quant_return-harian-kumulatif.md) |
| `1b` (`day001b`) | `002` | `sim` | `python` | [Estimasi π dengan Monte Carlo](day001b_sim_estimasi-dengan-monte-carlo.md) |
| `1c` (`day001c`) | `bonus` | `quant` | `python` | [Opsi Asia: Monte Carlo, variance reduction & Greeks (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) |
| `1d` (`day001d`) | `bonus` | `quant` | `rust` | [Opsi Asia (Rust): Monte Carlo, variance reduction & Greeks](day001d_quant_asian-options.md) |
| `1e` (`day001e`) | `bonus` | `quant` | `go` | [Opsi Asia (Go): Monte Carlo, variance reduction & Greeks](day001e_quant_asian-options.md) |
| `1f` (`day001f`) | `bonus` | `quant` | `typescript` | [Opsi Asia (TypeScript): Monte Carlo, variance reduction & Greeks](day001f_quant_asian-options.md) |
| `1g` (`day001g`) | `bonus` | `quant` | `cpp` | [Opsi Asia (C++): Monte Carlo, variance reduction & Greeks](day001g_quant_asian-options.md) |
| `1h` (`day001h`) | `bonus` | `quant` | `c` | [Opsi Asia (C): Monte Carlo, variance reduction & Greeks](day001h_quant_asian-options.md) |
| `1i` (`day001i`) | `bonus` | `quant` | `ruby` | [Opsi Asia (Ruby): Monte Carlo, variance reduction & Greeks](day001i_quant_asian-options.md) |
| `2b` (`day002b`) | `bonus` | `quant` | `python` | [Estimasi kovarians & portofolio varians minimum: Wishart, Marchenko–Pastur, Ledoit–Wolf](day002b_quant_estimasi-kovarians-portofolio-varians.md) |
| `2c` (`day002c`) | `bonus` | `quant` | `python` | [HMM dua rezim: forward-backward, Baum–Welch & Viterbi](day002c_quant_hmm-dua-rezim-forward-backward-baum.md) |
| `2d` (`day002d`) | `bonus` | `quant` | `python` | [Risiko kredit Vasicek: distribusi kerugian eksak, LHP & importance sampling](day002d_quant_risiko-kredit-vasicek-distribusi.md) |
| `2e` (`day002e`) | `bonus` | `quant` | `python` | [Model Heston: pricing Fourier, kalibrasi & Monte Carlo QE](day002e_quant_model-heston-pricing-fourier-kalibrasi.md) |
| `2f` (`day002f`) | `bonus` | `quant` | `rust` | [Inferensi Bayesian state-space: Kalman, particle filter & PMMH](day002f_quant_state-space-pmmh.md) |
| `2g` (`day002g`) | `bonus` | `quant` | `cpp` | [AAD untuk Greeks Monte Carlo: tape mode-balik, pathwise & keranjang 10 aset](day002g_quant_aad-greeks.md) |
| `2h` (`day002h`) | `bonus` | `quant` | `cpp` | [Hawkes + limit order book: simulasi eksak, MLE, mesin pencocokan & fakta stilistik](day002h_quant_hawkes-orderbook.md) |
| `2i` (`day002i`) | `bonus` | `quant` | `cpp` | [Hull–White: pohon trinomial, Jamshidian, kalibrasi & Bermudan](day002i_quant_hull-white.md) |
| `2j` (`day002j`) | `bonus` | `quant` | `ruby` | [Longevity risk: Lee–Carter, anuitas stokastik & modal solvabilitas](day002j_quant_longevity-risk.md) |
| `2k` (`day002k`) | `bonus` | `quant` | `ruby` | [Opsi Put Amerika (Ruby): CRR Binomial, LSM Monte Carlo & PSOR PDE](day002k_quant_american-put.md) |
| `3a` (`day003a`) | `bonus` | `quant` | `python` | [Volatilitas lokal Dupire dari permukaan SSVI](day003a_quant_volatilitas-lokal-dupire-ssvi.md) |
| `3b` (`day003b`) | `bonus` | `quant` | `rust` | [Eksekusi optimal: Almgren–Chriss, dampak transien, likuiditas stokastik](day003b_quant_optimal-execution.md) |
| `3c` (`day003c`) | `bonus` | `quant` | `go` | [Risiko sistemik jaringan: kliring Eisenberg–Noe, kaskade Watts & biaya kebangkrutan](day003c_quant_systemic-risk.md) |
| `3d` (`day003d`) | `bonus` | `quant` | `typescript` | [Risiko ekor: GARCH–EVT, VaR/ES & backtesting](day003d_quant_tail-risk.md) |
| `3e` (`day003e`) | `bonus` | `quant` | `go` | [CVA swap suku bunga: Hull–White, intensitas CIR dan wrong-way risk](day003e_quant_cva-wrong-way.md) |
| `3f` (`day003f`) | `bonus` | `quant` | `go` | [Teori ruin Cramér–Lundberg dan simulasi kejadian langka](day003f_quant_ruin-theory.md) |
| `4a` (`day004a`) | `bonus` | `quant` | `go` | [Struktur suku bunga: Nelson–Siegel, AFNS dan filter Kalman](day004a_quant_term-structure.md) |
| `4b` (`day004b`) | `bonus` | `quant` | `go` | [AFNS: bias sampel terbatas, rentang kalender (long-span), dan frekuensi observasi (infill)](day004b_quant_afns-inference.md) |
| `4c` (`day004c`) | `bonus` | `quant` | `go` | [AFNS: koreksi bias mean-reversion dan inferensi sampel terbatas (jackknife & bootstrap)](day004c_quant_afns-bias-correction.md) |
| `4d` (`day004d`) | `bonus` | `quant` | `go` | [AFNS: kalibrasi interval bootstrap yang peka bias dan batas parameter](day004d_quant_afns-interval-calibration.md) |
| `4e` (`day004e`) | `bonus` | `quant` | `go` | [Harga listrik negatif dan arbitrase baterai (DP stokastik)](day004e_quant_power-storage.md) |
| `4f` (`day004f`) | `bonus` | `quant` | `rust` | [Port lengkap AFNS term structure, kalibrasi interval, harga listrik & valuasi baterai (Rust)](day004f_rust-afns-energy-storage.md) |
| `4g` (`day004g`) | `bonus` | `quant` | `rust` | [Market making optimal: Guéant–Lehalle–Tapia, HJB dan salah spesifikasi (Rust)](day004g_quant_market-making.md) |
| `4h` (`day004h`) | `bonus` | `quant` | `rust` | [Volatilitas kasar: rough Bergomi, skema hybrid, skew ATM & identifiabilitas Hurst (Rust)](day004h_quant_rough-bergomi.md) |
| `4i` (`day004i`) | `bonus` | `quant` | `rust` | [Cointegration & trading OU optimal: Dickey–Fuller, Engle–Granger, Johansen & ambang Bertram (Rust)](day004i_quant_pairs-trading.md) |
| `4j` (`day004j`) | `bonus` | `quant` | `rust` | [Deep hedging dari nol: jaringan saraf, risiko entropik, biaya transaksi & pasar tak lengkap (Rust)](day004j_quant_deep-hedging.md) |
| `4k` (`day004k`) | `bonus` | `quant` | `rust` | [AMM: impermanent loss, LVR, biaya & likuiditas terkonsentrasi (Rust)](day004k_quant_amm-lvr.md) |
| `5a` (`day005a`) | `bonus` | `quant` | `rust` | [Model SABR: CEV eksak, ekspansi Hagan, kepadatan negatif, kalibrasi & delta Bartlett (Rust)](day005a_quant_sabr.md) |
| `5b` (`day005b`) | `bonus` | `quant` | `rust` | [Harga opsi Eropa Black–Scholes & put-call parity (Rust)](day005b_quant_black-scholes.md) |
| `5c` (`day005c`) | `bonus` | `quant` | `rust` | [Kriteria Kelly dan risiko drawdown: peluang sentuh eksak, Kelly fraksional, galat estimasi & ekor berat (Rust)](day005c_quant_kelly-drawdown.md) |

## Sesi yang Tidak Ada di Repositori

Sesuai pemeriksaan artefak repositori (`days/`, folder bahasa `projects/rust/`, `projects/go/`, `projects/typescript/`, `projects/cpp/`, `projects/c/`, `projects/ruby/`, serta `outputs/`), sesi berikut **tidak pernah ada** di dalam repositori ini dan oleh karena itu tidak dibuatkan dokumen penjelasan:
- **`day001j`, `day001k`, `day001l`** (implementasi Opsi Asia dalam Java, Kotlin, atau JavaScript murni): seri lintas bahasa Hari 1 berhenti pada `day001i` (Ruby).
- **`day002a`** (notebook/modul Python sesi pertama Hari 2, misalnya LSM Python): Hari 2 di repositori langsung dimulai dari `day002b` (*Estimasi kovarians & portofolio varians minimum*).
