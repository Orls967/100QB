# Penjelasan Proyek `day006c`: Perkolasi Bond di Kisi Persegi — Algoritma Newman–Ziff, Dualitas-Diri Eksak, dan Eksponen Kritis (Poliglot Python + Rust + Go)

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day006c` (Hari 6, sesi `c`, topik `bonus` di luar `ROADMAP.md`)
- **Tanggal**: `2026-11-24` (WITA)
- **Topik**: Perkolasi *bond* di kisi persegi — algoritma Newman–Ziff, dualitas-diri eksak $p_c = 1/2$, dan eksponen kritis ($\nu = 4/3$, $D = 91/48$)
- **Jenis**: `sim` (simulasi fisika statistik & transisi fasa geometrik)
- **Bahasa**: `python+rust+go` (inti Newman–Ziff di Rust murni tanpa *crate*, orkestrator kolam pekerja persisten di Go pustaka standar, serta enumerasi eksak, konvolusi kanonik, regresi penskalaan, dan visualisasi di Python)
- **Lokasi Proyek**: [`projects/rust/day006c_percolation_polyglot/`](../../projects/rust/day006c_percolation_polyglot/README.md) dengan metadata [`meta.json`](../../projects/rust/day006c_percolation_polyglot/meta.json) agar pemindai [`scripts/update_readme.py`](../../scripts/update_readme.py) mengenali sesi `6c` secara otomatis tanpa mengubah skrip pemeriksa.

## 2. Masalah dan Tujuan

Perkolasi *bond* pada kisi persegi membuka setiap ikatan (*bond*) antar-simpul bertetangga secara independen dengan peluang $p \in [0, 1]$. Pada batas termodinamika $L \to \infty$, sistem mengalami transisi fasa kontinu tepat di ambang kritis $p_c = 1/2$. Pada kisi berhingga, simetri dualitas-diri muncul secara eksak hanya apabila geometri persegi panjang dipilih sehingga graf dual planar isomorfik terhadap graf primal di bawah pertukaran batas kiri–kanan dan atas–bawah. Proyek ini menjawab empat pertanyaan kuantitatif:
1. Bentuk persegi panjang simpul $a \times b$ manakah yang memberikan peluang penyeberangan kiri–kanan eksak $P_{\text{cross}}(1/2) = 1/2$ tanpa koreksi ukuran berhingga orde utama, sebagaimana dibuktikan lewat enumerasi seluruh $2^E$ subhimpunan *bond* memakai aritmetika rasional eksak?
2. Apakah algoritma mikrokanonik Newman–Ziff di Rust yang diorkestrasi secara paralel oleh Go memproduksi kurva kanonik $P_{\text{cross}}(p)$ yang cocok dengan polinomial eksak pada kisi kecil serta menjaga $|P_{\text{cross}}(1/2) - 0{,}5| \le 3\,\text{SE}$ pada seluruh $L \in \{16, 32, 64, 128\}$?
3. Seberapa akurat kemiringan kritis $\mathrm{d}P_{\text{cross}}/\mathrm{d}p(1/2) \sim L^{1/\nu}$ memulihkan eksponen panjang korelasi 2D eksak $\nu = 4/3$ ($1/\nu = 3/4 = 0{,}75$)?
4. Seberapa dekat penskalaan rata-rata ukuran klaster fisik terbesar pada setengah *bond* $n = \lfloor N/2 \rfloor$ memulihkan dimensi fraktal kritis $D = 91/48 \approx 1{,}895833$?

## 3. Apa yang Dikerjakan Proyek

Sistem menggabungkan tiga komponen lintas bahasa yang saling memverifikasi. Pertama, Python mengevaluasi seluruh $2^E$ konfigurasi *bond* pada persegi panjang kandidat (seperti $2 \times 3$ vs $3 \times 2$ dengan $E = 7$, serta $3 \times 4$ vs $4 \times 3$ dengan $E = 17 \le 20$) untuk membuktikan secara aljabar bahwa orientasi $L$ baris $\times (L+1)$ kolom simpul ($N = 2L^2 - 1$ *bond*) memenuhi $P_{\text{cross}}(1/2) = 1/2$ persis. Kedua, biner Rust `percolate` mengeksekusi pencuplikan mikrokanonik Newman–Ziff pada kisi $L \times (L+1)$ melalui struktur Union-Find ganda (satu dengan dua simpul maya kiri/kanan untuk mendeteksi jumlah *bond* penyeberang terkecil $n^*$, dan satu pada simpul fisik untuk mengukur klaster terbesar di $n = \lfloor N/2 \rfloor$) serta menyediakan mode `sample` untuk menggambar klaster kritis pada $L = 96, p = 0{,}5$. Ketiga, orkestrator Go `orchestrate` membagi pekerjaan per ukuran $L \in \{16, 32, 64, 128\}$ menjadi $20$ potongan tetap dengan *seed* `SplitMix64` per potongan yang tak bergantung jumlah pekerja, mengelola kolam proses Rust persisten lewat pipa JSON per baris di bawah batas waktu `context.WithTimeout`, menjumlahkan histogram bilangan bulat secara deterministik, dan membuktikan invarian `1` lawan banyak *worker*. Terakhir, Python mengonversi histogram mikrokanonik menjadi kurva kanonik $P_{\text{cross}}(p)$ dan turunannya $\mathrm{d}P_{\text{cross}}/\mathrm{d}p$, mengevaluasi regresi log–log $\hat{\nu}$ dan $\hat{D}$, merender grafik `PNG` + `SVG`, dan mengunci waktu dinding lewat `finalize.py`.

## 4. Uraian Kode

- [`projects/rust/day006c_percolation_polyglot/rust/src/lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs): Mengimplementasikan PRNG `Xoshiro256Pp` (diinisialisasi `splitmix64`, batas bebas-bias Lemire `bounded`, pengacakan Fisher–Yates `shuffle`), struktur `Uf` (*union by size* dan *path compression*), konstruktor kisi `build_rect_bonds`, simulator mikrokanonik `run_newman_ziff`, pencuplik konfigurasi `sample_bonds`, serta 4 tes unit Rust.
- [`projects/rust/day006c_percolation_polyglot/rust/src/main.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/main.rs): Antarmuka CLI `percolate` untuk mode pekerja pipa JSON `job` (membaca `{"L","seed","samples"}` per baris dan menulis histogram $n^*$ beserta jumlah ukuran klaster terbesar dengan `flush` eksplisit) dan mode `sample` untuk konfigurasi *bond* terbuka.
- [`projects/rust/day006c_percolation_polyglot/go/main.go`](../../projects/rust/day006c_percolation_polyglot/go/main.go) & [`orchestrator_test.go`](../../projects/rust/day006c_percolation_polyglot/go/orchestrator_test.go): Fungsi `ChunkSeed` dan `BuildJobs` membentuk potongan tetap per $(L, \text{chunk})$, `RunPool` mendistribusikan potongan ke proses Rust persisten lewat pipa `stdin`/`stdout` berpenyangga `8 MB` dengan penghentian `context.WithTimeout` lalu mereduksi `uint64` sesuai urutan indeks potongan, dilengkapi uji invarian pekerja dan 3 tes unit Go.
- [`projects/rust/day006c_percolation_polyglot/py/percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py), [`finalize.py`](../../projects/rust/day006c_percolation_polyglot/py/finalize.py), & [`test_percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/test_percolation.py): Enumerasi eksak `fractions.Fraction` (`exact_enumeration`, `discover_self_dual_shapes`), versi acuan Newman–Ziff Python (`python_newman_ziff`) untuk uji silang kisi kecil, konvolusi bobot binomial berbasis `math.lgamma` (`binom_weights`, `canonical_eval`), estimasi OLS eksponen kritis $\hat{\nu}$ dan $\hat{D}$ (`ols_loglog`), pembuatan 5 pasang grafik `PNG`/`SVG` (`make_all_figures`), pencatat waktu tunggal `finalize.py`, serta 4 tes `unittest`.

## 5. Metode dan Alur Kerja

1. **Penemuan Geometri Swa-Dual Eksak**: Untuk kisi simpul $a$ baris $\times b$ kolom, penyeberangan kiri–kanan mengharuskan lintasan terbuka dari kolom $c = 0$ ke kolom $c = b - 1$. Graf dual planar memiliki $(a + 1)$ baris wajah (termasuk wajah luar atas dan bawah) serta $(b - 1)$ kolom wajah; rotasi $90^\circ$ memetakan penghalang atas–bawah dual ke penyeberangan kiri–kanan pada kisi $(b - 1) \times (a + 1)$ simpul. Keswadualan menuntut $(a, b) = (b - 1, a + 1) \iff b = a + 1$, yaitu $L$ baris $\times (L + 1)$ kolom dengan $N = 2L^2 - 1$ *bond*.
2. **Pencuplikan Mikrokanonik Newman–Ziff**: Dalam satu realisasi permutasi seragam $e_{\pi(1)}, \dots, e_{\pi(N)}$, keadaan ke-$n$ memiliki tepat $n$ *bond* terbuka. Karena sifat penyeberangan monoton naik terhadap penambahan *bond*, indikator seberang bernilai $0$ untuk $n < n^*$ dan $1$ untuk $n \ge n^*$. Akibatnya frekuensi kumulatif $C(n) = \sum_{k=1}^n H[k]$ dibagi jumlah sampel $M$ merupakan estimator tak-bias bagi peluang seberang mikrokanonik $\hat{Q}(n) = C(n)/M$.
3. **Transformasi Kanonik & Turunan Analitik**: Peluang seberang kanonik pada peluang buka $p$ dihitung lewat $\hat{P}(p) = \sum_{n=1}^N B_{N,n}(p)\,\hat{Q}(n)$ dengan log-bobot $\ln B_{N,n}(p) = \ln\Gamma(N+1) - \ln\Gamma(n+1) - \ln\Gamma(N-n+1) + n\ln p + (N-n)\ln(1-p)$. Turunan terhadap $p$ dihitung secara eksak tanpa beda hingga melalui identitas skor binomial $\frac{\mathrm{d}}{\mathrm{d}p} B_{N,n}(p) = \frac{n - Np}{p(1-p)} B_{N,n}(p)$.

## 6. Output yang Ditargetkan

- Bukti rasional eksak `P(1/2) = 1/2` pada kisi $2 \times 3$ dan $3 \times 4$ simpul (dan ketidakcocokan pada orientasi terbalik $3 \times 2$, $4 \times 3$).
- Kesesuaian Newman–Ziff terhadap polinomial eksak pada kisi kecil di seluruh 9 titik grid $p \in \{0{,}1, 0{,}2, \dots, 0{,}9\}$ dalam $3\,\text{SE}$.
- Kriteria *a priori* pada $L \in \{16, 32, 64, 128\}$ ($M = 20.000$ untuk $L \le 64$, $M = 6.000$ untuk $L = 128$, `seed = 20261124`):
  - $|P_{\text{cross}}(1/2) - 0{,}5| \le 3\,\text{SE}$ untuk setiap $L$;
  - $|\hat{\nu} - 4/3| < 0{,}15$ dari regresi $\ln \mathrm{d}P/\mathrm{d}p(1/2)$ terhadap $\ln L$;
  - $|\hat{D} - 91/48| < 0{,}06$ dari regresi $\ln \langle S_{\max}(N/2) \rangle$ terhadap $\ln L$;
  - Validitas 5 grafik `PNG` + `SVG` (`perc_cluster`, `perc_crossing_curves`, `perc_collapse`, `perc_derivative_scaling`, `perc_largest_cluster`), invarian pekerja `1` vs `N`, dan penyelesaian sebelum batas `finalize.py` (`e <= 570` detik).

## 7. Status Pencapaian

Tabel berikut menggunakan nama kriteria yang persis sama dengan [`outputs/day006c/results.json`](../../outputs/day006c/results.json):

| Kriteria | Status | Nilai Terukur | Ambang A Priori & Bukti |
|---|---|---|---|
| `c1_exact_self_dual_half` | `tercapai` | `2x3`: `1/2`, `3x4`: `1/2`, `3x2`: `7/8` (`4x3`: `201/256`) | `P(1/2) == 1/2` eksak hanya pada $L \times (L+1)$ (enumerasi rasional Python & `cargo test --release`) |
| `c2_nz_vs_exact_small_lattice` | `tercapai` | maks `\|Z\| = 0,4968` ($3 \times 4$, $M=15.000$, 9 titik $p$) | $\le 3\,\text{SE}$; terverifikasi di `percolation.py` & `test_nz_python_vs_exact` |
| `c3_pcross_half_all_L` | `tercapai` | maks `\|Z\| = 1,4839` ($L=128$: $0{,}5096 \pm 0{,}0065$; $L=16..64 \le 1{,}29$) | $\le 3\,\text{SE}$ untuk seluruh $L \in \{16,32,64,128\}$ |
| `c4_nu_hat` | `tercapai` | `nu_hat = 1,340556` ($\text{SE} = 0{,}007905$, `slope = 0,745959`) | $\|\hat{\nu} - 4/3\| = 0{,}0072 < 0{,}15$ via regresi log–log `ols_loglog` |
| `c5_fractal_dim_D_hat` | `tercapai` | `d_hat = 1,867779` ($\text{SE} = 0{,}007694$, acuan $91/48 \approx 1{,}895833$) | $\|\hat{D} - 91/48\| = 0{,}0281 < 0{,}06$ |
| `c6_figures_valid` | `tercapai` | `5` pasang `PNG`+`SVG` valid | Semua file non-kosong & lolos `ET.parse` (`day006c_perc_*.{png,svg}`) |
| `c7_worker_invariance` | `tercapai` | `8` *worker* identik bit-demi-bit | `worker_invariant: true` di `run_all.sh` & `TestDeterministicReductionOrder` (`go test`) |
| `c_wall_clock` | `belum tercapai` | `789` detik (`elapsed_sec` sesi asli) | Melewati batas `elapsed_sec <= 570` pada sesi asli 10 menit (dipertahankan jujur; *rerun* pasca-sesi 3 detik dicatat terpisah di `rerun_note`) |

Grafik keluaran (`PNG` + `SVG`):
- Klaster terbesar $L=96$: [day006c_perc_cluster.svg](../../outputs/figures/day006c_perc_cluster.svg) ([PNG](../../outputs/figures/day006c_perc_cluster.png))
- Kurva penyeberangan $P_{\text{cross}}(p)$: [day006c_perc_crossing_curves.svg](../../outputs/figures/day006c_perc_crossing_curves.svg) ([PNG](../../outputs/figures/day006c_perc_crossing_curves.png))
- *Data collapse* ($\nu = 4/3$): [day006c_perc_collapse.svg](../../outputs/figures/day006c_perc_collapse.svg) ([PNG](../../outputs/figures/day006c_perc_collapse.png))
- Penskalaan turunan $\mathrm{d}P/\mathrm{d}p(1/2)$: [day006c_perc_derivative_scaling.svg](../../outputs/figures/day006c_perc_derivative_scaling.svg) ([PNG](../../outputs/figures/day006c_perc_derivative_scaling.png))
- Penskalaan klaster terbesar $\langle S_{\max}\rangle$: [day006c_perc_largest_cluster.svg](../../outputs/figures/day006c_perc_largest_cluster.svg) ([PNG](../../outputs/figures/day006c_perc_largest_cluster.png))

## 8. Jika Belum Tercapai atau Belum Diketahui

### Pemisahan FAKTA vs HIPOTESIS
- **FAKTA (Rerun Pasca-Sesi — `day006d` Bagian A)**: Setelah sesi 10 menit berakhir, seluruh rangkaian tes dan *pipeline* penuh dijalankan ulang tanpa mengubah parameter, seed, maupun ambang:
  - `cargo test --release --manifest-path rust/Cargo.toml`: `4/4` lulus.
  - `go vet ./... && go test -v ./...`: `3/3` lulus (`0,385s`).
  - `python -m unittest -v test_percolation`: `4/4` lulus (`0,428s`).
  - `bash run_all.sh`: exit `0` dalam `3` detik, menghasilkan `outputs/day006c/results.json` dan `5` pasang grafik `outputs/figures/day006c_perc_*.{png,svg}`.
- **FAKTA (Batas Waktu Dinding Sesi Asli)**: `c_wall_clock` tetap `belum tercapai` (`789` detik vs `570` detik) karena pada sesi asli 10 menit `dl.sh` memasuki status `LOCK` (`>= 540` detik) sebelum `run_all.sh` dieksekusi. Nilai `timing` asli dipertahankan utuh, dan pelaksanaan ulang pasca-sesi dicatat pada kunci `rerun_note` serta `rerun_tests` di `results.json`.
- **HIPOTESIS**: Pada eksekusi *rerun* terukur (`3` detik), total komputasi Rust + Go + Python jauh di bawah anggaran 10 menit; keterlambatan pada sesi asli murni berasal dari latensi penulisan/perapihan kode awal sebelum `LOCK` menit ke-9, bukan dari beban komputasi simulasi.

### Deviasi dari Prompt
- **Jumlah Commit**: Sesi asli `day006c` menggunakan `2` commit, ditambah `1` commit perbaikan pasca-sesi (`fix(day006c): run full pipeline, figures and honest results`) pada sesi `day006d` Bagian A.
- **Ukuran Kode**: Seluruh kode sumber tetap di bawah batas baris (`Rust` 338, `Go` 224, `Python` 339 total termasuk `finalize.py` dan tes).

## 9. Landasan Teori

1. **Perkolasi Bond di Kisi Persegi**: Tiap *bond* terbuka secara independen dengan peluang $p$, dan ambang kritis transisi munculnya klaster tak-hingga bernilai $p_c = 1/2$ eksak (Kesten, 1980) — diimplementasikan pada pembangkit kisi dan simulator modul Rust [`lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs) serta CLI [`main.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/main.rs).
2. **Dualitas-Diri Kisi Persegi**: Kisi dual planar dari persegi panjang simpul $L \times (L+1)$ identik dengan kisi primal setelah rotasi $90^\circ$, sehingga setiap konfigurasi tanpa lintasan primal kiri–kanan berkorespondensi bijektif dengan lintasan penghalang dual atas–bawah dan menghasilkan $P_{\text{cross}}(1/2) = 1/2$ tepat — diimplementasikan pada enumerasi rasional eksak Python di [`percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py) dan diuji di [`lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs).
3. **Algoritma Newman–Ziff Mikrokanonik**: Satu permutasi seragam seluruh *bond* (Fisher–Yates dengan `Xoshiro256Pp`) memberikan sampel keadaan untuk semua jumlah *bond* $n = 0..N$ sekaligus dalam satu penyapuan — diimplementasikan pada `run_newman_ziff` di modul Rust [`lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs).
4. **Struktur Data Union-Find dengan Union by Size dan Path Compression**: Penggabungan pohon berdasarkan ukuran akar dan pelipatan jalur di `find` menekan biaya amortisasi per penyisipan *bond* menjadi hampir konstan $O(\alpha(V))$ — diimplementasikan pada tipe `Uf` di [`lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs).
5. **Konversi Mikrokanonik ke Kanonik Stabil Log-Gamma**: Peluang kanonik $P(p) = \sum_{n=0}^N \binom{N}{n} p^n (1-p)^{N-n} Q(n)$ dievaluasi tanpa *overflow*/*underflow* menggunakan log-gamma $\ln\Gamma$ dan sifat monoton indikator $Q(n) = \mathbb{P}(n^* \le n)$ — diimplementasikan pada `binom_weights` dan `canonical_eval` di [`percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py).
6. **Finite-Size Scaling Jendela Kritis**: Lebar daerah transisi menyempit sebagai $\Delta p(L) \sim L^{-1/\nu}$ dengan eksponen panjang korelasi eksak 2D $\nu = 4/3$, sehingga kecuraman di titik kritis tumbuh menurut hukum pangkat $\mathrm{d}P/\mathrm{d}p(p_c) \sim L^{3/4}$ — diekstraksi pada `ols_loglog` di [`percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py).
7. **Dimensi Fraktal Klaster Kritis**: Pada titik kritis $p_c = 1/2$, klaster perkolasi raksasa bersifat fraktal berongga dengan massa rata-rata $\langle S_{\max} \rangle \sim L^D$ dan dimensi eksak $D = 91/48 \approx 1{,}895833$ (hukum hiperskalar $D = d - \beta/\nu$ dengan $\beta = 5/36, \nu = 4/3$) — diukur lewat `sum_largest` dari [`lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs) dan [`percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py).
8. **Data Collapse Kurva Universal**: Penyekalaan ulang sumbu horizontal menjadi variabel tak-berdimensi $x = (p - p_c)L^{1/\nu}$ meruntuhkan kurva $P(p, L) = \mathcal{F}\bigl((p - p_c)L^{1/\nu}\bigr)$ dari berbagai ukuran kisi ke satu fungsi penskalaan universal $\mathcal{F}(x)$ — diimplementasikan pada `make_all_figures` di [`percolation.py`](../../projects/rust/day006c_percolation_polyglot/py/percolation.py).
9. **Determinisme Paralel Berbasis Potongan Tetap**: Penurunan *seed* `ChunkSeed` per pasangan $(L, \text{chunk})$ lewat `SplitMix64` dan reduksi bilangan bulat `uint64` menurut indeks potongan membuat hasil identik bit-demi-bit antara `1` dan banyak pekerja — diimplementasikan pada orkestrator Go [`go/main.go`](../../projects/rust/day006c_percolation_polyglot/go/main.go).
10. **Inferensi dan Galat Statistik**: Galat standar peluang di $p_c = 1/2$ dihitung dari SE binomial $\sqrt{\hat{P}(1-\hat{P})/M}$ (karena indikator seberang tiap replikasi independen bernilai biner dengan variansi $\le 1/4$), sedangkan galat standar eksponen $\hat{\nu}$ dan $\hat{D}$ diturunkan dari matriks kovarians regresi linier log–log pada potongan independen — diverifikasi pada rangkaian uji [`orchestrator_test.go`](../../projects/rust/day006c_percolation_polyglot/go/orchestrator_test.go).

## 10. Kesimpulan dan Batasan

Kerangka poliglot Rust + Go + Python memisahkan beban kombinatorik $O(N)$ per replikasi ke bahasa sistem bebas-*garbage-collector* (Rust), manajemen proses paralel berbatas waktu ke Go, serta pembuktian rasional eksak dan inferensi penskalaan kritis ke Python. Eksekusi penuh pasca-sesi memulihkan $\hat{\nu} = 1{,}3406 \pm 0{,}0079$ (vs $4/3$) dan $\hat{D} = 1{,}8678 \pm 0{,}0077$ (vs $91/48$) serta memenuhi seluruh kriteria `c1`–`c7`, dengan keterbatasan koreksi ukuran berhingga pada $L \le 128$ dan `c_wall_clock` sesi asli tetap tercatat `belum tercapai` (`789` detik).

## 11. Referensi dan Jejak Bukti

- Kesten, H. (1980). *The critical probability of bond percolation on the square lattice equals 1/2*. Communications in Mathematical Physics, 74(1), 41–59.
- Newman, M. E. J., & Ziff, R. M. (2000). *Efficient Monte Carlo algorithm and high-precision results for percolation*. Physical Review Letters, 85(19), 4104–4107.
- Stauffer, D., & Aharony, A. (1994). *Introduction to Percolation Theory* (2nd ed.). Taylor & Francis.
- Artefak kode terverifikasi: [`projects/rust/day006c_percolation_polyglot/rust/src/lib.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/lib.rs), [`projects/rust/day006c_percolation_polyglot/rust/src/main.rs`](../../projects/rust/day006c_percolation_polyglot/rust/src/main.rs), [`projects/rust/day006c_percolation_polyglot/go/main.go`](../../projects/rust/day006c_percolation_polyglot/go/main.go), dan [`projects/rust/day006c_percolation_polyglot/README.md`](../../projects/rust/day006c_percolation_polyglot/README.md).
