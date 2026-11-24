# Penjelasan Proyek `day006b`: Jalan Menghindari-Diri (SAW) Poliglot di Kisi Persegi — Enumerasi Eksak, Pivot Monte Carlo, dan Eksponen Kritis

## 1. Identitas dan Ringkasan Proyek

- **ID Sesi**: `day006b` (Hari 6, sesi `b`, topik `bonus` di luar `ROADMAP.md`)
- **Tanggal**: `2026-11-24` (WITA)
- **Topik**: Jalan menghindari-diri (*self-avoiding walk*, SAW) pada kisi persegi $\mathbb{Z}^2$
- **Jenis**: `sim` (simulasi fisika statistik & kombinatorika komputasional)
- **Bahasa**: `python+rust+go` (poliglot: inti komputasi Rust murni, orkestrator konkuren Go pustaka standar, serta referensi *brute-force*, regresi statistik, dan visualisasi Python)
- **Kebijakan Folder**: Skrip [`scripts/update_readme.py`](../../scripts/update_readme.py) hanya memindai `projects/{rust,go,typescript,ts,cpp,c,ruby}/*/meta.json` dan tidak memindai `projects/python/`. Agar `python scripts/update_readme.py --check` lulus tanpa memodifikasi skrip pemeriksa, folder induk ditempatkan di [`projects/rust/day006b_saw_polyglot/`](../../projects/rust/day006b_saw_polyglot/README.md) dengan `language: "python+rust+go"` dan `type: "sim"` pada [`meta.json`](../../projects/rust/day006b_saw_polyglot/meta.json). Sesi `b` dipakai sesuai aturan karena folder `day006a` dari percobaan sebelumnya sudah ada dan tidak boleh dihapus.

## 2. Masalah dan Tujuan

Jalan menghindari-diri sepanjang $n$ langkah di $\mathbb{Z}^2$ adalah lintasan kisi mulai dari titik asal $(0,0)$ yang tidak pernah mengunjungi satu titik kisi lebih dari sekali. Proyek ini menjawab tiga pertanyaan kuantitatif dalam batas disiplin waktu dinding ketat 10 menit:
1. Berapa jumlah eksak lintasan $c_n$ hingga $N_{\max} = 18$, dan apakah reduksi simetri $\times 8$ di Rust serta pemecahan awalan paralel di Go menghasilkan bilangan bulat yang persis sama dengan *brute-force* murni Python?
2. Seberapa dekat taksiran konstanta konektivitas $\hat{\mu}$ dan eksponen entropik $\hat{\gamma}$ dari metode rasio $r_n = c_n / c_{n-1}$ terhadap konstanta acuan kisi persegi ($\mu \approx 2{,}63816$, $\gamma = 43/32 = 1{,}34375$)?
3. Apakah pencuplikan MCMC algoritma pivot Madras–Sokal memulihkan eksponen metrik Flory 2D eksak $\nu = 3/4 = 0{,}75$ melalui penskalaan jarak kuadrat ujung-ke-ujung $\langle R^2(N) \rangle \sim A N^{2\nu}$ pada $N \in \{50, 100, 200, 400\}$?

## 3. Apa yang Dikerjakan Proyek

 Sistem membangun pipa poliglot tiga tahap yang diorkestrasi oleh [`run_all.sh`](../../projects/rust/day006b_saw_polyglot/run_all.sh). Tahap pertama mengompilasi biner Rust `saw`. Tahap kedua menjalankan orkestrator Go `orchestrate` yang memanggil `saw count` satu-proses (sekaligus mengakumulasi peta panas titik ujung $n = 14$), membagi pohon awalan panjang $K = 7$ menjadi $543$ tugas ke kolam pekerja Rust persisten untuk `1` dan `N` *worker*, serta mengeksekusi simulasi `saw pivot` pada $N \in \{50, 100, 200, 400\}$ dan $N = 300$ ($100$ sampel lintasan). Tahap ketiga di Python memverifikasi $c_n$ ($n \le 12$) dan ketiadaan tabrakan pada lintasan contoh, menjalankan regresi kuadrat terkecil biasa (*ordinary least squares*, OLS), menyimpan [`outputs/day006b/results.json`](../../outputs/day006b/results.json), dan merender 4 pasang grafik `PNG`/`SVG`.

## 4. Uraian Kode

- [`projects/rust/day006b_saw_polyglot/rust/src/lib.rs`](../../projects/rust/day006b_saw_polyglot/rust/src/lib.rs): Memuat PRNG `Xoshiro256Pp` (diinisialisasi lewat pencampur `SplitMix64` dan pencuplikan batas tak-bias Lemire `bounded`), papan bit `Bitmap([u64; 64])` berukuran $64 \times 64$, transformasi grup dihedral `apply_d4`, enumerator simetri `enumerate_symmetry` dengan jalan pintas penghitungan daun (*leaf-counting*), enumerator awalan `enumerate_from_prefix`, tabel hash bertanda bebas-alokasi `StampSet`, dan mesin MCMC `run_pivot` beserta 4 tes unit.
- [`projects/rust/day006b_saw_polyglot/rust/src/main.rs`](../../projects/rust/day006b_saw_polyglot/rust/src/main.rs): Antarmuka baris perintah `saw count` (mode langsung maupun pekerja pipa `--stdin` dengan `BufWriter::flush` tiap baris JSON) dan `saw pivot`.
- [`projects/rust/day006b_saw_polyglot/go/main.go`](../../projects/rust/day006b_saw_polyglot/go/main.go) & [`orchestrator_test.go`](../../projects/rust/day006b_saw_polyglot/go/orchestrator_test.go): Fungsi `GenPrefixes` membangkitkan $c_7 / 4 = 543$ awalan langkah pertama Timur berbobot $4$, `RunSharded` mendistribusikan tugas lewat kanal Go ke proses `saw count --stdin` persisten di bawah `context.WithTimeout` lalu menjumlahkan `uint64` menurut urutan indeks tugas `0..542`, dan `DeriveSeed` menurunkan *seed* deterministik per $N$ dari `20261124`. Dilengkapi 3 tes unit Go.
- [`projects/rust/day006b_saw_polyglot/py/saw_analysis.py`](../../projects/rust/day006b_saw_polyglot/py/saw_analysis.py) & [`test_saw.py`](../../projects/rust/day006b_saw_polyglot/py/test_saw.py): Fungsi `brute_force_saw` (DFS 4 arah murni tanpa simetri), `validate_saw`, `ols_fit`, `fit_ratio` (metode delta untuk SE $\hat{\gamma}$), `fit_scaling`, `make_figures`, dan `check_figures` (pemeriksa XML `ElementTree`), beserta 4 tes `unittest`.

## 5. Metode dan Alur Kerja

1. **Reduksi Simetri $\times 8$**: Langkah pertama dipaksa ke Timur $(1,0)$ sehingga menyumbang bobot $4$. Selama lintasan tetap lurus di sumbu $+x$, lintasan invarian terhadap pencerminan $y \mapsto -y$ (bobot tetap $4$). Saat belokan ortogonal pertama terjadi, arah dipaksa ke Utara $(+y)$ sehingga mencerminkan pasangan Utara/Selatan dan bobot seluruh perpanjangan berikutnya menjadi $4 \times 2 = 8$.
2. **Sharding Awalan Deterministik**: Go mengekspansi pohon hingga kedalaman $K = 7$, mencatat $c_1..c_6$ langsung, dan mengirim ke-543 daun sebagai JSON ke pekerja Rust. Penjumlahan dilakukan dalam urutan `id = 0..542` murni memakai `uint64`.
3. **Pivot Madras–Sokal & Regresi**: Pada tiap langkah MCMC dipilih simpul poros $k \in \{0,\dots,N-1\}$ dan simetri non-identitas $g \in D_4 \setminus \{I\}$ secara seragam; ekor $i > k$ diputar/dicerminkan $\omega'_i = \omega_k + g(\omega_i - \omega_k)$ dan diterima bila bebas tabrakan terhadap prefiks $0..k-1$ pada `StampSet`. Galat standar $\langle R^2 \rangle$ dihitung dari $50$ *batch* beranggotakan $4.000$ langkah pasca-*burn-in*. Python kemudian mencocokkan $r_n = \beta_0 + \beta_1 (1/n)$ pada $n \in [11, 18]$ ($\hat{\mu} = \beta_0$, $\hat{\gamma} = 1 + \beta_1 / \beta_0$) serta $\ln \langle R^2 \rangle = \ln A + 2\hat{\nu} \ln N$.

## 6. Output yang Ditargetkan

- Kecocokan bilangan bulat persis $c_n$ antara Python *brute-force* ($n \le 12$), Rust satu-proses, dan Go *sharded* ($n \le 18$), serta pembandingan jujur terhadap deret ingatan OEIS A001411.
- Kriteria eksponen *a priori*: $|\hat{\mu} - 2{,}638| < 0{,}03$, $|\hat{\gamma} - 1{,}34375| < 0{,}08$, dan $|\hat{\nu} - 0{,}75| < 0{,}03$.
- Validitas `100` contoh lintasan $N = 300$ (`0` tabrakan), invarian pekerja `1` vs `10`, empat grafik `PNG`+`SVG` sah (`day006b_saw_samples`, `day006b_saw_scaling`, `day006b_saw_ratio`, `day006b_saw_endpoint_heatmap`), dan penyelesaian seluruh sesi $\le 600$ detik dari `T0` skrip `dl.sh`.

## 7. Status Pencapaian

Seluruh nilai di bawah ini diambil langsung dari [`outputs/day006b/results.json`](../../outputs/day006b/results.json):

| Kriteria | Status | Nilai Terukur | Ambang A Priori & Bukti |
|---|---|---|---|
| `c1_bruteforce_matches_rust_n12` | `tercapai` | `c_12 = 324932` | Persis sama Python = Rust = Go pada `n = 1..12` |
| `c2_go_matches_rust_single_all_n` | `tercapai` | `c_18 = 124658732` | Go *sharded* ($543$ tugas) == Rust 1-proses pada `n = 1..18` |
| `c3_memory_sequence_a001411` | `tercapai` | `0` selisih | Cocok penuh terhadap 18 suku acuan ingatan A001411 |
| `c4_mu_hat` | `tercapai` | `2,628998 ± 0,017604` | `\|2,628998 - 2,638\| = 0,009002 < 0,03` |
| `c5_gamma_hat` | `tercapai` | `1,383377 ± 0,095926` | `\|1,383377 - 1,34375\| = 0,039627 < 0,08` |
| `c6_nu_hat` | `tercapai` | `0,744684 ± 0,001017` | `\|0,744684 - 0,75\| = 0,005316 < 0,03` |
| `c7_sample_walks_valid` | `tercapai` | `100 / 100` lintasan sah | `0` tabrakan pada seluruh sampel $N = 300$ (`validate_saw`) |
| `c8_figures_valid` | `tercapai` | `4` pasang PNG+SVG | Berkas ada, tak kosong, seluruh SVG lolos `ElementTree.parse` |
| `c9_worker_invariance` | `tercapai` | `1` vs `10` *worker* | Vektor `counts[0..18]` identik bit-demi-bit |
| `c10_ten_minute_wall_clock` | `tercapai` | `elapsed <= 600` detik | Diukur langsung dari `T0 = 1795506051` (`dl.sh`) |

## 8. Jika Belum Tercapai atau Belum Diketahui

- **FAKTA**: Seluruh item MUST dan SHOULD (`c1`–`c10`) pada sesi `day006b` telah dijalankan, diverifikasi otomatis oleh tes unit ketiga bahasa, dan berstatus `tercapai`. Item tingkat COULD (seperti penampil HTML interaktif atau peningkatan $N_{\max} > 18$) sengaja tidak dikerjakan pada sesi `day006b` karena spesifikasi tugas menetapkan bahwa COULD tidak dikerjakan demi disiplin batas 10 menit.
- **FAKTA Deviasi Kecil Asimtotik**: Taksiran $\hat{\nu} = 0{,}74468 \pm 0{,}00102$ sedikit di bawah nilai eksak Flory $0{,}75$, dan $\hat{\gamma} = 1{,}38338 \pm 0{,}09593$ memiliki galat standar regresi $\approx 0{,}096$ akibat osilasi paritas genap-ganjil pada $r_n$ kisi persegi bipartit.
- **HIPOTESIS**: Pergeseran kecil $\hat{\nu} - 0{,}75 \approx -0{,}0053$ berasal dari suku koreksi penskalaan berhingga (*corrections-to-scaling*) $\langle R^2(N) \rangle = A N^{2\nu}(1 + B N^{-\Delta_1} + \dots)$ dengan $\Delta_1 \approx 3/2$ yang masih terasa pada $N = 50..100$; pembuktiannya memerlukan pencuplikan $N \in \{800, 1600, 3200\}$ dan fit non-linear dua suku.

## 9. Landasan Teori

1. **Jalan Menghindari-Diri (Self-Avoiding Walk, SAW)**: Lintasan kombinatorial pada kisi $\mathbb{Z}^d$ tanpa perpotongan diri yang memodelkan volume terkecualikan (*excluded volume*) rantai polimer dalam pelarut baik, diimplementasikan pada `Bitmap` dan `validate_saw`.
2. **Konstanta Konektivitas ($\mu$) dan Subaditivitas Hammersley–Morton**: Karena $c_{n+m} \le c_n c_m$, barisan $\ln c_n$ bersifat subaditif sehingga menurut teorema Fekete limit $\mu = \lim_{n\to\infty} c_n^{1/n}$ ada dan berhingga, ditaksir lewat intersep `fit_ratio`.
3. **Asimtotik Jumlah Konformasi dan Eksponen Entropik ($\gamma$)**: Hukum penskalaan $c_n \sim \mathcal{A}\,\mu^n n^{\gamma - 1}$ menghubungkan pertumbuhan jumlah konfigurasi dengan eksponen kritis universal $\gamma = 43/32$ di 2D (Nienhuis), diekstraksi dari kemiringan rasio berurutan $r_n = c_n/c_{n-1} \approx \mu + \mu(\gamma - 1)/n$.
4. **Eksponen Metrik Flory ($\nu$)**: Jarak kuadrat ujung-ke-ujung memenuhi $\langle R^2(N) \rangle \sim A N^{2\nu}$ dengan $\nu = 3/4$ eksak di $d=2$ (jauh lebih mengembang dibandingkan gerak Brown biasa $\nu = 1/2$), diestimasi pada `fit_scaling` dan diplot pada [`day006b_saw_scaling.svg`](../../outputs/figures/day006b_saw_scaling.svg).
5. **Grup Dihedral $D_4$ pada Kisi Persegi**: Delapan isometri kisi (4 rotasi dan 4 refleksi) melestarikan ketetanggaan dan sifat bebas-tabrakan, dipakai di `apply_d4` untuk memangkas pohon DFS sebesar faktor $\times 8$ serta menyimetrisasikan peta panas [`day006b_saw_endpoint_heatmap.svg`](../../outputs/figures/day006b_saw_endpoint_heatmap.svg).
6. **Dekomposisi Awalan (*Prefix Sharding*) Pohon DFS**: Setiap lintasan panjang $n \ge K$ memiliki awalan tunggal panjang $K$, sehingga partisi ruang pencarian pada kedalaman $K = 7$ di `GenPrefixes` dan `enumerate_from_prefix` bersifat saling asing dan lengkap tanpa komunikasi antar-pekerja.
7. **Algoritma Pivot Madras–Sokal**: Operator MCMC global yang memilih simpul poros acak dan menerapkan transformasi grup kisi pada sub-rantai ekor; sifat simetris $g^{-1} \in D_4 \setminus \{I\}$ menjamin keseimbangan terperinci (*detailed balance*) terhadap distribusi seragam seluruh SAW panjang $N$ di `pivot_step`.
8. **Tabel Hash Bertanda (*Stamped Hash Table*)**: Struktur `StampSet` di Rust menggunakan pencacah generasi `stamp: u32` sehingga pengosongan tabel di tiap usulan pivot berbiaya $O(1)$ tanpa pengisian ulang memori `O(cap)`.
9. **Pembangkit Bilangan Acak Semu `xoshiro256++` dan Lemire Bounded**: Transformasi pergeseran-rotasi 256-bit Blackman–Vigna pada `Xoshiro256Pp` yang dipadukan dengan perkalian 128-bit Lemire menghasilkan indeks poros dan simetri seragam tanpa bias modulo.
10. **Metode *Batch Means* dan Metode Delta Kovarians OLS**: Pengelompokan rantai Markov menjadi $50$ blok berurutan di `run_pivot` menyerap autokorelasi MCMC ke dalam variansi rata-rata blok, sedangkan perambatan kovarians $(\beta_0, \beta_1)$ lewat gradien $\nabla g = (-\beta_1/\beta_0^2, 1/\beta_0)$ di `fit_ratio` memberi galat standar $\text{SE}(\hat{\gamma})$ yang konsisten.

## 10. Kesimpulan dan Batasan

Arsitektur poliglot Rust + Go + Python berhasil memverifikasi silang enumerasi eksak SAW hingga $N_{\max} = 18$ ($c_{18} = 124.658.732$) dengan kecocokan bulat sempurna antar-bahasa dan antar-jumlah pekerja, sekaligus menaksir ketiga konstanta kritis ($\hat{\mu} = 2{,}6290$, $\hat{\gamma} = 1{,}3834$, $\hat{\nu} = 0{,}7447$) di dalam ambang *a priori* serta menyelesaikan seluruh alur di bawah batas keras 10 menit. Batasan utama eksperimen ini adalah efek ukuran berhingga dan osilasi paritas pada $N_{\max} = 18$; perpanjangan alami berikutnya adalah pemisahan regresi rasio genap/ganjil (ekstrapolasi Richardson orde kedua) serta struktur data pohon biner *SAW-tree* Kennedy agar uji tabrakan pivot turun ke $O(\log N)$ sub-linear.

## 11. Referensi dan Jejak Bukti

- Madras, N., & Sokal, A. D. (1988). *The pivot algorithm: A highly efficient Monte Carlo method for the self-avoiding walk*. Journal of Statistical Physics, 50(1–2), 109–186.
- Nienhuis, B. (1982). *Exact critical point and critical exponents of $O(n)$ models in two dimensions*. Physical Review Letters, 49(15), 1062–1065.
- Blackman, D., & Vigna, S. (2021). *Scrambled linear pseudorandom number generators*. ACM TOMS, 47(4), 1–32.
- Artefak lokal terverifikasi: [`projects/rust/day006b_saw_polyglot/README.md`](../../projects/rust/day006b_saw_polyglot/README.md), [`outputs/day006b/results.json`](../../outputs/day006b/results.json), [`day006b_saw_samples.svg`](../../outputs/figures/day006b_saw_samples.svg), [`day006b_saw_scaling.svg`](../../outputs/figures/day006b_saw_scaling.svg), [`day006b_saw_ratio.svg`](../../outputs/figures/day006b_saw_ratio.svg), dan [`day006b_saw_endpoint_heatmap.svg`](../../outputs/figures/day006b_saw_endpoint_heatmap.svg).
