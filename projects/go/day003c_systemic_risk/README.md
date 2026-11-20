# Hari 003c — Risiko Sistemik Jaringan (Go Murni)

Implementasi Go 1.22+ murni (hanya standard library, nol modul eksternal) untuk tiga model penularan risiko sistemik finansial:
1. **Kliring Eisenberg–Noe (2001)**: pembayaran kliring ekuilibrium $p^* = \Phi(p^*) = \min(\bar{p}, e + \Pi^T p^*)$ via iterasi Picard monoton, algoritma default fiktif (solver LU buatan sendiri), dan enumerasi *brute-force* $2^n$ himpunan gagal.
2. **Kaskade Ambang Watts (2002)**: penularan kegagalan pada graf acak Erdős–Rényi $G(n, projects/c/n)$ ($n=100.000$ dibangun dalam $O(n+m)$ via metode lompatan geometrik Batagelj–Brandes), verifikasi rumus pendekatan pohon, batas transisi fase $(c_1, c_2)$, dan distribusi ukuran kaskade benih tunggal.
3. **Biaya Kebangkrutan Rogers–Veraart (2013) & Uji Stres Monte Carlo**: perluasan tak-kontinu dengan pemulihan aset eksternal $\alpha$ dan klaim antarbank $\beta$ saat gagal bayar, uji stres paralel deterministik `20.000` skenario per level guncangan $s \in \{0.2, 0.4, 0.6, 0.8\}$, serta perbandingan guncangan acak vs terarah ke *hub* pada jaringan Poisson vs Pangkat-Hukum.

## Turunan Rumus & Anchor Eksak
### 1. Rantai & Siklus Dua Bank Eisenberg–Noe
- **Rantai $0 \to 1 \to \dots \to n-1 \to \text{sink}$** ($L_{i,i+1} = 1 \implies \bar{p}_i = 1, \Pi_{i,i+1} = 1$):
  - Pola A ($e_i = a \in (0, 1/n]$ konstan): dengan induksi, $p_0^* = a$ dan $p_i^* = \min(1, a + p_{i-1}^*) = (i+1)a$.
  - Pola B ($e_0 = a_0, e_i = a$ untuk $i \ge 1$): selama belum mencapai 1, $p_i^*$ tumbuh aritmetik lalu terpotong pada 1: $p_i^* = \min(1, a_0 + i \cdot a)$.
- **Siklus dua bank saling berutang ($0 \leftrightarrow 1$)** dengan proporsi antarbank $\Pi_{01} = \gamma_0 \in (0,1), \Pi_{10} = \gamma_1 \in (0,1)$ dan sisanya ke sink eksternal ($\bar{p}_0 = \bar{p}_1 = 1$):
  Bila $e_0 + \gamma_1 < 1$ dan $e_1 + \gamma_0 < 1$, kedua bank gagal sehingga memenuhi sistem linear $p_0 = e_0 + \gamma_1 p_1$ dan $p_1 = e_1 + \gamma_0 p_0$, yang memiliki solusi tertutup unik:
  $$p_0^* = \frac{e_0 + \gamma_1 e_1}{1 - \gamma_0 \gamma_1}, \qquad p_1^* = \frac{e_1 + \gamma_0 e_0}{1 - \gamma_0 \gamma_1}.$$

### 2. Koreksi Rumus Pendekatan Pohon Watts (2002)
Dalam pendekatan pohon lokal $G(n, projects/c/n)$, peluang $q$ bahwa tetangga $v$ berderajat $k$ yang dicapai melalui satu sisi dari simpul $u$ (yang **belum** gagal) mengalami kegagalan mensyaratkan sedikitnya $m \ge \lceil \phi k \rceil$ dari $k-1$ tetangga $v$ lainnya telah gagal, yaitu indikator $\mathbf{1}\{m/k \ge \phi\}$ (Watts 2002 Eq. 3–4; Gleeson & Cahalane 2007 Eq. 1–2). Penulisan $\mathbf{1}\{(m+1)/k \ge \phi\}$ pada persamaan $q$ bersifat sirkular (mengandaikan $u$ sudah gagal agar $v$ gagal untuk menginfeksi $u$) dan menyebabkan simpul rentan gagal spontan pada $m=0$.

### 3. Model Rogers–Veraart (2013) & Rantai Tertutup
Ketika bank $i$ memiliki total nilai aset sebelum biaya $A_i(p) = e_i + \sum_j \Pi_{ji} p_j < \bar{p}_i$, likuidasi kebangkrutan menyisakan porsi $\alpha \in [0,1]$ dari aset eksternal dan $\beta \in [0,1]$ dari klaim antarbank:
$$\Phi_i^{\mathrm{RV}}(p) = \begin{cases} \bar{p}_i, & \text{jika } e_i + \sum_j \Pi_{ji} p_j \ge \bar{p}_i, \\ \alpha e_i + \beta \sum_j \Pi_{ji} p_j, & \text{jika } e_i + \sum_j \Pi_{ji} p_j < \bar{p}_i. \end{cases}$$
Pada rantai $0 \to 1 \to \dots \to n-1$ ($\bar{p}_i=1$) di mana setiap bank gagal ($e_i + p_{i-1}^* < 1$), rekursi linear $p_0^* = \alpha e_0$, $p_i^* = \alpha e_i + \beta p_{i-1}^*$ menghasilkan solusi tertutup $p_i^* = \alpha \sum_{k=0}^i \beta^{i-k} e_k$.

## Cara Build, Test, dan Run
```bash
cd projects/go/day003c_systemic_risk
gofmt -l .
go vet ./...
go test -race ./...
go run ./cmd/run
```

## Struktur Paket
- `internal/rng` — PRNG `xoshiro256++` + SplitMix64 + Box–Muller + Poisson Knuth.
- `internal/parallel` — worker pool goroutine deterministik berbasis chunk tetap (`RunChunks`).
- `internal/linalg` — dekomposisi LU dengan *partial pivoting* (`SolveLU`).
- `internal/graphs` — jaringan antarbank berarah (Poisson & Pangkat-Hukum) dan graf CSR tak berarah $G(n, projects/c/n)$ Batagelj–Brandes $O(n+m)$.
- `internal/clearing` — kliring Eisenberg–Noe (`Picard`, `PicardFromZero`, `FictitiousDefault`, `BruteForce`).
- `internal/cascade` — simulasi kaskade Watts (FIFO & LIFO), `MeanField`, dan `PhaseBoundaries`.
- `internal/rv` — kliring Rogers–Veraart (`Clear`, `ChainClosedForm`) dan `RunMonteCarloStress`.
- `internal/stats`, `internal/jsonout`, `internal/svg` — statistik VaR/ES, JSON 17 digit, dan grafik SVG.
- `cmd/run/main.go` — runner utama penghasil [../../../outputs/day003c/results.json](../../../outputs/day003c/results.json) dan 5 grafik SVG.
