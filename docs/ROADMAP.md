# 🗺️ ROADMAP 100QB

Daftar 100 topik, berurutan dari dasar ke lanjut, selang-seling **quant** (keuangan/kuantitatif)
dan **sim** (simulasi/perhitungan iseng). Topik yang sudah selesai dicentang `[x]`.

Format: `- [ ] **NNN** · `jenis` · **Topik** — tujuan.` (dibaca otomatis oleh `scripts/new_day.py`).

## Fase 1 — Fondasi (Hari 001–025)

- [x] **001** · `quant` · **Return harian & kumulatif** — Menghitung simple return dan return kumulatif dari deret harga.
- [x] **002** · `sim` · **Estimasi π dengan Monte Carlo** — Memperkirakan π dari titik acak di dalam persegi dan lingkaran.
- [ ] **003** · `quant` · **Log return vs simple return** — Memahami kapan memakai log return dan sifat aditifnya.
- [ ] **004** · `sim` · **Random walk 1D** — Mensimulasikan jalan acak dan melihat sebaran posisi akhirnya.
- [ ] **005** · `quant` · **Volatilitas & annualisasi** — Mengukur risiko dengan standar deviasi dan skala akar waktu.
- [ ] **006** · `sim` · **Lempar koin & hukum bilangan besar** — Melihat rata-rata sampel konvergen ke peluang sebenarnya.
- [ ] **007** · `quant` · **Ambil data harga dengan yfinance** — Mengunduh, menyimpan cache, dan membersihkan data harga saham.
- [ ] **008** · `sim` · **Lempar dadu & teorema limit pusat** — Melihat jumlah dadu mendekati distribusi normal.
- [ ] **009** · `quant` · **Distribusi return & normalitas** — Membandingkan histogram return dengan kurva normal (skew, kurtosis).
- [ ] **010** · `sim` · **Gerak peluru** — Mensimulasikan lintasan proyektil dengan dan tanpa hambatan udara.
- [ ] **011** · `quant` · **Sharpe ratio** — Mengukur return per unit risiko dan membandingkan beberapa aset.
- [ ] **012** · `sim` · **Peluruhan radioaktif** — Mensimulasikan peluruhan atom acak dan membandingkan dengan rumus eksponensial.
- [ ] **013** · `quant` · **Maximum drawdown** — Menghitung penurunan terdalam dari puncak dan durasi pemulihannya.
- [ ] **014** · `sim` · **Game of Life Conway** — Membuat automata seluler sederhana dan mengamati pola yang muncul.
- [ ] **015** · `quant` · **Moving average & crossover** — Menghitung SMA/EMA dan menandai sinyal golden/death cross.
- [ ] **016** · `sim` · **Fraktal Mandelbrot** — Menggambar himpunan Mandelbrot dengan iterasi bilangan kompleks.
- [ ] **017** · `quant` · **Bunga majemuk & nilai waktu uang** — Menghitung present/future value dan efek frekuensi compounding.
- [ ] **018** · `sim` · **Birthday paradox** — Mensimulasikan peluang dua orang berulang tahun sama.
- [ ] **019** · `quant` · **Korelasi antar saham** — Membuat matriks korelasi return dan heatmap-nya.
- [ ] **020** · `sim` · **Monty Hall** — Membuktikan lewat simulasi bahwa berganti pintu lebih menguntungkan.
- [ ] **021** · `quant` · **Rolling volatility** — Melihat volatilitas yang berubah seiring waktu dengan jendela bergulir.
- [ ] **022** · `sim` · **Jarum Buffon** — Mengestimasi π dari jarum yang dijatuhkan secara acak.
- [ ] **023** · `quant` · **Beta & CAPM** — Mengukur sensitivitas saham terhadap pasar dengan regresi.
- [ ] **024** · `sim` · **Random walk 2D** — Mensimulasikan gerak acak di bidang dan jarak rata-rata dari titik awal.
- [ ] **025** · `quant` · **Portofolio dua aset** — Menghitung return dan risiko portofolio dari bobot dan korelasi.

## Fase 2 — Menengah (Hari 026–050)

- [ ] **026** · `sim` · **Antrean kasir (M/M/1)** — Mensimulasikan kedatangan dan layanan pelanggan lalu mengukur waktu tunggu.
- [ ] **027** · `quant` · **Gerak Brown geometrik (GBM)** — Mensimulasikan jalur harga saham dengan model GBM.
- [ ] **028** · `sim` · **Pendulum sederhana** — Menyelesaikan persamaan gerak pendulum secara numerik.
- [ ] **029** · `quant` · **Monte Carlo harga saham** — Memproyeksikan sebaran harga masa depan dari banyak jalur GBM.
- [ ] **030** · `sim` · **Metode Euler vs Runge-Kutta** — Membandingkan akurasi dua metode integrasi numerik ODE.
- [ ] **031** · `quant` · **Value at Risk historis** — Menghitung VaR dari distribusi return masa lalu.
- [ ] **032** · `sim` · **Model predator-mangsa Lotka-Volterra** — Mensimulasikan dinamika populasi kelinci dan rubah.
- [ ] **033** · `quant` · **VaR parametrik & Monte Carlo** — Membandingkan tiga metode VaR dan asumsi masing-masing.
- [ ] **034** · `sim` · **Model epidemi SIR** — Mensimulasikan penyebaran penyakit dan efek angka reproduksi.
- [ ] **035** · `quant` · **Expected Shortfall (CVaR)** — Mengukur rata-rata kerugian di ekor distribusi.
- [ ] **036** · `sim` · **Segitiga Sierpinski (chaos game)** — Membentuk fraktal dari titik acak dengan aturan sederhana.
- [ ] **037** · `quant` · **Efficient frontier** — Membuat ribuan portofolio acak dan memetakan frontier efisien.
- [ ] **038** · `sim` · **Integrasi Monte Carlo** — Menghitung integral tentu dengan sampel acak dan membandingkan error-nya.
- [ ] **039** · `quant` · **Portofolio minimum variance** — Mencari bobot dengan risiko terendah menggunakan optimisasi.
- [ ] **040** · `sim` · **Difusi panas 1D** — Mensimulasikan perambatan panas pada batang dengan beda hingga.
- [ ] **041** · `quant` · **Portofolio Sharpe maksimum** — Mencari portofolio tangensial dengan scipy.optimize.
- [ ] **042** · `sim` · **Automata seluler Rule 30 & 110** — Menjelajahi automata 1D Wolfram dan pola kompleksnya.
- [ ] **043** · `quant` · **Backtest strategi SMA crossover** — Menguji strategi sederhana dan membandingkan dengan buy-and-hold.
- [ ] **044** · `sim` · **Himpunan Julia** — Menggambar fraktal Julia untuk berbagai konstanta kompleks.
- [ ] **045** · `quant` · **Biaya transaksi & slippage** — Melihat dampak biaya terhadap hasil backtest.
- [ ] **046** · `sim` · **Penyebaran kebakaran hutan** — Mensimulasikan api menjalar di grid dengan peluang tertentu.
- [ ] **047** · `quant` · **Strategi momentum** — Menguji strategi membeli aset dengan performa terbaik periode lalu.
- [ ] **048** · `sim` · **Perkolasi** — Mencari ambang peluang saat jalur terbuka menembus grid.
- [ ] **049** · `quant` · **Mean reversion & z-score** — Membuat sinyal beli/jual dari penyimpangan terhadap rata-rata.
- [ ] **050** · `sim` · **Gambler's ruin** — Mensimulasikan peluang bangkrut penjudi dengan modal terbatas.

## Fase 3 — Lanjut (Hari 051–075)

- [ ] **051** · `quant` · **Model Black-Scholes** — Menghitung harga opsi call/put Eropa dengan rumus tertutup.
- [ ] **052** · `sim` · **Simulasi lalu lintas Nagel-Schreckenberg** — Melihat kemacetan muncul dari aturan mobil sederhana.
- [ ] **053** · `quant` · **Greeks opsi** — Menghitung delta, gamma, vega, theta, rho dan memvisualisasikannya.
- [ ] **054** · `sim` · **Model Ising 2D** — Mensimulasikan magnetisasi dengan algoritma Metropolis.
- [ ] **055** · `quant` · **Harga opsi dengan Monte Carlo** — Membandingkan harga Monte Carlo dengan Black-Scholes.
- [ ] **056** · `sim` · **Simulated annealing (TSP)** — Mencari rute pendek untuk travelling salesman dengan annealing.
- [ ] **057** · `quant` · **Implied volatility** — Mencari volatilitas tersirat dari harga opsi dengan root finding.
- [ ] **058** · `sim` · **Boids (perilaku kawanan)** — Mensimulasikan kawanan burung dengan tiga aturan lokal.
- [ ] **059** · `quant` · **Model binomial pohon opsi** — Menilai opsi Eropa dan Amerika dengan pohon binomial.
- [ ] **060** · `sim` · **Lorenz attractor** — Menyelesaikan sistem Lorenz dan melihat sensitivitas kondisi awal.
- [ ] **061** · `quant` · **Volatility smile** — Memvisualisasikan implied volatility terhadap strike.
- [ ] **062** · `sim` · **Peta logistik & bifurkasi** — Menggambar diagram bifurkasi menuju chaos.
- [ ] **063** · `quant` · **Model GARCH(1,1) sederhana** — Mengestimasi volatilitas yang bergerombol (volatility clustering).
- [ ] **064** · `sim` · **Pola Turing reaksi-difusi** — Membentuk pola mirip kulit hewan dari model Gray-Scott.
- [ ] **065** · `quant` · **Proses Ornstein-Uhlenbeck** — Mensimulasikan proses mean-reverting untuk suku bunga/spread.
- [ ] **066** · `sim` · **Algoritma genetika sederhana** — Mengevolusikan solusi untuk masalah optimisasi sederhana.
- [ ] **067** · `quant` · **Model suku bunga Vasicek** — Mensimulasikan jalur suku bunga dan harga obligasi zero-coupon.
- [ ] **068** · `sim` · **Diffusion-limited aggregation** — Menumbuhkan struktur fraktal dari partikel yang berdifusi.
- [ ] **069** · `quant` · **Harga & durasi obligasi** — Menghitung harga, yield, durasi, dan konveksitas obligasi.
- [ ] **070** · `sim` · **Rantai Markov cuaca** — Mensimulasikan cuaca dengan matriks transisi dan distribusi stasioner.
- [ ] **071** · `quant` · **Pairs trading & kointegrasi** — Menguji kointegrasi dua saham dan membuat strategi spread.
- [ ] **072** · `sim` · **N-body gravitasi sederhana** — Mensimulasikan orbit beberapa benda dengan hukum Newton.
- [ ] **073** · `quant` · **Kriteria Kelly** — Menentukan ukuran taruhan/posisi optimal dan risikonya.
- [ ] **074** · `sim` · **Model segregasi Schelling** — Melihat segregasi muncul dari preferensi tetangga yang ringan.
- [ ] **075** · `quant` · **Jump diffusion Merton** — Menambahkan lompatan harga ke GBM dan melihat ekor tebal.

## Fase 4 — Mahir & Penutup (Hari 076–100)

- [ ] **076** · `sim` · **Gerak Brown partikel** — Mensimulasikan partikel yang ditumbuk molekul dan mean squared displacement.
- [ ] **077** · `quant` · **Model Heston (stochastic volatility)** — Mensimulasikan harga dengan volatilitas yang ikut acak.
- [ ] **078** · `sim` · **Gelombang 1D** — Menyelesaikan persamaan gelombang tali bergetar secara numerik.
- [ ] **079** · `quant` · **Opsi Asia & barrier** — Menilai opsi eksotis yang bergantung pada jalur dengan Monte Carlo.
- [ ] **080** · `sim` · **Pohon fraktal L-system** — Menggambar tanaman fraktal dari aturan tata bahasa sederhana.
- [ ] **081** · `quant` · **Variance reduction** — Mempercepat Monte Carlo dengan antithetic dan control variates.
- [ ] **082** · `sim` · **Ekosistem kebun (agent-based)** — Mensimulasikan tanaman, hama, dan pemangsa di satu kebun.
- [ ] **083** · `quant` · **Risk parity** — Menyusun portofolio dengan kontribusi risiko yang seimbang.
- [ ] **084** · `sim` · **Simulasi inventori toko** — Mencari kebijakan stok optimal dengan permintaan acak.
- [ ] **085** · `quant` · **Black-Litterman** — Menggabungkan ekuilibrium pasar dengan pandangan investor.
- [ ] **086** · `sim` · **Pasir Abelian (sandpile)** — Mengamati self-organized criticality dari tumpukan pasir.
- [ ] **087** · `quant` · **PCA pada return saham** — Menemukan faktor utama penggerak sekumpulan saham.
- [ ] **088** · `sim` · **Langton's ant** — Mensimulasikan semut sederhana yang akhirnya membangun "jalan raya".
- [ ] **089** · `quant` · **Model faktor Fama-French** — Menjelaskan return dengan faktor pasar, size, dan value.
- [ ] **090** · `sim` · **MCMC Metropolis-Hastings** — Mengambil sampel dari distribusi sulit dengan rantai Markov.
- [ ] **091** · `quant` · **Walk-forward backtest** — Menghindari overfitting dengan validasi bergulir.
- [ ] **092** · `sim` · **Bootstrap statistik** — Mengestimasi interval kepercayaan dengan resampling.
- [ ] **093** · `quant` · **Monte Carlo pensiun** — Mensimulasikan peluang dana pensiun cukup sampai akhir.
- [ ] **094** · `sim` · **Optimasi kawanan partikel (PSO)** — Mencari minimum fungsi dengan kawanan partikel.
- [ ] **095** · `quant` · **Stress test portofolio** — Menguji portofolio terhadap skenario krisis historis dan hipotetis.
- [ ] **096** · `sim` · **Hujan & genangan (cellular water flow)** — Mensimulasikan aliran air di permukaan tanah grid.
- [ ] **097** · `quant` · **Regime switching (Markov)** — Mendeteksi rezim pasar bull/bear dengan model dua keadaan.
- [ ] **098** · `sim` · **Kebun Simulasi: musim tanam** — Menggabungkan cuaca acak, pertumbuhan, dan panen dalam satu simulasi.
- [ ] **099** · `quant` · **Dashboard risiko portofolio** — Merangkum return, volatilitas, Sharpe, VaR, dan drawdown dalam satu laporan.
- [ ] **100** · `sim` · **Retrospektif 100 hari** — Memvisualisasikan perjalanan 100 notebook dan pelajaran utamanya.
