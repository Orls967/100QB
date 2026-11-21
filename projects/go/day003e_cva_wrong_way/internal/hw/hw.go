package hw

import (
	"math"
	"sort"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/parallel"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/special"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/stats"
)

// Params menyimpan parameter model Hull-White 1-faktor dengan kurva forward datar f(0,t) = R0.
type Params struct {
	R0    float64
	A     float64
	Sigma float64
}

// DefaultParams mengembalikan parameter baku Hull-White sesi day003e.
func DefaultParams() Params {
	return Params{R0: 0.03, A: 0.08, Sigma: 0.01}
}

// P0 mengembalikan faktor diskonto kurva awal P(0, T) = exp(-R0 * T).
func (p Params) P0(T float64) float64 {
	return math.Exp(-p.R0 * T)
}

// B mengembalikan koefisien eksposur afin B(t, T) = (1 - exp(-a*(T - t))) / a.
func (p Params) B(t, T float64) float64 {
	dt := T - t
	if dt <= 0 {
		return 0
	}
	u := p.A * dt
	if math.Abs(u) < 1e-6 {
		return dt * (1.0 - 0.5*u + u*u/6.0 - u*u*u/24.0)
	}
	return -math.Expm1(-u) / p.A
}

// expm1OverU menghitung (1 - exp(-u)) / u secara stabil untuk u kecil.
func expm1OverU(u float64) float64 {
	if math.Abs(u) < 1e-5 {
		return 1.0 - 0.5*u + u*u/6.0 - u*u*u/24.0 + u*u*u*u/120.0
	}
	return -math.Expm1(-u) / u
}

// intOneMinusExpSq menghitung integral_0^u (1 - exp(-s))^2 ds secara stabil bebas pembatalan angka.
func intOneMinusExpSq(u float64) float64 {
	if u <= 0 {
		return 0
	}
	if u < 0.05 {
		// Deret Taylor: sum_{k>=3} (-1)^{k-3} (2^{k-1} - 2) / k! * u^k
		u3 := u * u * u
		return u3 * (1.0/3.0 - u*(0.25-u*(7.0/60.0-u*(1.0/24.0-u*(31.0/2520.0-u/320.0)))))
	}
	return u - 2.0*(-math.Expm1(-u)) + 0.5*(-math.Expm1(-2.0*u))
}

// VarX mengembalikan varians eksak inovasi x_{t+dt} | x_t pada langkah dt.
func (p Params) VarX(dt float64) float64 {
	if dt <= 0 {
		return 0
	}
	return p.Sigma * p.Sigma * dt * expm1OverU(2.0*p.A*dt)
}

// CovXI mengembalikan kovarians eksak Cov(x_{t+dt} - e^{-a*dt} x_t, int_t^{t+dt} x_s ds).
func (p Params) CovXI(dt float64) float64 {
	if dt <= 0 {
		return 0
	}
	e := expm1OverU(p.A*dt) * dt
	return 0.5 * p.Sigma * p.Sigma * e * e
}

// VarI mengembalikan varians eksak inovasi I_dt = int_0^dt tilde_x_s ds.
func (p Params) VarI(dt float64) float64 {
	if dt <= 0 {
		return 0
	}
	u := p.A * dt
	return (p.Sigma * p.Sigma / (p.A * p.A * p.A)) * intOneMinusExpSq(u)
}

// Alpha mengembalikan pergeseran deterministik alpha(t) = r_0 + sigma^2/(2 a^2) (1 - exp(-a*t))^2.
func (p Params) Alpha(t float64) float64 {
	return p.R0 + p.CovXI(t)
}

// AInt mengembalikan integral_0^t alpha(s) ds = r_0*t + 0.5 * VarI(t).
func (p Params) AInt(t float64) float64 {
	return p.R0*t + 0.5*p.VarI(t)
}

// DiscountFactor mengembalikan faktor diskonto stokastik D(0, t) = exp(-AInt(0, t) - I_t).
func (p Params) DiscountFactor(t, integralX float64) float64 {
	return math.Exp(-p.AInt(t) - integralX)
}

// CondDiscountGivenX mengembalikan E[D(0, t) | x_t] secara analitik dari distribusi bivariat (x_t, I_t).
func (p Params) CondDiscountGivenX(t, xt float64) float64 {
	if t <= 0 {
		return 1.0
	}
	vx := p.VarX(t)
	cxi := p.CovXI(t)
	if vx <= 0 {
		return p.P0(t)
	}
	ratio := cxi / vx
	return p.P0(t) * math.Exp(-ratio*xt-0.5*ratio*cxi)
}

// BondPrice mengembalikan harga obligasi kupon nol afin P(t, T; x_t).
func (p Params) BondPrice(t, T, xt float64) float64 {
	if T <= t {
		return 1.0
	}
	b := p.B(t, T)
	vx := p.VarX(t)
	cxi := p.CovXI(t)
	return (p.P0(T) / p.P0(t)) * math.Exp(-b*cxi-0.5*vx*b*b-b*xt)
}

// StepParams menyimpan koefisien Cholesky dan pergeseran untuk transisi eksak (x_t, I_t) pada langkah dt.
type StepParams struct {
	Dt     float64
	ExpADt float64
	B0Dt   float64
	VarX   float64
	CovXI  float64
	VarI   float64
	L11    float64
	L21    float64
	L22    float64
}

// NewStepParams menghitung faktor Cholesky matriks kovarians 2x2 untuk langkah dt > 0.
func (p Params) NewStepParams(dt float64) StepParams {
	if dt <= 0 {
		return StepParams{Dt: 0, ExpADt: 1}
	}
	vx := p.VarX(dt)
	cxi := p.CovXI(dt)
	vi := p.VarI(dt)
	l11 := math.Sqrt(math.Max(vx, 0))
	l21 := 0.0
	if l11 > 0 {
		l21 = cxi / l11
	}
	rem := vi - l21*l21
	if rem < 0 {
		rem = 0
	}
	return StepParams{
		Dt:     dt,
		ExpADt: math.Exp(-p.A * dt),
		B0Dt:   p.B(0, dt),
		VarX:   vx,
		CovXI:  cxi,
		VarI:   vi,
		L11:    l11,
		L21:    l21,
		L22:    math.Sqrt(rem),
	}
}

// Advance melakukan satu transisi Gaussian bivariat eksak (x_{t+dt}, I_{t+dt}) dari (xt, it)
// memakai dua peubah normal standar independen (z1, z2).
func (sp StepParams) Advance(xt, it, z1, z2 float64) (nextX, nextI float64) {
	nextX = xt*sp.ExpADt + sp.L11*z1
	nextI = it + xt*sp.B0Dt + sp.L21*z1 + sp.L22*z2
	return nextX, nextI
}

// MartingaleCheckGH memverifikasi identitas E[D(0, t) P(t, T; x_t)] = P(0, T) via kuadratur Gauss-Hermite 1D & 2D.
func (p Params) MartingaleCheckGH(t, T float64, nGH int) (lh1D, lh2D, p0T, err1D, err2D float64) {
	p0T = p.P0(T)
	if t <= 0 {
		p0 := p.BondPrice(0, T, 0)
		return p0, p0, p0T, math.Abs(p0 - p0T), math.Abs(p0 - p0T)
	}
	nodes := special.GaussHermiteStandardNormal(nGH)
	sp := p.NewStepParams(t)
	for _, nx := range nodes {
		xt := sp.L11 * nx.X
		pt := p.BondPrice(t, T, xt)
		lh1D += nx.W * p.CondDiscountGivenX(t, xt) * pt
		for _, ni := range nodes {
			it := sp.L21*nx.X + sp.L22*ni.X
			lh2D += nx.W * ni.W * p.DiscountFactor(t, it) * pt
		}
	}
	err1D = math.Abs(lh1D - p0T)
	err2D = math.Abs(lh2D - p0T)
	return lh1D, lh2D, p0T, err1D, err2D
}

// MartingaleCheckMC memverifikasi E[D(0, t) P(t, T)] vs P(0, T) via Monte Carlo paralel.
func (p Params) MartingaleCheckMC(t, T float64, nPaths, nChunks, nWorkers int, seed uint64) (mcMean, mcSE, p0T, zScore float64) {
	p0T = p.P0(T)
	sp := p.NewStepParams(t)
	type chunkRes struct {
		sum, sumSq float64
		count      int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkRes {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x101, uint64(c)))
		var cr chunkRes
		cr.count = cnt
		for k := 0; k < cnt; k++ {
			z1 := rg.NormFloat64()
			z2 := rg.NormFloat64()
			xt, it := sp.Advance(0, 0, z1, z2)
			val := p.DiscountFactor(t, it) * p.BondPrice(t, T, xt)
			cr.sum += val
			cr.sumSq += val * val
		}
		return cr
	})
	sum, sumSq := 0.0, 0.0
	tot := 0
	for _, r := range res {
		sum += r.sum
		sumSq += r.sumSq
		tot += r.count
	}
	mcMean, mcSE = stats.SummaryMoments(sum, sumSq, tot)
	if mcSE > 0 {
		zScore = math.Abs(mcMean-p0T) / mcSE
	}
	return mcMean, mcSE, p0T, zScore
}

// SwapSpec mendefinisikan swap suku bunga tahunan (payer: +1, receiver: -1).
type SwapSpec struct {
	Maturity float64
	PayDates []float64
	Tau      []float64
	StrikeK  float64
	IsPayer  bool
}

// NewAnnualSwap membangun swap dengan jadwal pembayaran tahunan 1..MaturityYears.
// Bila strikeK <= 0, tingkat par K(0) dihitung otomatis dari kurva awal.
func (p Params) NewAnnualSwap(maturityYears int, strikeK float64, isPayer bool) SwapSpec {
	dates := make([]float64, maturityYears)
	tau := make([]float64, maturityYears)
	for i := 0; i < maturityYears; i++ {
		dates[i] = float64(i + 1)
		tau[i] = 1.0
	}
	k := strikeK
	if k <= 0 {
		ann := 0.0
		for i, tPay := range dates {
			ann += tau[i] * p.P0(tPay)
		}
		k = (1.0 - p.P0(float64(maturityYears))) / ann
	}
	return SwapSpec{
		Maturity: float64(maturityYears),
		PayDates: dates,
		Tau:      tau,
		StrikeK:  k,
		IsPayer:  isPayer,
	}
}

// Value menghitung nilai mark-to-market swap berjalan pada waktu t dan faktor x_t:
// V_payer(t) = 1 - P(t, T_n) - K sum_{T_i > t} tau_i P(t, T_i),
// dengan konvensi eksak V(T_n) = 0 pada jatuh tempo.
func (s SwapSpec) Value(p Params, t, xt float64) float64 {
	if t >= s.Maturity-1e-12 {
		return 0
	}
	pvFixed := 0.0
	for i, payDate := range s.PayDates {
		if payDate > t+1e-12 {
			pvFixed += s.Tau[i] * p.BondPrice(t, payDate, xt)
		}
	}
	vPayer := 1.0 - p.BondPrice(t, s.Maturity, xt) - s.StrikeK*pvFixed
	if s.IsPayer {
		return vPayer
	}
	return -vPayer
}

// PositiveBoundaryX mencari akar unik x^* di mana V_payer(t, x^*) = 0.
func (s SwapSpec) PositiveBoundaryX(p Params, t float64) float64 {
	if t >= s.Maturity-1e-12 {
		return 0
	}
	f := func(x float64) float64 {
		pv := 0.0
		for i, pd := range s.PayDates {
			if pd > t+1e-12 {
				pv += s.Tau[i] * p.BondPrice(t, pd, x)
			}
		}
		return 1.0 - p.BondPrice(t, s.Maturity, x) - s.StrikeK*pv
	}
	df := func(x float64) float64 {
		d := p.B(t, s.Maturity) * p.BondPrice(t, s.Maturity, x)
		for i, pd := range s.PayDates {
			if pd > t+1e-12 {
				d += s.StrikeK * s.Tau[i] * p.B(t, pd) * p.BondPrice(t, pd, x)
			}
		}
		return d
	}
	lo, hi := -0.5, 0.5
	for f(lo) > 0 {
		lo *= 2.0
	}
	for f(hi) < 0 {
		hi *= 2.0
	}
	x := 0.5 * (lo + hi)
	for it := 0; it < 60; it++ {
		fx := f(x)
		if math.Abs(fx) < 1e-14 || (hi-lo) < 1e-14 {
			break
		}
		if fx < 0 {
			lo = x
		} else {
			hi = x
		}
		cand := x - fx/df(x)
		if cand > lo && cand < hi {
			x = cand
		} else {
			x = 0.5 * (lo + hi)
		}
	}
	return x
}

// ExpectedExposureSmooth menghitung EE(t) dan EEDisc(t) bebas galat titik sudut (kink)
// dengan membelah domain integrasi tepat pada batas latihan x^*(t) via kuadratur Gauss-Legendre.
func (s SwapSpec) ExpectedExposureSmooth(p Params, t float64, nQuad int) (ee, eeDisc, eeDisc2D float64) {
	if t <= 1e-12 || t >= s.Maturity-1e-12 {
		v0 := math.Max(s.Value(p, t, 0), 0)
		return v0, v0, v0
	}
	sp := p.NewStepParams(t)
	if sp.L11 <= 0 {
		v0 := math.Max(s.Value(p, t, 0), 0)
		return v0, v0, v0
	}
	zStar := s.PositiveBoundaryX(p, t) / sp.L11
	const zBound = 10.0
	var a, b float64
	if s.IsPayer {
		a = math.Max(zStar, -zBound)
		b = zBound
	} else {
		a = -zBound
		b = math.Min(zStar, zBound)
	}
	if a >= b {
		return 0, 0, 0
	}
	glNodes := special.GaussLegendreInterval(nQuad, a, b)
	ghNodes := special.GaussHermiteStandardNormal(nQuad)
	for _, gz := range glNodes {
		z := gz.X
		wZ := gz.W * special.NormPDF(z)
		xt := sp.L11 * z
		val := s.Value(p, t, xt)
		if val <= 0 {
			continue
		}
		ee += wZ * val
		eeDisc += wZ * p.CondDiscountGivenX(t, xt) * val
		for _, gi := range ghNodes {
			itVal := sp.L21*z + sp.L22*gi.X
			eeDisc2D += wZ * gi.W * p.DiscountFactor(t, itVal) * val
		}
	}
	return ee, eeDisc, eeDisc2D
}

// ExposurePoint menyimpan ringkasan eksposur di satu tanggal kisi t_j.
type ExposurePoint struct {
	T          float64
	EEGH       float64
	EEDiscGH   float64
	EEDiscGH2D float64
	EEMC       float64
	EESEMC     float64
	EEDiscMC   float64
	EEDiscSEMC float64
	PFE95MC    float64
	PFE99MC    float64
}

// ComputeSwapExposureProfile menghitung EE & EEDisc via kuadratur halus bebas-kink (1D & 2D)
// dan profil penuh (EE, EEDisc, PFE 95%, PFE 99%) via Monte Carlo exact Hull-White pada kisi tGrid.
func (p Params) ComputeSwapExposureProfile(swap SwapSpec, tGrid []float64, nGH, nPaths, nChunks, nWorkers int, seed uint64) (pts []ExposurePoint, epeGH, epeMC float64) {
	nT := len(tGrid)
	pts = make([]ExposurePoint, nT)

	for j, t := range tGrid {
		pts[j].T = t
		ee1, disc1, disc2 := swap.ExpectedExposureSmooth(p, t, nGH)
		pts[j].EEGH = ee1
		pts[j].EEDiscGH = disc1
		pts[j].EEDiscGH2D = disc2
	}

	steps := make([]StepParams, nT)
	for j := 1; j < nT; j++ {
		steps[j] = p.NewStepParams(tGrid[j] - tGrid[j-1])
	}

	type chunkExp struct {
		sumEE    []float64
		sumSqEE  []float64
		sumDisc  []float64
		sumSqDis []float64
		expMat   [][]float64
		count    int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkExp {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		ce := chunkExp{
			sumEE:    make([]float64, nT),
			sumSqEE:  make([]float64, nT),
			sumDisc:  make([]float64, nT),
			sumSqDis: make([]float64, nT),
			expMat:   make([][]float64, nT),
			count:    cnt,
		}
		for j := 0; j < nT; j++ {
			ce.expMat[j] = make([]float64, cnt)
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x201, uint64(c)))
		for k := 0; k < cnt; k++ {
			xt, it := 0.0, 0.0
			for j := 0; j < nT; j++ {
				if j > 0 {
					z1 := rg.NormFloat64()
					z2 := rg.NormFloat64()
					xt, it = steps[j].Advance(xt, it, z1, z2)
				}
				v := math.Max(swap.Value(p, tGrid[j], xt), 0)
				d := 1.0
				if j > 0 {
					d = p.DiscountFactor(tGrid[j], it)
				}
				dv := d * v
				ce.expMat[j][k] = v
				ce.sumEE[j] += v
				ce.sumSqEE[j] += v * v
				ce.sumDisc[j] += dv
				ce.sumSqDis[j] += dv * dv
			}
		}
		return ce
	})

	allExp := make([]float64, nPaths)
	for j := 0; j < nT; j++ {
		sEE, sqEE, sDis, sqDis := 0.0, 0.0, 0.0, 0.0
		offset := 0
		for _, r := range res {
			sEE += r.sumEE[j]
			sqEE += r.sumSqEE[j]
			sDis += r.sumDisc[j]
			sqDis += r.sumSqDis[j]
			copy(allExp[offset:offset+r.count], r.expMat[j])
			offset += r.count
		}
		pts[j].EEMC, pts[j].EESEMC = stats.SummaryMoments(sEE, sqEE, nPaths)
		pts[j].EEDiscMC, pts[j].EEDiscSEMC = stats.SummaryMoments(sDis, sqDis, nPaths)
		sort.Float64s(allExp)
		pts[j].PFE95MC = sortedQuantile(allExp, 0.95)
		pts[j].PFE99MC = sortedQuantile(allExp, 0.99)
	}

	totalT := tGrid[nT-1] - tGrid[0]
	if totalT > 0 {
		areaGH, areaMC := 0.0, 0.0
		for j := 1; j < nT; j++ {
			dt := tGrid[j] - tGrid[j-1]
			areaGH += 0.5 * (pts[j-1].EEGH + pts[j].EEGH) * dt
			areaMC += 0.5 * (pts[j-1].EEMC + pts[j].EEMC) * dt
		}
		epeGH = areaGH / totalT
		epeMC = areaMC / totalT
	}
	return pts, epeGH, epeMC
}

func sortedQuantile(sorted []float64, q float64) float64 {
	n := len(sorted)
	if n == 0 {
		return 0
	}
	pos := q * float64(n-1)
	lo := int(math.Floor(pos))
	hi := int(math.Ceil(pos))
	w := pos - float64(lo)
	return (1.0-w)*sorted[lo] + w*sorted[hi]
}

// SwaptionResult menyimpan hasil verifikasi swaption Jamshidian 2y x 8y terhadap GH dan MC.
type SwaptionResult struct {
	Strike        float64
	Label         string
	XStar         float64
	CheckSumCiXi  float64
	Jamshidian    float64
	QuadGH1D      float64
	QuadGH2D      float64
	AbsErrJamVsGH float64
	MCMean        float64
	MCSE          float64
	ZScoreMC      float64
}

// ForwardSwapRate menghitung tingkat swap forward K_fwd untuk swap mulai tStart dan bayar payDates.
func (p Params) ForwardSwapRate(tStart float64, payDates []float64) float64 {
	ann := 0.0
	for _, tPay := range payDates {
		ann += p.P0(tPay)
	}
	return (p.P0(tStart) - p.P0(payDates[len(payDates)-1])) / ann
}

// ZeroBondPut menghitung harga opsi put obligasi kupon nol P_put(0, te, T_i, X_i) di bawah Hull-White.
func (p Params) ZeroBondPut(te, ti, strikeBond float64) float64 {
	b := p.B(te, ti)
	sigP := b * math.Sqrt(p.VarX(te))
	p0Te := p.P0(te)
	p0Ti := p.P0(ti)
	h := math.Log(p0Ti/(p0Te*strikeBond))/sigP + 0.5*sigP
	return strikeBond*p0Te*special.NormCDF(-h+sigP) - p0Ti*special.NormCDF(-h)
}

// VerifySwaption menghitung harga swaption Eropa payer via dekomposisi Jamshidian, kuadratur GH 1D/2D, dan MC.
func (p Params) VerifySwaption(te float64, payDates []float64, strikes []float64, labels []string, nGH, nPaths, nChunks, nWorkers int, seed uint64) []SwaptionResult {
	nK := len(strikes)
	out := make([]SwaptionResult, nK)
	tn := payDates[len(payDates)-1]
	nodes := special.GaussHermiteStandardNormal(nGH)
	sp := p.NewStepParams(te)

	swapVal := func(kRate, x float64) float64 {
		pv := 0.0
		for _, tp := range payDates {
			pv += p.BondPrice(te, tp, x)
		}
		return 1.0 - p.BondPrice(te, tn, x) - kRate*pv
	}

	swapValDerivX := func(kRate, x float64) float64 {
		bn := p.B(te, tn)
		dv := bn * p.BondPrice(te, tn, x)
		for _, tp := range payDates {
			dv += kRate * p.B(te, tp) * p.BondPrice(te, tp, x)
		}
		return dv
	}

	for idx, kRate := range strikes {
		lo, hi := -0.5, 0.5
		for swapVal(kRate, lo) > 0 {
			lo *= 2.0
		}
		for swapVal(kRate, hi) < 0 {
			hi *= 2.0
		}
		x := 0.5 * (lo + hi)
		for it := 0; it < 60; it++ {
			fx := swapVal(kRate, x)
			if math.Abs(fx) < 1e-14 || (hi-lo) < 1e-14 {
				break
			}
			if fx < 0 {
				lo = x
			} else {
				hi = x
			}
			dfx := swapValDerivX(kRate, x)
			cand := x - fx/dfx
			if cand > lo && cand < hi {
				x = cand
			} else {
				x = 0.5 * (lo + hi)
			}
		}

		jam := 0.0
		checkSum := 0.0
		for j, tp := range payDates {
			ci := kRate
			if j == len(payDates)-1 {
				ci += 1.0
			}
			xi := p.BondPrice(te, tp, x)
			checkSum += ci * xi
			jam += ci * p.ZeroBondPut(te, tp, xi)
		}

		gh1, gh2 := 0.0, 0.0
		zStar := x / sp.L11
		const zBound = 10.0
		a := math.Max(zStar, -zBound)
		if a < zBound {
			glZ := special.GaussLegendreInterval(nGH, a, zBound)
			for _, gz := range glZ {
				z := gz.X
				wZ := gz.W * special.NormPDF(z)
				xt := sp.L11 * z
				payoff := swapVal(kRate, xt)
				if payoff <= 0 {
					continue
				}
				gh1 += wZ * p.CondDiscountGivenX(te, xt) * payoff
				for _, ni := range nodes {
					itVal := sp.L21*z + sp.L22*ni.X
					gh2 += wZ * ni.W * p.DiscountFactor(te, itVal) * payoff
				}
			}
		}

		out[idx] = SwaptionResult{
			Strike:        kRate,
			Label:         labels[idx],
			XStar:         x,
			CheckSumCiXi:  checkSum,
			Jamshidian:    jam,
			QuadGH1D:      gh1,
			QuadGH2D:      gh2,
			AbsErrJamVsGH: math.Abs(jam - gh2),
		}
	}

	type chunkSum struct {
		sums   []float64
		sumSqs []float64
		count  int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkSum {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cs := chunkSum{sums: make([]float64, nK), sumSqs: make([]float64, nK), count: cnt}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x301, uint64(c)))
		for k := 0; k < cnt; k++ {
			z1 := rg.NormFloat64()
			z2 := rg.NormFloat64()
			xt, it := sp.Advance(0, 0, z1, z2)
			disc := p.DiscountFactor(te, it)
			pn := p.BondPrice(te, tn, xt)
			pvFixed := 0.0
			for _, tp := range payDates {
				pvFixed += p.BondPrice(te, tp, xt)
			}
			base := 1.0 - pn
			for idx, kRate := range strikes {
				payoff := math.Max(base-kRate*pvFixed, 0)
				val := disc * payoff
				cs.sums[idx] += val
				cs.sumSqs[idx] += val * val
			}
		}
		return cs
	})

	for idx := range strikes {
		sum, sumSq := 0.0, 0.0
		tot := 0
		for _, r := range res {
			sum += r.sums[idx]
			sumSq += r.sumSqs[idx]
			tot += r.count
		}
		mean, se := stats.SummaryMoments(sum, sumSq, tot)
		out[idx].MCMean = mean
		out[idx].MCSE = se
		if se > 0 {
			out[idx].ZScoreMC = math.Abs(mean-out[idx].Jamshidian) / se
		}
	}
	return out
}

// CollateralCurvePoint menyimpan profil eksposur terkolateralisasi pada tanggal t_j untuk threshold H.
type CollateralCurvePoint struct {
	T          float64
	UncollatEE float64
	CollatEE   map[float64]float64
}

// ComputeCollateralExposureProfile menghitung EE terkolateralisasi dengan MPOR delta dan daftar threshold H
// memakai simulasi eksak dua-waktu (0 -> t_j - delta -> t_j).
func (p Params) ComputeCollateralExposureProfile(swap SwapSpec, tGrid []float64, delta float64, thresholds []float64, nPaths, nChunks, nWorkers int, seed uint64) (eeByH map[float64][]float64, eeDiscByH map[float64][]float64, epeByH map[float64]float64) {
	nT := len(tGrid)
	nH := len(thresholds)
	eeByH = make(map[float64][]float64, nH)
	eeDiscByH = make(map[float64][]float64, nH)
	epeByH = make(map[float64]float64, nH)
	for _, h := range thresholds {
		eeByH[h] = make([]float64, nT)
		eeDiscByH[h] = make([]float64, nT)
	}

	type pairStep struct {
		needLag bool
		sp1     StepParams
		sp2     StepParams
		spFull  StepParams
	}
	// Simulasi eksak per tanggal t_j via transisi eksak 0 -> t_j - delta -> t_j dengan CRN antar j dan H.
	pairs := make([]pairStep, nT)
	for j, t := range tGrid {
		if t <= 1e-12 {
			continue
		}
		if t <= delta+1e-12 {
			pairs[j] = pairStep{needLag: false, spFull: p.NewStepParams(t)}
		} else {
			pairs[j] = pairStep{
				needLag: true,
				sp1:     p.NewStepParams(t - delta),
				sp2:     p.NewStepParams(delta),
			}
		}
	}

	type chunkCol struct {
		sumEE   [][]float64
		sumDisc [][]float64
		count   int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkCol {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cc := chunkCol{
			sumEE:   make([][]float64, nH),
			sumDisc: make([][]float64, nH),
			count:   cnt,
		}
		for ih := 0; ih < nH; ih++ {
			cc.sumEE[ih] = make([]float64, nT)
			cc.sumDisc[ih] = make([]float64, nT)
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x401, uint64(c)))
		for k := 0; k < cnt; k++ {
			for j, t := range tGrid {
				if t <= 1e-12 || t >= swap.Maturity-1e-12 {
					continue
				}
				z1 := rg.NormFloat64()
				z2 := rg.NormFloat64()
				var xt, it, vLag float64
				if !pairs[j].needLag {
					xt, it = pairs[j].spFull.Advance(0, 0, z1, z2)
					vLag = swap.Value(p, 0, 0)
				} else {
					z3 := rg.NormFloat64()
					z4 := rg.NormFloat64()
					xLag, iLag := pairs[j].sp1.Advance(0, 0, z1, z2)
					vLag = swap.Value(p, t-delta, xLag)
					xt, it = pairs[j].sp2.Advance(xLag, iLag, z3, z4)
				}
				vt := swap.Value(p, t, xt)
				dt := p.DiscountFactor(t, it)
				for ih, h := range thresholds {
					cAmt := 0.0
					if !math.IsInf(h, 1) {
						cAmt = math.Max(vLag-h, 0)
					}
					exp := math.Max(vt-cAmt, 0)
					cc.sumEE[ih][j] += exp
					cc.sumDisc[ih][j] += dt * exp
				}
			}
		}
		return cc
	})

	totalT := tGrid[nT-1] - tGrid[0]
	for ih, h := range thresholds {
		for j := 0; j < nT; j++ {
			sEE, sDisc := 0.0, 0.0
			for _, r := range res {
				sEE += r.sumEE[ih][j]
				sDisc += r.sumDisc[ih][j]
			}
			eeByH[h][j] = sEE / float64(nPaths)
			eeDiscByH[h][j] = sDisc / float64(nPaths)
		}
		if totalT > 0 {
			area := 0.0
			for j := 1; j < nT; j++ {
				area += 0.5 * (eeByH[h][j-1] + eeByH[h][j]) * (tGrid[j] - tGrid[j-1])
			}
			epeByH[h] = area / totalT
		}
	}
	return eeByH, eeDiscByH, epeByH
}
