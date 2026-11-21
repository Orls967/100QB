package cir

import (
	"math"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/parallel"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/stats"
)

// Params menyimpan parameter intensitas gagal Cox-Ingersoll-Ross (CIR):
// d lambda_t = kappa (theta - lambda_t) dt + xi sqrt(lambda_t) dW_t.
type Params struct {
	Lambda0  float64
	Kappa    float64
	Theta    float64
	Xi       float64
	Recovery float64
}

// DefaultParams mengembalikan parameter baku intensitas gagal CIR untuk sesi day003e.
func DefaultParams() Params {
	return Params{
		Lambda0:  0.02,
		Kappa:    0.5,
		Theta:    0.025,
		Xi:       0.1,
		Recovery: 0.4,
	}
}

// LGD mengembalikan loss given default (1 - Recovery).
func (p Params) LGD() float64 {
	return 1.0 - p.Recovery
}

// FellerCheck memeriksa syarat Feller 2*kappa*theta >= xi^2.
func (p Params) FellerCheck() (lhs, rhs float64, satisfied bool) {
	lhs = 2.0 * p.Kappa * p.Theta
	rhs = p.Xi * p.Xi
	return lhs, rhs, lhs >= rhs
}

// AffineCoeffs menghitung koefisien afin A(T) dan B(T) untuk Q(T) = A(T) exp(-B(T) lambda_0).
func (p Params) AffineCoeffs(T float64) (A, B float64) {
	if T <= 0 {
		return 1.0, 0.0
	}
	gamma := math.Sqrt(p.Kappa*p.Kappa + 2.0*p.Xi*p.Xi)
	eGamma := math.Exp(gamma * T)
	em1 := eGamma - 1.0
	denom := (p.Kappa+gamma)*em1 + 2.0*gamma
	power := 2.0 * p.Kappa * p.Theta / (p.Xi * p.Xi)
	// Faktor eksponensial ditulis dalam bentuk log agar bebas luapan numerik.
	logBase := math.Log(2.0*gamma) + 0.5*(p.Kappa+gamma)*T - math.Log(denom)
	A = math.Exp(power * logBase)
	B = 2.0 * em1 / denom
	return A, B
}

// SurvivalProb menghitung probabilitas survival Q(T) = E[exp(-int_0^T lambda_s ds)] secara analitik.
func (p Params) SurvivalProb(T float64) float64 {
	if T <= 0 {
		return 1.0
	}
	A, B := p.AffineCoeffs(T)
	return A * math.Exp(-B*p.Lambda0)
}

// SurvivalDLambda0 menghitung turunan analitik dQ(T)/d(lambda_0) = -B(T) Q(T).
func (p Params) SurvivalDLambda0(T float64) float64 {
	A, B := p.AffineCoeffs(T)
	return -B * A * math.Exp(-B*p.Lambda0)
}

// SurvivalRiccatiRK4 mengintegrasikan sistem ODE Riccati secara independen memakai RK4
// pada nSteps langkah seragam untuk verifikasi rumus afin.
func (p Params) SurvivalRiccatiRK4(T float64, nSteps int) float64 {
	if T <= 0 {
		return 1.0
	}
	h := T / float64(nSteps)
	bVal := 0.0
	lnA := 0.0
	dB := func(b float64) float64 {
		return 1.0 - p.Kappa*b - 0.5*p.Xi*p.Xi*b*b
	}
	dLnA := func(b float64) float64 {
		return -p.Kappa * p.Theta * b
	}
	for step := 0; step < nSteps; step++ {
		k1b := dB(bVal)
		k1a := dLnA(bVal)

		b2 := bVal + 0.5*h*k1b
		k2b := dB(b2)
		k2a := dLnA(b2)

		b3 := bVal + 0.5*h*k2b
		k3b := dB(b3)
		k3a := dLnA(b3)

		b4 := bVal + h*k3b
		k4b := dB(b4)
		k4a := dLnA(b4)

		bVal += (h / 6.0) * (k1b + 2.0*k2b + 2.0*k3b + k4b)
		lnA += (h / 6.0) * (k1a + 2.0*k2a + 2.0*k3a + k4a)
	}
	return math.Exp(lnA - bVal*p.Lambda0)
}

// ExactStepPrecomp menyimpan konstanta transisi chi-kuadrat tak-sentral untuk langkah tetap dt.
type ExactStepPrecomp struct {
	Dt          float64
	C           float64
	HalfD       float64
	HalfExpRate float64 // c * exp(-kappa * dt)
}

// NewExactStepPrecomp menyiapkan konstanta campuran Poisson-Gamma untuk sampling eksak CIR.
func (p Params) NewExactStepPrecomp(dt float64) ExactStepPrecomp {
	expK := math.Exp(-p.Kappa * dt)
	c := 2.0 * p.Kappa / (p.Xi * p.Xi * (-math.Expm1(-p.Kappa * dt)))
	d := 4.0 * p.Kappa * p.Theta / (p.Xi * p.Xi)
	return ExactStepPrecomp{
		Dt:          dt,
		C:           c,
		HalfD:       0.5 * d,
		HalfExpRate: c * expK,
	}
}

// SampleNext menghasilkan lambda_{t+dt} | lambda_t secara eksak melalui campuran Poisson-Gamma.
func (es ExactStepPrecomp) SampleNext(lambdaT float64, rg *rng.RNG) float64 {
	if lambdaT < 0 {
		lambdaT = 0
	}
	poisMean := es.HalfExpRate * lambdaT
	nPois := rg.Poisson(poisMean)
	g := rg.Gamma(es.HalfD + float64(nPois))
	return g / es.C
}

// TransitionMoments mengembalikan momen bersyarat eksak E[lambda_{t+dt} | lambda_t] dan Var(lambda_{t+dt} | lambda_t).
func (p Params) TransitionMoments(lambdaT, dt float64) (mean, variance float64) {
	ek := math.Exp(-p.Kappa * dt)
	e2k := ek * ek
	oneMinusEk := -math.Expm1(-p.Kappa * dt)
	mean = p.Theta + (lambdaT-p.Theta)*ek
	variance = (lambdaT*p.Xi*p.Xi/p.Kappa)*(ek-e2k) + (p.Theta*p.Xi*p.Xi/(2.0*p.Kappa))*oneMinusEk*oneMinusEk
	return mean, variance
}

// StepFullTruncationEuler melakukan satu langkah skema Euler full-truncation (Lord dkk., 2010).
func (p Params) StepFullTruncationEuler(lambdaT, dt, sqrtDt, z float64) float64 {
	lPos := math.Max(lambdaT, 0)
	next := lambdaT + p.Kappa*(p.Theta-lPos)*dt + p.Xi*math.Sqrt(lPos)*sqrtDt*z
	if next < 0 {
		return 0
	}
	return next
}

// TrapezoidBiasRow merangkum estimasi MC survival & bias diskretisasi pada resolusi sub-langkah tertentu.
type TrapezoidBiasRow struct {
	Horizon         float64
	SubStepsPerYear int
	AnalyticQ       float64
	ExactCIRMC      float64
	ExactCIRSE      float64
	ExactCIRBias    float64
	EulerFTMC       float64
	EulerFTSE       float64
	EulerFTBias     float64
}

// VerifySurvivalMC menguji survival MC dengan integrasi trapesium int_0^T lambda_s ds
// untuk skema transisi eksak CIR dan Euler full-truncation pada daftar sub-langkah/tahun.
func (p Params) VerifySurvivalMC(horizons []float64, subStepsList []int, nPaths, nChunks, nWorkers int, seed uint64) []TrapezoidBiasRow {
	rows := make([]TrapezoidBiasRow, 0, len(horizons)*len(subStepsList))
	maxT := 0.0
	for _, T := range horizons {
		if T > maxT {
			maxT = T
		}
	}
	// Gunakan sub-grid bersarang (nested Brownian increments) berdasarkan resolusi tertinggi agar komparasi antarskala ketat.
	maxM := 0
	for _, m := range subStepsList {
		if m > maxM {
			maxM = m
		}
	}

	for _, m := range subStepsList {
		dt := 1.0 / float64(m)
		sqrtDt := math.Sqrt(dt)
		nTotalSteps := int(math.Round(maxT * float64(m)))
		es := p.NewExactStepPrecomp(dt)

		horizonSteps := make([]int, len(horizons))
		for hi, T := range horizons {
			horizonSteps[hi] = int(math.Round(T * float64(m)))
		}

		type chunkRes struct {
			sumEx, sumSqEx []float64
			sumEu, sumSqEu []float64
			count          int
		}
		perChunk := nPaths / nChunks
		rem := nPaths % nChunks
		res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkRes {
			cnt := perChunk
			if c < rem {
				cnt++
			}
			nH := len(horizons)
			cr := chunkRes{
				sumEx:   make([]float64, nH),
				sumSqEx: make([]float64, nH),
				sumEu:   make([]float64, nH),
				sumSqEu: make([]float64, nH),
				count:   cnt,
			}
			// Benih terikat pada resolusi m agar dua skema memakai chunk seed independen yang terkontrol,
			// atau memakai sub-langkah bersarang bila m membagi maxM.
			rgEx := rng.New(rng.DeriveChunkSeed(seed, 0x501^uint64(m), uint64(c)))
			rgEu := rng.New(rng.DeriveChunkSeed(seed, 0x502, uint64(c)))
			groupSize := 1
			if maxM%m == 0 {
				groupSize = maxM / m
			}

			for k := 0; k < cnt; k++ {
				lamEx := p.Lambda0
				lamEu := p.Lambda0
				intEx, intEu := 0.0, 0.0
				hIdx := 0
				for s := 1; s <= nTotalSteps; s++ {
					nextEx := es.SampleNext(lamEx, rgEx)
					intEx += 0.5 * (lamEx + nextEx) * dt
					lamEx = nextEx

					zAgg := 0.0
					for g := 0; g < groupSize; g++ {
						zAgg += rgEu.NormFloat64()
					}
					if groupSize > 1 {
						zAgg /= math.Sqrt(float64(groupSize))
					}
					nextEu := p.StepFullTruncationEuler(lamEu, dt, sqrtDt, zAgg)
					intEu += 0.5 * (lamEu + nextEu) * dt
					lamEu = nextEu

					for hIdx < nH && s == horizonSteps[hIdx] {
						qEx := math.Exp(-intEx)
						qEu := math.Exp(-intEu)
						cr.sumEx[hIdx] += qEx
						cr.sumSqEx[hIdx] += qEx * qEx
						cr.sumEu[hIdx] += qEu
						cr.sumSqEu[hIdx] += qEu * qEu
						hIdx++
					}
				}
			}
			return cr
		})

		for hi, T := range horizons {
			sEx, sqEx, sEu, sqEu := 0.0, 0.0, 0.0, 0.0
			for _, r := range res {
				sEx += r.sumEx[hi]
				sqEx += r.sumSqEx[hi]
				sEu += r.sumEu[hi]
				sqEu += r.sumSqEu[hi]
			}
			mEx, seEx := stats.SummaryMoments(sEx, sqEx, nPaths)
			mEu, seEu := stats.SummaryMoments(sEu, sqEu, nPaths)
			qAnal := p.SurvivalProb(T)
			rows = append(rows, TrapezoidBiasRow{
				Horizon:         T,
				SubStepsPerYear: m,
				AnalyticQ:       qAnal,
				ExactCIRMC:      mEx,
				ExactCIRSE:      seEx,
				ExactCIRBias:    mEx - qAnal,
				EulerFTMC:       mEu,
				EulerFTSE:       seEu,
				EulerFTBias:     mEu - qAnal,
			})
		}
	}
	return rows
}

// CDSResult merangkum par spread CDS kuartalan analitik vs Monte Carlo.
type CDSResult struct {
	Maturity   float64
	AnalPVProt float64
	AnalRPV01  float64
	AnalSpread float64
	MCPVProt   float64
	MCRPV01    float64
	MCSpread   float64
	SEPVProt   float64
	SERPV01    float64
	ZScoreProt float64
}

// VerifyParCDS menghitung par spread CDS kuartalan (termasuk premi akrual di tengah kuartal)
// untuk jatuh tempo T dengan kurva diskonto datar exp(-r0 * t).
func (p Params) VerifyParCDS(maturities []float64, r0 float64, nPaths, nChunks, nWorkers int, seed uint64) []CDSResult {
	out := make([]CDSResult, len(maturities))
	maxT := 0.0
	for _, T := range maturities {
		if T > maxT {
			maxT = T
		}
	}
	nQ := int(math.Round(maxT * 4.0))
	dtQ := 0.25
	es := p.NewExactStepPrecomp(dtQ)

	qAnal := make([]float64, nQ+1)
	d0Q := make([]float64, nQ+1)
	d0Mid := make([]float64, nQ+1)
	for j := 0; j <= nQ; j++ {
		t := float64(j) * dtQ
		qAnal[j] = p.SurvivalProb(t)
		d0Q[j] = math.Exp(-r0 * t)
		if j > 0 {
			d0Mid[j] = math.Exp(-r0 * (t - 0.5*dtQ))
		}
	}

	type chunkCDS struct {
		sumProt, sumSqProt []float64
		sumRPV, sumSqRPV   []float64
		count              int
	}
	nM := len(maturities)
	targetQ := make([]int, nM)
	for im, T := range maturities {
		targetQ[im] = int(math.Round(T * 4.0))
	}

	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	dtSub := dtQ / 5.0
	esSub := p.NewExactStepPrecomp(dtSub)
	_ = es
	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkCDS {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		cc := chunkCDS{
			sumProt:   make([]float64, nM),
			sumSqProt: make([]float64, nM),
			sumRPV:    make([]float64, nM),
			sumSqRPV:  make([]float64, nM),
			count:     cnt,
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x601, uint64(c)))
		for k := 0; k < cnt; k++ {
			lam := p.Lambda0
			intLam := 0.0
			sPrev := 1.0
			pvProt, rpv01 := 0.0, 0.0
			mIdx := 0
			for j := 1; j <= nQ; j++ {
				for sub := 0; sub < 5; sub++ {
					nextLam := esSub.SampleNext(lam, rg)
					intLam += 0.5 * (lam + nextLam) * dtSub
					lam = nextLam
				}
				sCurr := math.Exp(-intLam)
				dpd := sPrev - sCurr
				pvProt += p.LGD() * d0Mid[j] * dpd
				rpv01 += dtQ * d0Q[j] * 0.5 * (sPrev + sCurr)
				sPrev = sCurr

				for mIdx < nM && j == targetQ[mIdx] {
					cc.sumProt[mIdx] += pvProt
					cc.sumSqProt[mIdx] += pvProt * pvProt
					cc.sumRPV[mIdx] += rpv01
					cc.sumSqRPV[mIdx] += rpv01 * rpv01
					mIdx++
				}
			}
		}
		return cc
	})

	for im, T := range maturities {
		mQ := targetQ[im]
		analProt, analRPV := 0.0, 0.0
		for j := 1; j <= mQ; j++ {
			analProt += p.LGD() * d0Mid[j] * (qAnal[j-1] - qAnal[j])
			analRPV += dtQ * d0Q[j] * 0.5 * (qAnal[j-1] + qAnal[j])
		}
		sP, sqP, sR, sqR := 0.0, 0.0, 0.0, 0.0
		for _, r := range res {
			sP += r.sumProt[im]
			sqP += r.sumSqProt[im]
			sR += r.sumRPV[im]
			sqR += r.sumSqRPV[im]
		}
		mcProt, seProt := stats.SummaryMoments(sP, sqP, nPaths)
		mcRPV, seRPV := stats.SummaryMoments(sR, sqR, nPaths)
		zP := 0.0
		if seProt > 0 {
			zP = math.Abs(mcProt-analProt) / seProt
		}
		out[im] = CDSResult{
			Maturity:   T,
			AnalPVProt: analProt,
			AnalRPV01:  analRPV,
			AnalSpread: analProt / analRPV,
			MCPVProt:   mcProt,
			MCRPV01:    mcRPV,
			MCSpread:   mcProt / mcRPV,
			SEPVProt:   seProt,
			SERPV01:    seRPV,
			ZScoreProt: zP,
		}
	}
	return out
}
