package reinsurance

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/panjer"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/special"
)

// ReinsurerLoading adalah loading premi reasuradur excess-of-loss (xi = 0.4).
const ReinsurerLoading = 0.4

// RetentionParams menyimpan parameter premi bersih dan retensi untuk klaim tertahan min(X, M).
type RetentionParams struct {
	Dist           claims.Dist
	M              float64
	MuRetained     float64
	ExcessMean     float64
	ReinsPrem      float64
	NetPremiumC    float64
	RhoRetained    float64
	SafetyLoading  float64
	ProfitPositive bool
}

// EvaluateRetention menghitung E[min(X, M)], E[(X-M)_+], premi bersih cedent c_M,
// dan syarat kondisi profit c_M > lambda * E[min(X, M)].
func EvaluateRetention(d claims.Dist, m float64) RetentionParams {
	var muM, excM float64
	if math.IsInf(m, 1) || m > 500.0 {
		muM = claims.Mu
		excM = 0
	} else {
		switch d.Kind {
		case claims.KindExponential:
			em := math.Exp(-d.ExpRate * m)
			muM = (1.0 - em) / d.ExpRate
			excM = em / d.ExpRate
		case claims.KindHyperexp2:
			p1 := d.HyperP1
			p2 := 1.0 - p1
			muM = p1*(1.0-math.Exp(-d.HyperB1*m))/d.HyperB1 + p2*(1.0-math.Exp(-d.HyperB2*m))/d.HyperB2
			excM = p1*math.Exp(-d.HyperB1*m)/d.HyperB1 + p2*math.Exp(-d.HyperB2*m)/d.HyperB2
		default:
			muM = claims.Mu * (1.0 - d.LadderSurvival(m))
			excM = claims.Mu * d.LadderSurvival(m)
		}
	}
	reinsPrem := (1.0 + ReinsurerLoading) * claims.Lambda * excM
	netC := claims.Premium - reinsPrem
	rhoM := claims.Lambda * muM / netC
	netLoad := (netC - claims.Lambda*muM) / (claims.Lambda * muM)
	return RetentionParams{
		Dist:           d,
		M:              m,
		MuRetained:     muM,
		ExcessMean:     excM,
		ReinsPrem:      reinsPrem,
		NetPremiumC:    netC,
		RhoRetained:    rhoM,
		SafetyLoading:  netLoad,
		ProfitPositive: netC > claims.Lambda*muM && rhoM > 0 && rhoM < 1.0,
	}
}

// MinValidRetention mencari batas bawah M_min sedemikian sehingga c_M = lambda * E[min(X, M)],
// yaitu integral_M^inf F_bar(x) dx = ((xi - theta)/xi) * mu = 0.5 * mu.
func MinValidRetention(d claims.Dist) (float64, error) {
	if d.Kind == claims.KindExponential {
		return math.Log(2.0), nil
	}
	return special.BrentRoot(func(m float64) float64 {
		rp := EvaluateRetention(d, m)
		return rp.NetPremiumC - claims.Lambda*rp.MuRetained
	}, 0.01, 50.0, 1e-13)
}

// RetainedLadderSurvival mengembalikan survival tinggi tangga untuk klaim tertahan X_M = min(X, M):
// F_bar_L^(M)(x) = (1 / mu_M) * integral_x^M F_bar(y) dy untuk 0 <= x < M, dan 0 untuk x >= M.
func RetainedLadderSurvival(d claims.Dist, m, muM, x float64) float64 {
	if x <= 0 {
		return 1.0
	}
	if !math.IsInf(m, 1) && x >= m {
		return 0.0
	}
	raw := claims.Mu * (d.LadderSurvival(x) - d.LadderSurvival(m)) / muM
	if raw < 0 {
		return 0
	}
	if raw > 1 {
		return 1
	}
	return raw
}

// evalGm menghitung L^{-1}{1 / (s * (s + alpha)^k)}(x) = integral_0^x t^{k-1}/(k-1)! e^{-alpha t} dt
// secara numerik akurat untuk alpha > 0, alpha == 0, maupun alpha < 0 (untuk u/M moderat).
func evalGm(k int, x, alpha float64) float64 {
	if k == 0 {
		return 1.0
	}
	if x <= 0 {
		return 0.0
	}
	ax := alpha * x
	if math.Abs(ax) < 0.8 {
		// Deret positif / teratur: x^k * e^{-alpha x} * sum_{j=0}^{35} (alpha x)^j / (k+j)!
		term := 1.0
		for j := 1; j <= k; j++ {
			term *= x / float64(j)
		}
		acc := term
		for j := 1; j <= 35; j++ {
			term *= ax / float64(k+j)
			acc += term
			if math.Abs(term) < 1e-17*math.Abs(acc) {
				break
			}
		}
		return math.Exp(-ax) * acc
	}
	// Identitas eksak (1/alpha^k) * [1 - e^{-alpha x} sum_{j=0}^{k-1} (alpha x)^j / j!]
	sumPoly := 1.0
	powTerm := 1.0
	for j := 1; j < k; j++ {
		powTerm *= ax / float64(j)
		sumPoly += powTerm
	}
	return (1.0 - math.Exp(-ax)*sumPoly) / math.Pow(alpha, float64(k))
}

// ExactRetainedExpRuin mengevaluasi rumus tertutup metode langkah berhingga psi_M(u)
// untuk klaim Eksponensial(1) yang dipotong pada retensi M (akurat sangat tinggi saat u/M <= 15).
func ExactRetainedExpRuin(u, m float64) (float64, error) {
	if u < 0 {
		return 1.0, nil
	}
	if math.IsInf(m, 1) || m > 500.0 || m >= u+50.0 {
		rp := EvaluateRetention(claims.Exponential(), m)
		if !rp.ProfitPositive {
			return 1.0, fmt.Errorf("negative net safety loading at M=%g", m)
		}
		return claims.Rho * math.Exp(-(claims.Theta/((1.0+claims.Theta)*claims.Mu))*u), nil
	}
	rp := EvaluateRetention(claims.Exponential(), m)
	if !rp.ProfitPositive {
		return 1.0, fmt.Errorf("net safety loading not positive for M=%g (c_M=%g, lam*muM=%g)", m, rp.NetPremiumC, rp.MuRetained)
	}
	aM := claims.Lambda / rp.NetPremiumC
	alpha := 1.0 - aM
	kMax := int(math.Floor(u / m))
	sumPhi := 0.0
	emNeg := math.Exp(-m)
	weight := 1.0
	for k := 0; k <= kMax; k++ {
		xk := u - float64(k)*m
		gk := evalGm(k, xk, alpha)
		gk1 := evalGm(k+1, xk, alpha)
		sumPhi += weight * (gk + aM*gk1)
		weight *= -aM * emNeg
	}
	phi := (1.0 - rp.RhoRetained) * sumPhi
	psi := 1.0 - phi
	if psi < 0 {
		psi = 0
	} else if psi > 1 {
		psi = 1
	}
	return psi, nil
}

// PanjerPsi evaluasi [low, mid, high] untuk psi_M(u) dengan langkah h (standar h=0.01 atau 0.02).
func PanjerPsi(d claims.Dist, m, u, h float64) (low, mid, high float64) {
	rp := EvaluateRetention(d, m)
	if !rp.ProfitPositive {
		return 1.0, 1.0, 1.0
	}
	ser := panjer.ComputeBounds(rp.RhoRetained, func(x float64) float64 {
		return RetainedLadderSurvival(d, m, rp.MuRetained, x)
	}, u, h)
	lo, md, hi, _, _, _ := ser.At(u)
	return lo, md, hi
}

// OptimalRetentionResult menyimpan hasil pencarian retensi M^* optimal pada u = 10 beserta pindai kisi.
type OptimalRetentionResult struct {
	DistCode        string
	DistName        string
	U               float64
	MMin            float64
	OptimalM        float64
	MinPsiMid       float64
	MinPsiClosed    float64
	PsiNoReins      float64
	MaxClosedVsBand bool
	ScanM           []float64
	ScanPsiMid      []float64
	ScanPsiLow      []float64
	ScanPsiHigh     []float64
	ScanPsiClosed   []float64
	InfLimitDiff    float64
}

// FindOptimalRetention meminimalkan psi_M(u) di atas M > M_min dengan pencarian Golden Section + pindai kisi halus.
func FindOptimalRetention(d claims.Dist, u, h float64) (OptimalRetentionResult, error) {
	mMin, err := MinValidRetention(d)
	if err != nil {
		return OptimalRetentionResult{}, err
	}
	evalObj := func(m float64) float64 {
		if d.Kind == claims.KindExponential {
			v, errC := ExactRetainedExpRuin(u, m)
			if errC == nil {
				return v
			}
		}
		_, mid, _ := PanjerPsi(d, m, u, h)
		return mid
	}

	// Pindai kasar lebih dulu untuk mengurung lembah minimum global secara pasti, lalu pertajam dengan GoldenSectionMin.
	gridM := []float64{
		mMin + 0.08, mMin + 0.15, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5,
		2.75, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 30.0,
	}
	bestIdx := 0
	bestVal := math.Inf(1)
	vals := make([]float64, len(gridM))
	for i, m := range gridM {
		if m <= mMin+0.02 {
			vals[i] = 1.0
			continue
		}
		vals[i] = evalObj(m)
		if vals[i] < bestVal {
			bestVal = vals[i]
			bestIdx = i
		}
	}
	loIdx := bestIdx - 1
	if loIdx < 0 {
		loIdx = 0
	}
	hiIdx := bestIdx + 1
	if hiIdx >= len(gridM) {
		hiIdx = len(gridM) - 1
	}
	mStar, fStar := special.GoldenSectionMin(evalObj, gridM[loIdx], gridM[hiIdx], 1e-4)

	scanM := make([]float64, 0, len(gridM)+1)
	scanMid := make([]float64, 0, len(gridM)+1)
	scanLow := make([]float64, 0, len(gridM)+1)
	scanHigh := make([]float64, 0, len(gridM)+1)
	scanClosed := make([]float64, 0, len(gridM)+1)
	bracketOK := true

	for _, m := range gridM {
		if m <= mMin+0.02 {
			continue
		}
		lo, md, hi := PanjerPsi(d, m, u, h)
		cl := math.NaN()
		if d.Kind == claims.KindExponential {
			cl, _ = ExactRetainedExpRuin(u, m)
			if !(lo <= cl+1e-12 && cl <= hi+1e-12) {
				bracketOK = false
			}
		}
		scanM = append(scanM, m)
		scanMid = append(scanMid, md)
		scanLow = append(scanLow, lo)
		scanHigh = append(scanHigh, hi)
		scanClosed = append(scanClosed, cl)
	}

	_, psiNoReinsMid, _ := PanjerPsi(d, math.Inf(1), u, h)
	_, psiHugeM, _ := PanjerPsi(d, 1000.0, u, h)
	minClosed := math.NaN()
	if d.Kind == claims.KindExponential {
		minClosed, _ = ExactRetainedExpRuin(u, mStar)
	}

	return OptimalRetentionResult{
		DistCode:        d.Code,
		DistName:        d.Name,
		U:               u,
		MMin:            mMin,
		OptimalM:        mStar,
		MinPsiMid:       fStar,
		MinPsiClosed:    minClosed,
		PsiNoReins:      psiNoReinsMid,
		MaxClosedVsBand: bracketOK,
		ScanM:           scanM,
		ScanPsiMid:      scanMid,
		ScanPsiLow:      scanLow,
		ScanPsiHigh:     scanHigh,
		ScanPsiClosed:   scanClosed,
		InfLimitDiff:    math.Abs(psiHugeM - psiNoReinsMid),
	}, nil
}

// CapitalRequirementRow merangkum modal awal minimum u^*(eps) dengan vs tanpa reinsurans.
type CapitalRequirementRow struct {
	DistCode       string
	DistName       string
	Epsilon        float64
	UStarNoReins   float64
	UStarWithReins float64
	RetentionM     float64
	CapitalSaving  float64
	SavingPercent  float64
	PsiCheckNo     float64
	PsiCheckWith   float64
}

// SolveMinimumCapital mencari u^*(eps) pemecah psi_M(u) = eps memakai Brent di atas interpolasi linear
// batas Panjer beresolusi halus h = 0.02 (dan fungsi tertutup eksak bila M = Inf).
func SolveMinimumCapital(d claims.Dist, mOpt float64, epsilons []float64, h float64) ([]CapitalRequirementRow, error) {
	rpOpt := EvaluateRetention(d, mOpt)
	const uCeiling = 110.0
	serNo := panjer.ComputeBounds(claims.Rho, d.LadderSurvival, uCeiling, h)
	serOpt := panjer.ComputeBounds(rpOpt.RhoRetained, func(x float64) float64 {
		return RetainedLadderSurvival(d, mOpt, rpOpt.MuRetained, x)
	}, uCeiling, h)

	interpMid := func(ser panjer.Series, u float64) float64 {
		if u <= 0 {
			return ser.High[0]
		}
		idx := int(math.Floor(u / ser.H))
		if idx >= ser.MMax {
			_, md, _, _, _, _ := ser.At(float64(ser.MMax) * ser.H)
			return md
		}
		f := (u - float64(idx)*ser.H) / ser.H
		v0 := 0.5 * (ser.Low[idx] + ser.High[idx])
		v1 := 0.5 * (ser.Low[idx+1] + ser.High[idx+1])
		return (1.0-f)*v0 + f*v1
	}

	out := make([]CapitalRequirementRow, len(epsilons))
	for i, eps := range epsilons {
		uNo, err := special.BrentRoot(func(u float64) float64 {
			return interpMid(serNo, u) - eps
		}, 0.0, uCeiling, 1e-8)
		if err != nil {
			return nil, err
		}
		uWith, err := special.BrentRoot(func(u float64) float64 {
			return interpMid(serOpt, u) - eps
		}, 0.0, uCeiling, 1e-8)
		if err != nil {
			return nil, err
		}
		sav := uNo - uWith
		pct := 0.0
		if uNo > 0 {
			pct = 100.0 * sav / uNo
		}
		out[i] = CapitalRequirementRow{
			DistCode:       d.Code,
			DistName:       d.Name,
			Epsilon:        eps,
			UStarNoReins:   uNo,
			UStarWithReins: uWith,
			RetentionM:     mOpt,
			CapitalSaving:  sav,
			SavingPercent:  pct,
			PsiCheckNo:     interpMid(serNo, uNo),
			PsiCheckWith:   interpMid(serOpt, uWith),
		}
	}
	return out, nil
}
