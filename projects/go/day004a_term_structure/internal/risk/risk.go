package risk

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
)

// StandardKeyRateNodes adalah simpul zero-rate standar untuk analisis Key-Rate Duration.
var StandardKeyRateNodes = []float64{1, 2, 3, 5, 7, 10, 20, 30}

// BondCashFlows membangun jadwal waktu dan arus kas obligasi kupon tetap (nosional = 1).
func BondCashFlows(maturity, couponRate, freq float64) (times, cfs []float64) {
	n := int(math.Round(maturity * freq))
	dt := 1.0 / freq
	times = make([]float64, n)
	cfs = make([]float64, n)
	for i := 0; i < n; i++ {
		times[i] = float64(i+1) * dt
		cfs[i] = couponRate * dt
	}
	if n > 0 {
		cfs[n-1] += 1.0
	}
	return times, cfs
}

// PriceBondWithShift menghitung harga obligasi dengan pergeseran paralel atau bobot kurva zero.
func PriceBondWithShift(curve ns.Params, times, cfs []float64, shiftFn func(tau float64) float64) float64 {
	pv := 0.0
	for i, t := range times {
		y := curve.Yield(t)
		if shiftFn != nil {
			y += shiftFn(t)
		}
		pv += cfs[i] * math.Exp(-t*y)
	}
	return pv
}

// SolveYTM mencari yield-to-maturity majemuk dengan frekuensi freq melalui metode Newton.
func SolveYTM(price float64, times, cfs []float64, freq float64) float64 {
	y := 0.04
	for iter := 0; iter < 60; iter++ {
		pv, dpv := 0.0, 0.0
		base := 1.0 + y/freq
		for i, t := range times {
			m := t * freq
			df := math.Pow(base, -m)
			pv += cfs[i] * df
			dpv -= cfs[i] * t * df / base
		}
		diff := pv - price
		if math.Abs(diff) < 1e-14 {
			break
		}
		y -= diff / dpv
	}
	return y
}

// BondRiskMetrics merangkum harga, durasi, konveksitas, dan perbandingan bump-and-reprice.
type BondRiskMetrics struct {
	Price                  float64
	YTM                    float64
	MacDurationYTM         float64
	ModDurationYTMClosed   float64
	ModDurationYTMFD       float64
	ModDurationYTMRelErr   float64
	ConvexityYTMClosed     float64
	ConvexityYTMFD         float64
	ConvexityYTMRelErr     float64
	EffDurationZeroClosed  float64
	EffDurationZeroFD      float64
	EffDurationZeroRelErr  float64
	EffConvexityZeroClosed float64
	EffConvexityZeroFD     float64
	EffConvexityZeroRelErr float64
}

// AnalyzeBondRisk menghitung metrik risiko bentuk tertutup vs bump-and-reprice pusat.
// Langkah bumpDur (~1e-6) dipakai untuk turunan pertama (durasi) dan bumpConv (~1e-4)
// dipakai untuk turunan kedua (konveksitas) guna menghindari pembatalan roundoff O(eps/h^2).
func AnalyzeBondRisk(curve ns.Params, maturity, couponRate, freq, bumpDur, bumpConv float64) BondRiskMetrics {
	times, cfs := BondCashFlows(maturity, couponRate, freq)
	price := PriceBondWithShift(curve, times, cfs, nil)

	// 1. Durasi & konveksitas efektif terhadap pergeseran paralel kurva zero kontinu
	sumT, sumT2 := 0.0, 0.0
	for i, t := range times {
		pvCF := cfs[i] * math.Exp(-t*curve.Yield(t))
		sumT += t * pvCF
		sumT2 += t * t * pvCF
	}
	effDurClosed := sumT / price
	effConvClosed := sumT2 / price

	pUp := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return bumpDur })
	pDn := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return -bumpDur })
	effDurFD := -(pUp - pDn) / (2.0 * bumpDur * price)
	pUpC := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return bumpConv })
	pDnC := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return -bumpConv })
	effConvFD := (pUpC - 2.0*price + pDnC) / (bumpConv * bumpConv * price)

	// 2. Durasi Macaulay/modifikasi & konveksitas terhadap YTM majemuk freq
	ytm := SolveYTM(price, times, cfs, freq)
	ytmPrice := func(y float64) float64 {
		base := 1.0 + y/freq
		s := 0.0
		for i, t := range times {
			s += cfs[i] * math.Pow(base, -t*freq)
		}
		return s
	}
	macSum, convSum := 0.0, 0.0
	base := 1.0 + ytm/freq
	for i, t := range times {
		pvCF := cfs[i] * math.Pow(base, -t*freq)
		macSum += t * pvCF
		convSum += t * (t + 1.0/freq) * pvCF
	}
	macDur := macSum / price
	modDurClosed := macDur / base
	convYTMClosed := convSum / (price * base * base)

	pyUp := ytmPrice(ytm + bumpDur)
	pyDn := ytmPrice(ytm - bumpDur)
	modDurFD := -(pyUp - pyDn) / (2.0 * bumpDur * price)
	pyUpC := ytmPrice(ytm + bumpConv)
	pyDnC := ytmPrice(ytm - bumpConv)
	convYTMFD := (pyUpC - 2.0*price + pyDnC) / (bumpConv * bumpConv * price)

	return BondRiskMetrics{
		Price:                  price,
		YTM:                    ytm,
		MacDurationYTM:         macDur,
		ModDurationYTMClosed:   modDurClosed,
		ModDurationYTMFD:       modDurFD,
		ModDurationYTMRelErr:   math.Abs(modDurFD-modDurClosed) / math.Abs(modDurClosed),
		ConvexityYTMClosed:     convYTMClosed,
		ConvexityYTMFD:         convYTMFD,
		ConvexityYTMRelErr:     math.Abs(convYTMFD-convYTMClosed) / math.Abs(convYTMClosed),
		EffDurationZeroClosed:  effDurClosed,
		EffDurationZeroFD:      effDurFD,
		EffDurationZeroRelErr:  math.Abs(effDurFD-effDurClosed) / math.Abs(effDurClosed),
		EffConvexityZeroClosed: effConvClosed,
		EffConvexityZeroFD:     effConvFD,
		EffConvexityZeroRelErr: math.Abs(effConvFD-effConvClosed) / math.Abs(effConvClosed),
	}
}

// TentWeight mengembalikan bobot segitiga KRD untuk simpul ke-idx pada grid simpul terurut.
// Jumlah seluruh bobot pada setiap tau > 0 adalah identik 1.
func TentWeight(nodes []float64, idx int, tau float64) float64 {
	m := len(nodes)
	if idx == 0 {
		if tau <= nodes[0] {
			return 1.0
		}
		if m > 1 && tau < nodes[1] {
			return (nodes[1] - tau) / (nodes[1] - nodes[0])
		}
		return 0.0
	}
	if idx == m-1 {
		if tau >= nodes[m-1] {
			return 1.0
		}
		if tau > nodes[m-2] {
			return (tau - nodes[m-2]) / (nodes[m-1] - nodes[m-2])
		}
		return 0.0
	}
	if tau <= nodes[idx-1] || tau >= nodes[idx+1] {
		return 0.0
	}
	if tau <= nodes[idx] {
		return (tau - nodes[idx-1]) / (nodes[idx] - nodes[idx-1])
	}
	return (nodes[idx+1] - tau) / (nodes[idx+1] - nodes[idx])
}

// KRDReport merangkum profil key-rate duration dan konsistensi penjumlahan dengan durasi paralel.
type KRDReport struct {
	Nodes           []float64
	KRDFD           []float64
	KRDClosed       []float64
	DollarKRDFD     []float64
	DollarKRDClosed []float64
	SumKRDFD        float64
	SumKRDClosed    float64
	EffParallelFD   float64
	EffParallelCl   float64
	SumVsParFDErr   float64
	SumVsParClErr   float64
}

// ComputeKeyRateDurations menghitung KRD tertutup dan bump-and-reprice pusat untuk arus kas obligasi.
func ComputeKeyRateDurations(curve ns.Params, times, cfs, nodes []float64, bump float64) KRDReport {
	price := PriceBondWithShift(curve, times, cfs, nil)
	krdFD := make([]float64, len(nodes))
	krdCl := make([]float64, len(nodes))
	dkrdFD := make([]float64, len(nodes))
	dkrdCl := make([]float64, len(nodes))
	sumFD, sumCl := 0.0, 0.0

	for idx := range nodes {
		pUp := PriceBondWithShift(curve, times, cfs, func(t float64) float64 {
			return bump * TentWeight(nodes, idx, t)
		})
		pDn := PriceBondWithShift(curve, times, cfs, func(t float64) float64 {
			return -bump * TentWeight(nodes, idx, t)
		})
		krdFD[idx] = -(pUp - pDn) / (2.0 * bump * price)
		dkrdFD[idx] = krdFD[idx] * price

		num := 0.0
		for j, t := range times {
			num += t * TentWeight(nodes, idx, t) * cfs[j] * math.Exp(-t*curve.Yield(t))
		}
		krdCl[idx] = num / price
		dkrdCl[idx] = num
		sumFD += krdFD[idx]
		sumCl += krdCl[idx]
	}

	pParUp := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return bump })
	pParDn := PriceBondWithShift(curve, times, cfs, func(float64) float64 { return -bump })
	effParFD := -(pParUp - pParDn) / (2.0 * bump * price)
	sumT := 0.0
	for j, t := range times {
		sumT += t * cfs[j] * math.Exp(-t*curve.Yield(t))
	}
	effParCl := sumT / price

	return KRDReport{
		Nodes:           nodes,
		KRDFD:           krdFD,
		KRDClosed:       krdCl,
		DollarKRDFD:     dkrdFD,
		DollarKRDClosed: dkrdCl,
		SumKRDFD:        sumFD,
		SumKRDClosed:    sumCl,
		EffParallelFD:   effParFD,
		EffParallelCl:   effParCl,
		SumVsParFDErr:   math.Abs(sumFD - effParFD),
		SumVsParClErr:   math.Abs(sumCl - effParCl),
	}
}

// SwapPV01Report merangkum tingkat par dan PV01/DV01 swap 10y semi-tahunan (nosional 1).
type SwapPV01Report struct {
	Maturity         float64
	ParSwapRate      float64
	Annuity          float64
	AnnuityPV01Per1  float64
	AnnuityPV01Per1M float64
	PayerCurveDV01Ex float64
	PayerCurveDV01FD float64
	PayerCurveDV011M float64
}

// AnalyzeSwapPV01 menghitung tingkat par swap dan PV01/DV01 per 1 bp untuk swap suku bunga par.
func AnalyzeSwapPV01(curve ns.Params, maturity, freq float64) SwapPV01Report {
	n := int(math.Round(maturity * freq))
	dt := 1.0 / freq
	annuity := 0.0
	dAnnDy := 0.0
	for i := 1; i <= n; i++ {
		t := float64(i) * dt
		df := curve.Discount(t)
		annuity += dt * df
		dAnnDy -= dt * t * df
	}
	pMat := curve.Discount(maturity)
	parRate := (1.0 - pMat) / annuity
	dPMatDy := -maturity * pMat
	// Nilai payer swap V(y) = 1 - P(T) - K * Annuity
	// dV/dy = -dP(T)/dy - K * dAnn/dy = T*P(T) + K * sum(dt * t * P(t))
	dVdy := -dPMatDy - parRate*dAnnDy

	evalPayer := func(shift float64) float64 {
		ann := 0.0
		for i := 1; i <= n; i++ {
			t := float64(i) * dt
			ann += dt * math.Exp(-t*(curve.Yield(t)+shift))
		}
		pT := math.Exp(-maturity * (curve.Yield(maturity) + shift))
		return 1.0 - pT - parRate*ann
	}
	dv01FD := (evalPayer(0.5e-4) - evalPayer(-0.5e-4))
	return SwapPV01Report{
		Maturity:         maturity,
		ParSwapRate:      parRate,
		Annuity:          annuity,
		AnnuityPV01Per1:  annuity * 1e-4,
		AnnuityPV01Per1M: annuity * 1e-4 * 1_000_000.0,
		PayerCurveDV01Ex: dVdy * 1e-4,
		PayerCurveDV01FD: dv01FD,
		PayerCurveDV011M: dVdy * 1e-4 * 1_000_000.0,
	}
}

// TwoBondHedgeReport merangkum lindung nilai KRD menggunakan dua obligasi pada dua simpul kunci.
type TwoBondHedgeReport struct {
	TargetName        string
	Hedge1Name        string
	Hedge2Name        string
	TargetNode1       float64
	TargetNode2       float64
	WeightHedge1      float64
	WeightHedge2      float64
	TargetDollarKRDs  []float64
	HedgedDollarKRDs  []float64
	ResidualNode1     float64
	ResidualNode2     float64
	ResidualTotalParD float64
	TargetTotalParD   float64
}

// HedgeTwoKRDs mencari bobot dua obligasi lindung nilai (5y 4% dan 10y 2%) untuk menetralkan
// dua simpul KRD utama (nodeIdx1=3 [5y], nodeIdx2=5 [10y]) dari obligasi target 10y 4%.
func HedgeTwoKRDs(curve ns.Params, nodes []float64, nodeIdx1, nodeIdx2 int, bump float64) (TwoBondHedgeReport, error) {
	tTgt, cfTgt := BondCashFlows(10.0, 0.04, 2.0)
	tH1, cfH1 := BondCashFlows(5.0, 0.04, 2.0)
	tH2, cfH2 := BondCashFlows(10.0, 0.02, 2.0)

	kTgt := ComputeKeyRateDurations(curve, tTgt, cfTgt, nodes, bump)
	kH1 := ComputeKeyRateDurations(curve, tH1, cfH1, nodes, bump)
	kH2 := ComputeKeyRateDurations(curve, tH2, cfH2, nodes, bump)

	mat := [][]float64{
		{kH1.DollarKRDClosed[nodeIdx1], kH2.DollarKRDClosed[nodeIdx1]},
		{kH1.DollarKRDClosed[nodeIdx2], kH2.DollarKRDClosed[nodeIdx2]},
	}
	rhs := []float64{
		-kTgt.DollarKRDClosed[nodeIdx1],
		-kTgt.DollarKRDClosed[nodeIdx2],
	}
	w, err := lag.SolveLU(mat, rhs)
	if err != nil {
		return TwoBondHedgeReport{}, fmt.Errorf("HedgeTwoKRDs solve failed: %w", err)
	}

	hedged := make([]float64, len(nodes))
	sumTgt, sumHedged := 0.0, 0.0
	for i := range nodes {
		hedged[i] = kTgt.DollarKRDClosed[i] + w[0]*kH1.DollarKRDClosed[i] + w[1]*kH2.DollarKRDClosed[i]
		sumTgt += kTgt.DollarKRDClosed[i]
		sumHedged += hedged[i]
	}

	return TwoBondHedgeReport{
		TargetName:        "10y_4pct_bond",
		Hedge1Name:        "5y_4pct_bond",
		Hedge2Name:        "10y_2pct_bond",
		TargetNode1:       nodes[nodeIdx1],
		TargetNode2:       nodes[nodeIdx2],
		WeightHedge1:      w[0],
		WeightHedge2:      w[1],
		TargetDollarKRDs:  kTgt.DollarKRDClosed,
		HedgedDollarKRDs:  hedged,
		ResidualNode1:     hedged[nodeIdx1],
		ResidualNode2:     hedged[nodeIdx2],
		ResidualTotalParD: sumHedged,
		TargetTotalParD:   sumTgt,
	}, nil
}
