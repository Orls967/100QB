package panjer

import (
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/closedform"
)

// Series menyimpan batas bawah dan batas atas Panjer pada kisi langkah h sampai mMax*h.
type Series struct {
	H          float64
	MMax       int
	Low        []float64
	High       []float64
	LowCdfSum  []float64 // via 1 - sum g_k untuk mengukur kehilangan akurasi pada u besar
	HighCdfSum []float64
}

// ComputeBounds menghitung [psi_low(m*h), psi_high(m*h)] untuk distribusi F_L dengan rho = Lambda*Mu/c
// hingga m = ceil(uMax / h) menggunakan rekursi ekor positif (tanpa pengurangan 1 - sum g_k)
// sekaligus menghitung versi penjumlahan CDF g_k untuk uji underflow/pembatalan floating-point.
func ComputeBounds(rho float64, fBarL func(float64) float64, uMax, h float64) Series {
	mMax := int(math.Round(uMax / h))
	if mMax < 1 {
		mMax = 1
	}

	// Hitung array fBar[0..mMax+1]
	fBar := make([]float64, mMax+2)
	for j := 0; j <= mMax+1; j++ {
		fBar[j] = fBarL(float64(j) * h)
	}

	// 1. Pembulatan ke bawah Y^- = h * floor(Y/h) -> batas bawah psi_low.
	fLow := make([]float64, mMax+1)
	for j := 0; j <= mMax; j++ {
		diff := fBar[j] - fBar[j+1]
		if diff < 0 {
			diff = 0
		}
		fLow[j] = diff
	}

	// 2. Pembulatan ke atas Y^+ = h * ceil(Y/h) -> batas atas psi_high (fHigh[0] = 0).
	fHigh := make([]float64, mMax+1)
	for j := 1; j <= mMax; j++ {
		diff := fBar[j-1] - fBar[j]
		if diff < 0 {
			diff = 0
		}
		fHigh[j] = diff
	}

	psiLow := make([]float64, mMax+1)
	psiHigh := make([]float64, mMax+1)
	scaleLow := rho / (1.0 - rho*fLow[0])

	for m := 0; m <= mMax; m++ {
		convL := 0.0
		for j := 1; j <= m; j++ {
			convL += fLow[j] * psiLow[m-j]
		}
		psiLow[m] = scaleLow * (fBar[m+1] + convL)

		convH := 0.0
		for j := 1; j <= m; j++ {
			convH += fHigh[j] * psiHigh[m-j]
		}
		psiHigh[m] = rho * (fBar[m] + convH)
	}

	// Versi alternatif via PMF majemuk geometrik g_m lalu psi = 1 - sum g_k.
	gLow := make([]float64, mMax+1)
	gHigh := make([]float64, mMax+1)
	lowCdf := make([]float64, mMax+1)
	highCdf := make([]float64, mMax+1)

	gLow[0] = (1.0 - rho) / (1.0 - rho*fLow[0])
	sumL := gLow[0]
	lowCdf[0] = 1.0 - sumL

	gHigh[0] = 1.0 - rho
	sumH := gHigh[0]
	highCdf[0] = 1.0 - sumH

	for m := 1; m <= mMax; m++ {
		accL := 0.0
		accH := 0.0
		for j := 1; j <= m; j++ {
			accL += fLow[j] * gLow[m-j]
			accH += fHigh[j] * gHigh[m-j]
		}
		gLow[m] = scaleLow * accL
		sumL += gLow[m]
		lowCdf[m] = 1.0 - sumL

		gHigh[m] = rho * accH
		sumH += gHigh[m]
		highCdf[m] = 1.0 - sumH
	}

	return Series{
		H:          h,
		MMax:       mMax,
		Low:        psiLow,
		High:       psiHigh,
		LowCdfSum:  lowCdf,
		HighCdfSum: highCdf,
	}
}

// At mengembalikan batas bawah, titik tengah (referensi), batas atas, dan lebar selang pada modal u = m*h.
func (s Series) At(u float64) (low, mid, high, width, lowCdf, highCdf float64) {
	m := int(math.Round(u / s.H))
	if m < 0 {
		m = 0
	}
	if m > s.MMax {
		m = s.MMax
	}
	low = s.Low[m]
	high = s.High[m]
	mid = 0.5 * (low + high)
	width = high - low
	return low, mid, high, width, s.LowCdfSum[m], s.HighCdfSum[m]
}

// BoundRow merangkum perbandingan batas Panjer terhadap nilai eksak/asimtotik di satu (dist, h, u).
type BoundRow struct {
	DistCode  string
	DistName  string
	H         float64
	U         float64
	PsiLow    float64
	PsiMid    float64
	PsiHigh   float64
	Width     float64
	LowCdfSum float64
	PsiExact  float64
	Bracketed bool
	MidAbsErr float64
	Asymptote float64
	RatioAsym float64
}

// EvaluateGrid menghitung seluruh baris tabel Panjer untuk satu distribusi melintasi hList dan uList.
func EvaluateGrid(d claims.Dist, hList, uList []float64) []BoundRow {
	uMax := 0.0
	for _, u := range uList {
		if u > uMax {
			uMax = u
		}
	}
	cf, errCF := closedform.Derive(d)
	hasCF := errCF == nil
	rows := make([]BoundRow, 0, len(hList)*len(uList))

	for _, h := range hList {
		ser := ComputeBounds(claims.Rho, d.LadderSurvival, uMax, h)
		for _, u := range uList {
			lo, mid, hi, width, loCdf, _ := ser.At(u)
			psiEx := math.NaN()
			errMid := math.NaN()
			bracket := false
			if hasCF {
				psiEx = cf.Eval(u)
				bracket = lo <= psiEx+1e-15 && psiEx <= hi+1e-15
				errMid = math.Abs(mid - psiEx)
			}
			asym := d.SubexpAsymptotic(u)
			ratio := 0.0
			if asym > 0 {
				ratio = mid / asym
			}
			rows = append(rows, BoundRow{
				DistCode:  d.Code,
				DistName:  d.Name,
				H:         h,
				U:         u,
				PsiLow:    lo,
				PsiMid:    mid,
				PsiHigh:   hi,
				Width:     width,
				LowCdfSum: loCdf,
				PsiExact:  psiEx,
				Bracketed: bracket,
				MidAbsErr: errMid,
				Asymptote: asym,
				RatioAsym: ratio,
			})
		}
	}
	return rows
}
