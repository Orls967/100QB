package bootstrap

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
)

// SemiAnnualGrid mengembalikan grid maturitas semi-tahunan 0.5, 1.0, ..., maxYear.
func SemiAnnualGrid(maxYear float64) []float64 {
	n := int(math.Round(maxYear * 2.0))
	out := make([]float64, n)
	for i := 0; i < n; i++ {
		out[i] = 0.5 * float64(i+1)
	}
	return out
}

// ParRatesFromCurve menghitung tingkat kupon par obligasi semi-tahunan pada grid maturitas seragam
// c_k = (1 - P(tau_k)) / (delta * sum_{j=1}^k P(tau_j)).
func ParRatesFromCurve(curve ns.Params, maturities []float64, delta float64) (parRates, trueDF, trueYields []float64) {
	n := len(maturities)
	parRates = make([]float64, n)
	trueDF = make([]float64, n)
	trueYields = make([]float64, n)
	annuity := 0.0
	for i, tau := range maturities {
		y := curve.Yield(tau)
		df := math.Exp(-tau * y)
		annuity += delta * df
		trueYields[i] = y
		trueDF[i] = df
		parRates[i] = (1.0 - df) / annuity
	}
	return parRates, trueDF, trueYields
}

// BootstrapFromParRates memulihkan faktor diskonto P(tau_k) dan zero yield y(tau_k) secara berurutan
// dari tingkat par pada grid berkala dengan interval delta.
func BootstrapFromParRates(maturities, parRates []float64, delta float64) (dfs, zeroYields []float64, err error) {
	n := len(maturities)
	if n == 0 || len(parRates) != n {
		return nil, nil, fmt.Errorf("BootstrapFromParRates: invalid slice lengths")
	}
	dfs = make([]float64, n)
	zeroYields = make([]float64, n)
	prevAnnuity := 0.0
	for i := 0; i < n; i++ {
		c := parRates[i]
		denom := 1.0 + delta*c
		if denom <= 0 {
			return nil, nil, fmt.Errorf("BootstrapFromParRates: non-positive denominator at k=%d", i)
		}
		p := (1.0 - c*prevAnnuity) / denom
		if p <= 0 {
			return nil, nil, fmt.Errorf("BootstrapFromParRates: non-positive discount factor %g at k=%d", p, i)
		}
		dfs[i] = p
		zeroYields[i] = -math.Log(p) / maturities[i]
		prevAnnuity += delta * p
	}
	return dfs, zeroYields, nil
}
