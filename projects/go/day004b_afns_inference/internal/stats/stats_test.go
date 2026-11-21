package stats

import (
	"math"
	"testing"
)

func TestWilsonCIAndLinearTrend(t *testing.T) {
	// Edge case: 50/50 -> Wilson CI harus ketat di dalam (0, 1] dan lower < 1.0
	p, lo, hi := WilsonCI(50, 50, Z95)
	if p != 1.0 || !(lo > 0.90 && lo < 1.0) || math.Abs(hi-1.0) > 1e-12 {
		t.Fatalf("unexpected WilsonCI(50,50): p=%v lo=%v hi=%v", p, lo, hi)
	}
	// Edge case: 0/50 -> lower = 0, upper < 0.08
	p0, lo0, hi0 := WilsonCI(0, 50, Z95)
	if p0 != 0.0 || math.Abs(lo0) > 1e-12 || !(hi0 > 0.05 && hi0 < 0.09) {
		t.Fatalf("unexpected WilsonCI(0,50): p=%v lo=%v hi=%v", p0, lo0, hi0)
	}
	// Known example: 34/50 = 0.68 -> CI sekitar [0.542, 0.792]
	p34, lo34, hi34 := WilsonCI(34, 50, Z95)
	if math.Abs(p34-0.68) > 1e-12 || lo34 >= 0.68 || hi34 <= 0.68 || math.Abs(lo34-0.5419) > 1e-3 || math.Abs(hi34-0.7924) > 1e-3 {
		t.Fatalf("unexpected WilsonCI(34,50): p=%v lo=%v hi=%v", p34, lo34, hi34)
	}

	// Exact linear trend test: y = 0.02 + 4.0 * x
	xs := []float64{0.1, 0.05, 0.025, 0.0125}
	ys := make([]float64, len(xs))
	ses := []float64{0.02, 0.01, 0.005, 0.0025}
	for i, x := range xs {
		ys[i] = 0.02 + 4.0*x
	}
	fit := FitLinearTrend(xs, ys, ses)
	if math.Abs(fit.OLSSlope-4.0) > 1e-12 || math.Abs(fit.OLSIntercept-0.02) > 1e-12 {
		t.Fatalf("OLS failed exact line: slope=%v intercept=%v", fit.OLSSlope, fit.OLSIntercept)
	}
	if !fit.WLSWeightsAvailable || math.Abs(fit.WLSSlope-4.0) > 1e-12 || math.Abs(fit.WLSIntercept-0.02) > 1e-12 {
		t.Fatalf("WLS failed exact line: slope=%v intercept=%v", fit.WLSSlope, fit.WLSIntercept)
	}
	// Verify linear combination property: sum_i w_i * y_i == Slope
	chkOLS, chkWLS := 0.0, 0.0
	for i := range ys {
		chkOLS += fit.OLSWeightsX[i] * ys[i]
		chkWLS += fit.WLSWeightsX[i] * ys[i]
	}
	if math.Abs(chkOLS-fit.OLSSlope) > 1e-12 || math.Abs(chkWLS-fit.WLSSlope) > 1e-12 {
		t.Fatalf("weight linear combination mismatch: %v vs %v, %v vs %v", chkOLS, fit.OLSSlope, chkWLS, fit.WLSSlope)
	}
}
