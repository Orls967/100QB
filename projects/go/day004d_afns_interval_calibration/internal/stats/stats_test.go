package stats

import (
	"math"
	"testing"
)

func TestWilsonCIQuantileAndInvalidDenominator(t *testing.T) {
	p, lo, hi := WilsonCI(95, 100, Z95)
	if math.Abs(p-0.95) > 1e-12 || lo <= 0.88 || hi >= 0.99 {
		t.Fatalf("unexpected WilsonCI(95,100): p=%g lo=%g hi=%g", p, lo, hi)
	}
	qs := []float64{10, 20, 30, 40, 50}
	if q := Quantile(qs, 0.25); math.Abs(q-20.0) > 1e-12 {
		t.Fatalf("unexpected Quantile(0.25): %g", q)
	}
	vCnt, invCnt, covCnt, tot, rate, _, _ := CoverageWithValidity(
		[]bool{true, false, true, true},
		[]bool{true, true, false, true},
	)
	if vCnt != 3 || invCnt != 1 || covCnt != 2 || tot != 4 || math.Abs(rate-0.5) > 1e-12 {
		t.Fatalf("unexpected CoverageWithValidity: v=%d inv=%d cov=%d tot=%d rate=%g", vCnt, invCnt, covCnt, tot, rate)
	}
}
