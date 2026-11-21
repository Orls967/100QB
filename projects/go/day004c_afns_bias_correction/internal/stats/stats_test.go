package stats

import (
	"math"
	"testing"
)

func TestWilsonCIQuantileAndInvalidDenominator(t *testing.T) {
	p, lo, hi := WilsonCI(50, 50, Z95)
	if p != 1.0 || !(lo > 0.92 && lo < 0.94) || math.Abs(hi-1.0) > 1e-12 {
		t.Fatalf("unexpected WilsonCI(50,50): p=%v lo=%v hi=%v", p, lo, hi)
	}
	p0, lo0, hi0 := WilsonCI(0, 50, Z95)
	if p0 != 0.0 || math.Abs(lo0) > 1e-12 || !(hi0 > 0.06 && hi0 < 0.08) {
		t.Fatalf("unexpected WilsonCI(0,50): p=%v lo=%v hi=%v", p0, lo0, hi0)
	}

	// Uji cakupan dengan replikasi luar tidak valid tetap masuk ke denominator
	valid := []bool{true, true, false, true, false}
	covered := []bool{true, false, true, true, true} // indeks 2 & 4 tidak valid walau bertanda true -> wajib dihitung tidak tercakup!
	vc, ic, cc, tot, rate, wlo, whi := CoverageWithValidity(valid, covered)
	if vc != 3 || ic != 2 || cc != 2 || tot != 5 || math.Abs(rate-0.4) > 1e-12 || !(wlo < 0.4 && whi > 0.4) {
		t.Fatalf("unexpected CoverageWithValidity: vc=%d ic=%d cc=%d tot=%d rate=%v [%v,%v]", vc, ic, cc, tot, rate, wlo, whi)
	}

	xs := []float64{5, 1, 4, 2, 3}
	if q50 := Quantile(xs, 0.5); math.Abs(q50-3.0) > 1e-12 {
		t.Fatalf("unexpected median: %v", q50)
	}
	if qEmpty := Quantile(nil, 0.5); !math.IsNaN(qEmpty) {
		t.Fatalf("expected NaN quantile for empty slice, got %v", qEmpty)
	}
}
