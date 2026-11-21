package ak

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/panjer"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
)

func TestAsmussenKroeseHeavyTailsAndWorkers(t *testing.T) {
	for _, d := range claims.HeavyTailedDists() {
		pSer := panjer.ComputeBounds(claims.Rho, d.LadderSurvival, 200.0, 0.05)
		for _, u := range []float64{10, 50, 200} {
			r1 := Simulate(d, u, pSer, 25000, 8, 1, rng.MainSeed)
			r4 := Simulate(d, u, pSer, 25000, 8, 4, rng.MainSeed)
			if math.Float64bits(r1.PsiAK) != math.Float64bits(r4.PsiAK) {
				t.Fatalf("AK worker invariance failed for %s u=%g", d.Name, u)
			}
			if !r1.WithinBand {
				t.Fatalf("AK vs Panjer failed for %s u=%g: ak=%v se=%v band=[%v,%v]",
					d.Name, u, r1.PsiAK, r1.SEAK, r1.PanjerLow, r1.PanjerHigh)
			}
			t.Logf("%s u=%.0f: ak=%.5g (relSE=%.3f) vs panjer=[%.5g, %.5g], ak/asym=%.3f",
				d.Name, u, r1.PsiAK, r1.RelSEAK, r1.PanjerLow, r1.PanjerHigh, r1.AKOverAsym)
		}
	}
}
