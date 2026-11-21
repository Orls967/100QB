package is

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/lundberg"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
)

func TestImportanceSamplingAndMisspecifiedR(t *testing.T) {
	for _, d := range claims.LightTailedDists() {
		resL, err := lundberg.Solve(d)
		if err != nil {
			t.Fatal(err)
		}
		row100_1 := Simulate(d, 100.0, resL.R, 1.0, 20000, 8, 1, rng.MainSeed)
		row100_4 := Simulate(d, 100.0, resL.R, 1.0, 20000, 8, 4, rng.MainSeed)
		if math.Float64bits(row100_1.PsiIS) != math.Float64bits(row100_4.PsiIS) {
			t.Fatalf("IS worker invariance failed for %s", d.Name)
		}
		if row100_1.ZScore > 3.5 || row100_1.RelSEIS > 0.05 {
			t.Fatalf("IS failed at u=100 for %s: z=%v relSE=%v psi=%v exact=%v", d.Name, row100_1.ZScore, row100_1.RelSEIS, row100_1.PsiIS, row100_1.PsiExact)
		}

		// Eksponen salah 0.7R dan 1.3R: tetap tak-bias, variansi lebih besar dari R optimal.
		r07 := Simulate(d, 15.0, 0.7*resL.R, 0.7, 25000, 8, 4, rng.MainSeed)
		r10 := Simulate(d, 15.0, 1.0*resL.R, 1.0, 25000, 8, 4, rng.MainSeed)
		r13 := Simulate(d, 15.0, 1.3*resL.R, 1.3, 25000, 8, 4, rng.MainSeed)
		if r07.ZScore > 3.5 || r10.ZScore > 3.5 || r13.ZScore > 3.5 {
			t.Fatalf("%s tilted unbiased check failed: z=(%v,%v,%v)", d.Name, r07.ZScore, r10.ZScore, r13.ZScore)
		}
		if r10.SDIS >= r07.SDIS || r10.SDIS >= r13.SDIS {
			t.Fatalf("%s optimal R expected min variance: sd=(%v, %v, %v)", d.Name, r07.SDIS, r10.SDIS, r13.SDIS)
		}
	}

	div := VerifyHeavyTailMGFFailure(claims.HeavyTailedDists(), 0.1, []float64{10, 50, 200})
	if len(div) != 6 || div[2].PartialMGF <= div[0].PartialMGF || div[5].PartialMGF <= div[3].PartialMGF {
		t.Fatalf("heavy-tail partial MGF divergence check failed: %+v", div)
	}
}
