package process

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
)

func TestProcessMomentsDeficitAndWorkerDeterminism(t *testing.T) {
	mom := VerifyMoments(claims.AllDists(), 10.0, 10.0, 40000, 8, 4, rng.MainSeed)
	for _, r := range mom {
		if r.ZN > 3.5 || r.ZUt > 3.5 {
			t.Fatalf("%s moment z too high: zN=%v zUt=%v", r.DistName, r.ZN, r.ZUt)
		}
		// Pareto(alpha=2.5) memiliki E[X^4] = inf, sehingga varians sampel dari Var(S_t) tidak berhingga (konvergensi lambat ber-bias bawah).
		if r.DistCode != "D" && math.Abs(r.VarStDiff) > 0.05 {
			t.Fatalf("%s Var(S_t) rel diff=%v too large", r.DistName, r.VarStDiff)
		}
		t.Logf("%s: zN=%.2f zUt=%.2f VarSt sim=%.3f theo=%.3f (relDiff=%.4f)",
			r.DistName, r.ZN, r.ZUt, r.SimVarSt, r.TheoVarSt, r.VarStDiff)
	}

	defs := VerifyExponentialDeficit([]float64{2.0, 5.0}, 100000, 1000, 8, 4, rng.MainSeed)
	for _, df := range defs {
		meanZ := math.Abs(df.DeficitMean-1.0) / df.DeficitMeanSE
		t.Logf("Deficit u=%.0f: psiMC=%.5f (z=%.2f) mean=%.5f+-%.5f (z=%.2f) SD=%.5f KS=%.5f crit95=%.5f p=%.3f",
			df.U, df.PsiMC, df.PsiZScore, df.DeficitMean, df.DeficitMeanSE, meanZ, df.DeficitSD, df.KSStat, df.KSCrit95, df.KSPValue)
		if df.PsiZScore > 3.5 || meanZ > 3.5 || !df.KSPass {
			t.Fatalf("Exponential deficit test failed at u=%g: %+v", df.U, df)
		}
	}

	mc1 := SimulateCrudeMC(claims.LightTailedDists(), []float64{0, 5, 10}, 25000, 8, 1, rng.MainSeed)
	mc4 := SimulateCrudeMC(claims.LightTailedDists(), []float64{0, 5, 10}, 25000, 8, 4, rng.MainSeed)
	for i := range mc1 {
		if math.Float64bits(mc1[i].PsiInfMC) != math.Float64bits(mc4[i].PsiInfMC) {
			t.Fatalf("CrudeMC worker determinism failed at idx %d", i)
		}
		if mc1[i].ZScore > 3.5 {
			t.Fatalf("%s u=%g CrudeMC z=%v > 3.5", mc1[i].DistName, mc1[i].U, mc1[i].ZScore)
		}
		if !(mc1[i].PsiT50 <= mc1[i].PsiT200 && mc1[i].PsiT200 <= mc1[i].PsiInfMC) {
			t.Fatalf("%s u=%g horizon monotonicity violated: %v <= %v <= %v", mc1[i].DistName, mc1[i].U, mc1[i].PsiT50, mc1[i].PsiT200, mc1[i].PsiInfMC)
		}
	}
}
