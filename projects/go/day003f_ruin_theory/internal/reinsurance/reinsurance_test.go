package reinsurance

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
)

func TestReinsuranceClosedFormOptimalAndCapital(t *testing.T) {
	expD := claims.Exponential()
	optA, err := FindOptimalRetention(expD, 10.0, 0.02)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(optA.MMin-math.Log(2.0)) > 1e-12 {
		t.Fatalf("Exp M_min=%v != ln(2)", optA.MMin)
	}
	if !optA.MaxClosedVsBand {
		t.Fatalf("Exp closed-form retention formula outside Panjer band")
	}
	if optA.MinPsiMid >= optA.PsiNoReins {
		t.Fatalf("Expected optimal reinsurance to beat no-reins at u=10: min=%v vs noReins=%v", optA.MinPsiMid, optA.PsiNoReins)
	}
	if optA.InfLimitDiff > 1e-10 {
		t.Fatalf("M->inf limit diff=%v", optA.InfLimitDiff)
	}

	optC, err := FindOptimalRetention(claims.Hyperexp2(), 10.0, 0.02)
	if err != nil {
		t.Fatal(err)
	}
	if optC.MinPsiMid >= optC.PsiNoReins {
		t.Fatalf("Expected optimal reinsurance to beat no-reins for Hyperexp2: min=%v vs noReins=%v", optC.MinPsiMid, optC.PsiNoReins)
	}

	capA, err := SolveMinimumCapital(expD, optA.OptimalM, []float64{1e-2, 1e-3, 1e-4}, 0.05)
	if err != nil {
		t.Fatal(err)
	}
	for _, row := range capA {
		if row.UStarWithReins >= row.UStarNoReins {
			t.Fatalf("expected capital saving at eps=%g: with=%v vs without=%v", row.Epsilon, row.UStarWithReins, row.UStarNoReins)
		}
		t.Logf("Exp eps=%.0e: uNo=%.3f uWith(M=%.2f)=%.3f (hemat %.2f%%)",
			row.Epsilon, row.UStarNoReins, row.RetentionM, row.UStarWithReins, row.SavingPercent)
	}
}
