package interp

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
)

func TestInterpolatorsExactKnotsAndWorkerDeterminism(t *testing.T) {
	curve := ns.DefaultNS()
	mats := bootstrap.SemiAnnualGrid(30.0)
	_, trueDF, trueY := bootstrap.ParRatesFromCurve(curve, mats, 0.5)

	ll, err := NewLogLinearP(mats, trueDF)
	if err != nil {
		t.Fatalf("NewLogLinearP failed: %v", err)
	}
	sp, err := NewNaturalCubicSplineY(mats, trueY)
	if err != nil {
		t.Fatalf("NewNaturalCubicSplineY failed: %v", err)
	}
	for i, tau := range mats {
		yll, _ := ll.Eval(tau)
		ysp, _ := sp.Eval(tau)
		if math.Abs(yll-trueY[i]) > 1e-13 {
			t.Fatalf("LogLinearP knot mismatch at %g: %v vs %v", tau, yll, trueY[i])
		}
		if math.Abs(ysp-trueY[i]) > 1e-13 {
			t.Fatalf("NaturalCubicSplineY knot mismatch at %g: %v vs %v", tau, ysp, trueY[i])
		}
	}

	// Uji determinisme 1 vs 4 pekerja
	res1 := RunNoisyBootstrapComparison(curve, 0.5e-4, 40, 8, 1, rng.MainSeed)
	res4 := RunNoisyBootstrapComparison(curve, 0.5e-4, 40, 8, 4, rng.MainSeed)
	if math.Float64bits(res1.LogLinearP.ForwardTVMean) != math.Float64bits(res4.LogLinearP.ForwardTVMean) ||
		math.Float64bits(res1.CubicSpline.ForwardTVMean) != math.Float64bits(res4.CubicSpline.ForwardTVMean) ||
		math.Float64bits(res1.CubicSpline.YieldRMSEMean) != math.Float64bits(res4.CubicSpline.YieldRMSEMean) {
		t.Fatalf("Worker invariance violated in noisy bootstrap comparison")
	}
	t.Logf("Noisy R=40 LogLin TV=%.4f yRMSE=%.3e | Spline TV=%.4f yRMSE=%.3e",
		res1.LogLinearP.ForwardTVMean, res1.LogLinearP.YieldRMSEMean,
		res1.CubicSpline.ForwardTVMean, res1.CubicSpline.YieldRMSEMean)
}
