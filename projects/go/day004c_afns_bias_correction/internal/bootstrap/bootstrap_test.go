package bootstrap

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/afns"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/stats"
)

func TestEvalParamBootstrapBiasQuantilesAndInvalidHandling(t *testing.T) {
	// Kasus normal: 100 sampel seragam/bergeser dengan bias +0.20 relatif terhadap outerEst=0.25
	draws := make([]float64, 100)
	for i := 0; i < 100; i++ {
		draws[i] = 0.05 + 0.004*float64(i) + 0.20 // range [0.25, 0.646], mean = 0.448
	}
	res := EvalParamBootstrap(1, 0.25, 0.05, draws, 100)
	if !res.IntervalValid || res.FiniteCount != 100 || res.InvalidCount != 0 {
		t.Fatalf("expected valid bootstrap summary, got %+v", res)
	}
	// BootBias = 0.448 - 0.25 = 0.198 -> BootCorrected = 0.25 - 0.198 = 0.052
	if math.Abs(res.BootMean-0.448) > 1e-12 || math.Abs(res.BootBias-0.198) > 1e-12 || math.Abs(res.BootCorrected-0.052) > 1e-12 {
		t.Fatalf("unexpected boot bias correction: %+v", res)
	}
	if res.Percentile95Covered {
		t.Fatalf("expected trueVal=0.05 outside [%.4f, %.4f]", res.Percentile95Lower, res.Percentile95Upper)
	}

	// Kasus invalid (>5% NaN -> 90 finite < 95): walaupun mencakup trueVal secara numerik,
	// IntervalValid=false dan Percentile95Covered=false serta outer rep tetap di denominator cakupan!
	drawsBad := make([]float64, 100)
	for i := 0; i < 90; i++ {
		drawsBad[i] = 0.04 + 0.0002*float64(i) // mencakup 0.05
	}
	for i := 90; i < 100; i++ {
		drawsBad[i] = math.NaN()
	}
	badRes := EvalParamBootstrap(1, 0.05, 0.05, drawsBad, 100)
	if badRes.IntervalValid || badRes.Percentile95Covered || badRes.FiniteCount != 90 || badRes.InvalidCount != 10 {
		t.Fatalf("expected invalid bootstrap interval with covered=false, got %+v", badRes)
	}

	// Cek agregasi cakupan dengan 1 invalid rep + 1 covered rep -> rate = 1/2 = 0.50 (bukan 1/1)
	_, invCnt, covCnt, tot, rate, _, _ := stats.CoverageWithValidity(
		[]bool{res.IntervalValid, badRes.IntervalValid, true},
		[]bool{true, true, false},
	)
	if invCnt != 1 || covCnt != 1 || tot != 3 || math.Abs(rate-1.0/3.0) > 1e-12 {
		t.Fatalf("unexpected coverage denominator with invalid outer rep: inv=%d cov=%d tot=%d rate=%v", invCnt, covCnt, tot, rate)
	}
}

func TestBootstrapStudyWorkerDeterminism(t *testing.T) {
	spec := kalman.DefaultSpec()
	var outerFits [2]kalman.MLEResult
	for r := 0; r < 2; r++ {
		rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xBEEF, uint64(r)))
		_, y := afns.SimulateStationaryPath(spec.AFNS, kalman.StandardTenors, 80, kalman.WeeklyDt, spec.SigmaEps, rg)
		init := kalman.PerturbInitialGuess(spec, rg)
		outerFits[r] = kalman.EstimateMLE(init, spec, kalman.StandardTenors, kalman.WeeklyDt, y)
	}

	b1 := RunParametricBootstrapStudy(spec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, outerFits[:], 3, 1, rng.MainSeed, 0xB100)
	b4 := RunParametricBootstrapStudy(spec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, outerFits[:], 3, 4, rng.MainSeed, 0xB100)

	for i := 0; i < 11; i++ {
		if math.Float64bits(b1.Params[i].BootCorrectedMean) != math.Float64bits(b4.Params[i].BootCorrectedMean) {
			t.Fatalf("param %d BootCorrectedMean diverged across 1 vs 4 workers", i)
		}
	}
	for r := 0; r < 2; r++ {
		for b := 0; b < 3; b++ {
			if math.Float64bits(b1.OuterReplications[r].InnerDraws[b].EstLogLik) != math.Float64bits(b4.OuterReplications[r].InnerDraws[b].EstLogLik) {
				t.Fatalf("outer %d boot %d EstLogLik diverged across 1 vs 4 workers", r, b)
			}
		}
	}
}
