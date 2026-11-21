package kalman

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/afns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
)

func TestJointGaussianT5AnchorSmootherAndFastLik(t *testing.T) {
	spec := DefaultSpec()
	shortTenors := []float64{1.0, 5.0, 10.0}
	rg := rng.New(rng.DeriveChunkSeed(rng.MainSeed, 0x3101, 0))
	_, yields := afns.SimulateStationaryPath(spec.AFNS, shortTenors, 5, WeeklyDt, spec.SigmaEps, rg)

	anchor, err := VerifyJointGaussianAnchor(spec, shortTenors, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("VerifyJointGaussianAnchor failed: %v", err)
	}
	t.Logf("T=5 Joint anchor: llDiff=%.3e smoothMeanDiff=%.3e smoothCovDiff=%.3e",
		anchor.LogLikAbsDiff, anchor.MaxSmoothMeanAbsDiff, anchor.MaxSmoothCovAbsDiff)
	if anchor.LogLikAbsDiff >= 1e-9 || anchor.MaxSmoothMeanAbsDiff >= 1e-9 || anchor.MaxSmoothCovAbsDiff >= 1e-9 {
		t.Fatalf("Joint anchor tolerance 1e-9 violated: %+v", anchor)
	}

	// Uji FastLogLikelihood vs FilterAndSmooth pada 10 tenor T=200
	_, yFull := afns.SimulateStationaryPath(spec.AFNS, StandardTenors, 200, WeeklyDt, spec.SigmaEps, rg)
	fastLL, err := FastLogLikelihood(spec, StandardTenors, WeeklyDt, yFull)
	if err != nil {
		t.Fatalf("FastLogLikelihood failed: %v", err)
	}
	fullOut, err := FilterAndSmooth(spec, StandardTenors, WeeklyDt, yFull)
	if err != nil {
		t.Fatalf("FilterAndSmooth failed: %v", err)
	}
	if d := math.Abs(fastLL - fullOut.LogLikelihood); d >= 1e-11 {
		t.Fatalf("FastLogLikelihood vs FilterAndSmooth mismatch: %g", d)
	}
}

func TestSingleMLEAndWorkerDeterminism(t *testing.T) {
	spec := DefaultSpec()
	rg := rng.New(rng.DeriveChunkSeed(rng.MainSeed, 0x3201, 0))
	_, yields := afns.SimulateStationaryPath(spec.AFNS, StandardTenors, 1000, WeeklyDt, spec.SigmaEps, rg)
	initSpec := PerturbInitialGuess(spec, rg)

	mle := EstimateMLE(initSpec, spec, StandardTenors, WeeklyDt, yields)
	t.Logf("Single MLE: llInit=%.2f llTrue=%.2f llEst=%.2f gradErrInit=%.2e gradErrTrue=%.2e gradAbsOpt=%.2e minEig=%.2f conv=%v",
		mle.InitLogLik, mle.TrueLogLik, mle.EstLogLik, mle.GradRelErrInit, mle.GradRelErrTrue, mle.GradAbsDiffOpt, mle.HessianMinEig, mle.Converged)
	for i, name := range ParamNames {
		t.Logf("  %-9s true=%.6f est=%.6f se=%.6f", name, mle.TrueParams[i], mle.EstParams[i], mle.HessianSE[i])
	}
	if mle.GradRelErrInit >= 1e-5 || mle.GradRelErrTrue >= 1e-5 || mle.GradAbsDiffOpt >= 1e-4 {
		t.Fatalf("Score FD vs Richardson error too large: init=%g true=%g absOpt=%g", mle.GradRelErrInit, mle.GradRelErrTrue, mle.GradAbsDiffOpt)
	}
	if mle.EstLogLik < mle.TrueLogLik-1e-6 {
		t.Fatalf("MLE logLik < true logLik: est=%.6f true=%.6f", mle.EstLogLik, mle.TrueLogLik)
	}

	// Uji determinisme studi replikasi kecil 1 vs 4 pekerja
	rep1 := RunReplicationStudy(spec, StandardTenors, WeeklyDt, 300, 4, 4, 1, rng.MainSeed)
	rep4 := RunReplicationStudy(spec, StandardTenors, WeeklyDt, 300, 4, 4, 4, rng.MainSeed)
	for i := 0; i < 11; i++ {
		if math.Float64bits(rep1.Params[i].MeanEstimate) != math.Float64bits(rep4.Params[i].MeanEstimate) {
			t.Fatalf("Replication worker invariance violated at param %d", i)
		}
	}
}
