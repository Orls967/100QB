package kalman

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/afns"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/optim"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
)

func TestJointGaussianAnchorExactT5(t *testing.T) {
	spec := DefaultSpec()
	taus := []float64{0.5, 2.0, 10.0}
	rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xA001, 0))
	_, yields := afns.SimulateStationaryPath(spec.AFNS, taus, 5, WeeklyDt, spec.SigmaEps, rg)

	anc, err := VerifyJointGaussianAnchor(spec, taus, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("VerifyJointGaussianAnchor error: %v", err)
	}
	if anc.LogLikAbsDiff >= 1e-9 || anc.MaxSmoothMeanAbsDiff >= 1e-9 || anc.MaxSmoothCovAbsDiff >= 1e-9 {
		t.Fatalf("anchor mismatch: ll=%g mean=%g cov=%g", anc.LogLikAbsDiff, anc.MaxSmoothMeanAbsDiff, anc.MaxSmoothCovAbsDiff)
	}
	fastLL, err := FastLogLikelihood(spec, taus, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("FastLogLikelihood error: %v", err)
	}
	if math.Abs(fastLL-anc.JointLogLik) >= 1e-9 {
		t.Fatalf("FastLogLikelihood vs Joint mismatch: %g", math.Abs(fastLL-anc.JointLogLik))
	}
}

func TestFastLogLikelihoodMatchesFullFilterWithRiccatiSteadyState(t *testing.T) {
	spec := DefaultSpec()
	rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xB002, 1040))
	_, yields := afns.SimulateStationaryPath(spec.AFNS, StandardTenors, 1040, WeeklyDt, spec.SigmaEps, rg)

	fastLL, err := FastLogLikelihood(spec, StandardTenors, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("FastLogLikelihood error: %v", err)
	}
	fullOut, err := FilterAndSmooth(spec, StandardTenors, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("FilterAndSmooth error: %v", err)
	}
	if diff := math.Abs(fastLL - fullOut.LogLikelihood); diff >= 1e-9 {
		t.Fatalf("Riccati steady-state diverged from FilterAndSmooth: diff=%g", diff)
	}
}

func TestScaledScoreGradientVsRichardsonAndNonFiniteHandling(t *testing.T) {
	trueSpec := DefaultSpec()
	rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xC003, 200))
	_, yields := afns.SimulateStationaryPath(trueSpec.AFNS, StandardTenors, 200, WeeklyDt, trueSpec.SigmaEps, rg)
	initSpec := PerturbInitialGuess(trueSpec, rg)

	negLL := func(u []float64) float64 {
		ll, err := FastLogLikelihood(VecToSpec(UnpackUnconstrained(u), true), StandardTenors, WeeklyDt, yields)
		if err != nil {
			return 1e30
		}
		return -ll
	}
	_, _, errInit := optim.VerifyGradVsRichardson(negLL, PackUnconstrained(SpecToVec(initSpec)), nil)
	_, _, errTrue := optim.VerifyGradVsRichardson(negLL, PackUnconstrained(SpecToVec(trueSpec)), nil)
	if errInit >= 1e-5 || errTrue >= 1e-5 {
		t.Fatalf("score vs Richardson too large: init=%g true=%g", errInit, errTrue)
	}

	// Uji penanganan NaN pada data observasi -> MLE harus menandai Finite=false & TerminationReason=invalid_objective_or_gradient
	badYields := make([][]float64, len(yields))
	for i := range yields {
		badYields[i] = append([]float64(nil), yields[i]...)
	}
	badYields[0][0] = math.NaN()
	badRes := EstimateMLEWithOptions(initSpec, trueSpec, StandardTenors, WeeklyDt, badYields, false)
	if badRes.Finite || badRes.Converged || badRes.TerminationReason != optim.TermInvalidObjective {
		t.Fatalf("expected Finite=false invalid_objective_or_gradient, got Finite=%v Term=%s", badRes.Finite, badRes.TerminationReason)
	}
}
