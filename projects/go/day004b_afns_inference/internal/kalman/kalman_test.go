package kalman

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/afns"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/optim"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
)

func TestJointGaussianAnchorExactT5(t *testing.T) {
	spec := DefaultSpec()
	taus := []float64{0.5, 2.0, 10.0}
	rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xA001, 0))
	_, yields := afns.SimulateStationaryPath(spec.AFNS, taus, 5, WeeklyDt, spec.SigmaEps, rg)

	anc, err := VerifyJointGaussianAnchor(spec, taus, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("VerifyJointGaussianAnchor failed: %v", err)
	}
	if anc.LogLikAbsDiff >= 1e-9 {
		t.Fatalf("Kalman vs Joint log-likelihood mismatch: %g (>= 1e-9)", anc.LogLikAbsDiff)
	}
	if anc.MaxSmoothMeanAbsDiff >= 1e-9 {
		t.Fatalf("RTS smooth mean vs Joint conditional mismatch: %g (>= 1e-9)", anc.MaxSmoothMeanAbsDiff)
	}
	if anc.MaxSmoothCovAbsDiff >= 1e-9 {
		t.Fatalf("RTS smooth cov vs Joint conditional mismatch: %g (>= 1e-9)", anc.MaxSmoothCovAbsDiff)
	}

	fastLL, err := FastLogLikelihood(spec, taus, WeeklyDt, yields)
	if err != nil {
		t.Fatalf("FastLogLikelihood failed: %v", err)
	}
	if math.Abs(fastLL-anc.JointLogLik) >= 1e-9 {
		t.Fatalf("FastLogLikelihood vs Joint mismatch: %g", math.Abs(fastLL-anc.JointLogLik))
	}
}

func TestFastLogLikelihoodMatchesFullFilterWithRiccatiSteadyState(t *testing.T) {
	cases := []struct {
		name  string
		spec  ModelSpec
		dt    float64
		steps int
	}{
		{name: "weekly_1000", spec: DefaultSpec(), dt: WeeklyDt, steps: 1000},
		{name: "daily_slow_kp11_5040", spec: VecToSpec([11]float64{0.6, 5e-4, 0.5, 0.4, 0.045, -0.02, 0.01, 0.006, 0.01, 0.02, 0.0005}, true), dt: DailyBusinessDt, steps: 5040},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xB002, uint64(tc.steps)))
			_, yields := afns.SimulateStationaryPath(tc.spec.AFNS, StandardTenors, tc.steps, tc.dt, tc.spec.SigmaEps, rg)

			fastLL, err := FastLogLikelihood(tc.spec, StandardTenors, tc.dt, yields)
			if err != nil {
				t.Fatalf("FastLogLikelihood failed: %v", err)
			}
			fullOut, err := FilterAndSmooth(tc.spec, StandardTenors, tc.dt, yields)
			if err != nil {
				t.Fatalf("FilterAndSmooth failed: %v", err)
			}
			diff := math.Abs(fastLL - fullOut.LogLikelihood)
			if diff >= 1e-9 {
				t.Fatalf("Riccati steady-state FastLL diverged from full FilterAndSmooth: diff=%g (fast=%v, full=%v)", diff, fastLL, fullOut.LogLikelihood)
			}
		})
	}
}

func TestScaledGradientVsRichardsonOnDailyAndWeekly(t *testing.T) {
	trueSpec := DefaultSpec()
	for _, tc := range []struct {
		name  string
		dt    float64
		steps int
	}{
		{name: "weekly_400", dt: WeeklyDt, steps: 400},
		{name: "daily_1000", dt: DailyBusinessDt, steps: 1000},
	} {
		t.Run(tc.name, func(t *testing.T) {
			rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0xC003, uint64(tc.steps)))
			_, yields := afns.SimulateStationaryPath(trueSpec.AFNS, StandardTenors, tc.steps, tc.dt, trueSpec.SigmaEps, rg)
			initSpec := PerturbInitialGuess(trueSpec, rg)

			negLL := func(u []float64) float64 {
				ll, err := FastLogLikelihood(VecToSpec(UnpackUnconstrained(u), true), StandardTenors, tc.dt, yields)
				if err != nil {
					return 1e30
				}
				return -ll
			}
			_, _, errInit := optim.VerifyGradVsRichardson(negLL, PackUnconstrained(SpecToVec(initSpec)), nil)
			_, _, errTrue := optim.VerifyGradVsRichardson(negLL, PackUnconstrained(SpecToVec(trueSpec)), nil)
			if errInit >= 1e-5 {
				t.Fatalf("scoreFD vs Richardson at initSpec too large: %g", errInit)
			}
			if errTrue >= 1e-5 {
				t.Fatalf("scoreFD vs Richardson at trueSpec too large: %g", errTrue)
			}
		})
	}
}
