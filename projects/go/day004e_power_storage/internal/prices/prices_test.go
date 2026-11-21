package prices

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

func TestOneStepConditionalMomentsAndModuloHour(t *testing.T) {
	spec := DefaultSpec()

	// Verifikasi m(r, h) periodik modulo 24.
	for r := 0; r < regime.NumRegimes; r++ {
		if math.Abs(spec.SeasonalMean(r, 24)-spec.SeasonalMean(r, 0)) > 1e-14 {
			t.Fatalf("seasonal mean not periodic mod 24 for regime %d", r)
		}
		if math.Abs(spec.SeasonalMean(r, 19)-(spec.Mu[r]+15.0)) > 1e-14 {
			t.Fatalf("peak hour 19 seasonal mean mismatch for regime %d", r)
		}
	}

	const nDraws = 200000
	testCases := []struct {
		rNow int
		h    int
		pNow float64
	}{
		{regime.Normal, 18, 55.0},
		{regime.Scarcity, 23, 210.0},
		{regime.Oversupply, 4, -85.0},
	}

	for idx, tc := range testCases {
		rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x3100+uint64(idx), 0))
		var sum, sum2 float64
		for i := 0; i < nDraws; i++ {
			_, pNext, _, _ := spec.Step(tc.rNow, tc.h, tc.pNow, rg)
			sum += pNext
			sum2 += pNext * pNext
		}
		empM := sum / nDraws
		empV := (sum2 - nDraws*empM*empM) / (nDraws - 1)
		theoM, theoV := spec.OneStepConditionalMoments(tc.rNow, tc.h, tc.pNow)
		mcseM := math.Sqrt(theoV / nDraws)

		// Hitung SE varians via momen keempat empiris.
		rg2 := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x3100+uint64(idx), 0))
		var m4 float64
		for i := 0; i < nDraws; i++ {
			_, pNext, _, _ := spec.Step(tc.rNow, tc.h, tc.pNow, rg2)
			d := pNext - empM
			d2 := d * d
			m4 += d2 * d2
		}
		m4 /= nDraws
		mcseV := math.Sqrt(math.Max(0, (m4-empV*empV)/nDraws))

		if math.Abs(empM-theoM) > 4.0*mcseM {
			t.Fatalf("case %d mean mismatch: emp=%.4f theo=%.4f mcse=%.4f", idx, empM, theoM, mcseM)
		}
		if math.Abs(empV-theoV) > 4.0*mcseV {
			t.Fatalf("case %d var mismatch: emp=%.4f theo=%.4f mcse=%.4f", idx, empV, theoV, mcseV)
		}
	}
}

func TestSimulateEvalDaysDeterminismAndUnbounded(t *testing.T) {
	spec := DefaultSpec()
	pi, err := regime.StationaryDistribution(spec.TransitionMatrix)
	if err != nil {
		t.Fatalf("stationary dist error: %v", err)
	}
	d1 := SimulateEvalDays(40, 1, spec, pi, rng.EvalSeed, 0x4E01)
	d4 := SimulateEvalDays(40, 4, spec, pi, rng.EvalSeed, 0x4E01)
	for i := range d1 {
		for h := 0; h < HoursPerDay; h++ {
			if math.Float64bits(d1[i].Prices[h]) != math.Float64bits(d4[i].Prices[h]) ||
				d1[i].Regimes[h] != d4[i].Regimes[h] {
				t.Fatalf("day %d hour %d diverged across 1 vs 4 workers", i, h)
			}
		}
	}
}
