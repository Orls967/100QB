package bootstrap

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
)

func TestCalibrationStudyWorkerDeterminism(t *testing.T) {
	spec := kalman.DefaultSpec()
	outer1 := RunOuterFits(spec, kalman.StandardTenors, kalman.WeeklyDt, 80, 2, 1, rng.MainSeed, 0xD101)
	outer4 := RunOuterFits(spec, kalman.StandardTenors, kalman.WeeklyDt, 80, 2, 4, rng.MainSeed, 0xD101)

	res1 := RunCalibrationStudy(spec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, outer1, 0.0, 3, 1, rng.MainSeed, 0xD102)
	res4 := RunCalibrationStudy(spec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, outer4, 0.0, 3, 4, rng.MainSeed, 0xD102)

	for r := 0; r < 2; r++ {
		if math.Float64bits(res1.OuterReplications[r].OuterMLE.EstLogLik) != math.Float64bits(res4.OuterReplications[r].OuterMLE.EstLogLik) {
			t.Fatalf("rep %d outer EstLogLik diverged across 1 vs 4 workers", r)
		}
		for b := 0; b < 3; b++ {
			d1 := res1.OuterReplications[r].InnerDraws[b]
			d4 := res4.OuterReplications[r].InnerDraws[b]
			if math.Float64bits(d1.EstLogLik) != math.Float64bits(d4.EstLogLik) ||
				math.Float64bits(d1.HessianSE[1]) != math.Float64bits(d4.HessianSE[1]) ||
				math.Float64bits(d1.HessianSE[4]) != math.Float64bits(d4.HessianSE[4]) {
				t.Fatalf("rep %d boot %d diverged across 1 vs 4 workers", r, b)
			}
		}
	}
}
