package bootstrap

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/interval"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
)

func TestSummarizeMethodRobustAndPhysicalCounts(t *testing.T) {
	results := []interval.IntervalResult{
		{Method: interval.MethodBasic, Lower: -0.4, Upper: 0.6, Width: 1.0, Valid: true, Covered: true, NonPhysical: true, OutOfBounds: true},
		{Method: interval.MethodBasic, Lower: 0.02, Upper: 0.52, Width: 0.5, Valid: true, Covered: true, NonPhysical: false, OutOfBounds: false},
		{Method: interval.MethodBasic, Lower: 0.10, Upper: 2.10, Width: 2.0, Valid: true, Covered: false, NonPhysical: false, OutOfBounds: false},
		{Method: interval.MethodBasic, Lower: 0.01, Upper: 100.01, Width: 100.0, Valid: true, Covered: true, NonPhysical: false, OutOfBounds: false},
		{Method: interval.MethodBasic, Lower: math.NaN(), Upper: math.NaN(), Width: math.NaN(), Valid: false, InvalidReason: interval.ReasonInsufficientValidSE, Covered: false},
	}
	ms := summarizeMethod(interval.MethodBasic, true, results)
	if ms.TotalOuterReps != 5 || ms.ValidCount != 4 || ms.InvalidCount != 1 {
		t.Fatalf("denominator/validity mismatch: %+v", ms)
	}
	// Coverage tetap menggunakan denominator penuh 5 (3/5 = 0.60), bukan 3/4 = 0.75.
	if ms.CoveredCount != 3 || math.Abs(ms.Coverage95-0.60) > 1e-12 {
		t.Fatalf("coverage mismatch: got covered=%d rate=%.6f", ms.CoveredCount, ms.Coverage95)
	}
	// Interval ke-0 mencakup nilai benar tetapi nonfisik, sehingga hanya 2 interval yang valid fisik & tercakup.
	if ms.NonPhysicalCount != 1 || ms.PhysicallyValidCoveredCount != 2 {
		t.Fatalf("physical validity mismatch: nonphys=%d physCov=%d", ms.NonPhysicalCount, ms.PhysicallyValidCoveredCount)
	}
	if ms.InvalidInsuffValidSE != 1 {
		t.Fatalf("expected InvalidInsuffValidSE=1, got %d", ms.InvalidInsuffValidSE)
	}
	// Lebar valid terurut: [0.5, 1.0, 2.0, 100.0] -> mean=25.875, med=1.5, q25=0.875, q75=26.5, max=100.0
	if math.Abs(ms.MeanWidth-25.875) > 1e-12 || math.Abs(ms.MedianWidth-1.5) > 1e-12 ||
		math.Abs(ms.Q25Width-0.875) > 1e-12 || math.Abs(ms.Q75Width-26.5) > 1e-12 ||
		math.Abs(ms.IQRWidth-(26.5-0.875)) > 1e-12 || math.Abs(ms.MaxWidth-100.0) > 1e-12 {
		t.Fatalf("robust width summary mismatch: %+v", ms)
	}
}

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
