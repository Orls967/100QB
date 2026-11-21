package main

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
)

func TestQuickRunnerPipelineAndSVGs(t *testing.T) {
	figDir := filepath.Join(t.TempDir(), "figs")
	spec := kalman.DefaultSpec()
	outer := bootstrap.RunOuterFits(spec, kalman.StandardTenors, kalman.WeeklyDt, 80, 3, 2, rng.MainSeed, 0x7701)
	st := bootstrap.RunCalibrationStudy(spec, kalman.StandardTenors, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, outer, 0.1, 4, 2, rng.MainSeed, 0x7702)

	if !checkComplete(st, 3, 4) {
		t.Fatalf("expected complete calibration study")
	}
	if !checkAlgebra(st) {
		t.Fatalf("summary algebra check failed")
	}
	if err := buildFigures(figDir, st); err != nil {
		t.Fatalf("buildFigures error: %v", err)
	}
	for _, fn := range []string{
		"day004d_coverage-by-method-kp11-theta1.svg",
		"day004d_kp11-endpoints-raw-vs-log.svg",
		"day004d_theta1-interval-widths.svg",
		"day004d_kp11-vs-theta1-width-coverage.svg",
	} {
		st, err := os.Stat(filepath.Join(figDir, fn))
		if err != nil || st.Size() < 500 {
			t.Fatalf("expected non-empty SVG %s: err=%v", fn, err)
		}
	}
}
