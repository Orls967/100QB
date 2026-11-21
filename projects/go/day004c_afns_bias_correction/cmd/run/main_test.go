package main

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/jackknife"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
)

func TestQuickRunnerPipelineAndSVGs(t *testing.T) {
	figDir := filepath.Join(t.TempDir(), "figs")
	spec := kalman.DefaultSpec()
	paths := jackknife.GenerateSharedPrefixOuterPaths(spec, kalman.StandardTenors, 80, kalman.WeeklyDt, 2, 2, rng.MainSeed, 0x7701)
	jk20 := jackknife.RunJackknifeHorizonStudy("quick_jk20", spec, kalman.StandardTenors, paths, 40.0*kalman.WeeklyDt, kalman.WeeklyDt, 40, 2)
	jk80 := jackknife.RunJackknifeHorizonStudy("quick_jk80", spec, kalman.StandardTenors, paths, 80.0*kalman.WeeklyDt, kalman.WeeklyDt, 80, 2)
	outer := []kalman.MLEResult{jk20.Replications[0].Full, jk20.Replications[1].Full}
	boot := bootstrap.RunParametricBootstrapStudy(spec, kalman.StandardTenors, 40.0*kalman.WeeklyDt, kalman.WeeklyDt, 40, outer, 4, 2, rng.MainSeed, 0x7702)

	if !checkAllFiniteJK(jk20, 2) || !checkAllFiniteJK(jk80, 2) || !checkAllFiniteBoot(boot, 2, 4) {
		t.Fatalf("expected finite complete summaries")
	}
	if !checkSummaryAlgebra(jk20, jk80, boot) {
		t.Fatalf("summary algebra check failed")
	}
	if err := buildFigures(figDir, jk20, jk80, boot, 0); err != nil {
		t.Fatalf("buildFigures error: %v", err)
	}
	for _, fn := range []string{
		"day004c_jackknife-vs-mle-bias-rmse.svg",
		"day004c_bootstrap-vs-jackknife-bias-20y.svg",
		"day004c_wald-vs-bootstrap-coverage.svg",
		"day004c_bootstrap-distributions-rep0.svg",
	} {
		st, err := os.Stat(filepath.Join(figDir, fn))
		if err != nil || st.Size() < 500 {
			t.Fatalf("expected non-empty SVG %s: err=%v", fn, err)
		}
	}
}
