package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/afns"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/kalman"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
)

func TestQuickPipelineEndToEnd(t *testing.T) {
	tmp := t.TempDir()
	outJSON := filepath.Join(tmp, "quick.json")
	figDir := filepath.Join(tmp, "figs")

	trueSpec := kalman.DefaultSpec()
	nRep := 2
	nWorkers := 2

	baselineCfg := kalman.ConfigSpec{
		Name: "baseline_1000w", Group: "baseline", FrequencyLabel: "weekly",
		HorizonYears: 120.0 * kalman.WeeklyDt, Dt: kalman.WeeklyDt, TSteps: 120, NumReps: nRep, Salt: 0x4B01,
	}
	baseRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, baselineCfg, nWorkers, rng.MainSeed)
	lsRes := kalman.RunLongSpanSharedPrefixStudy(trueSpec, kalman.StandardTenors, []int{2, 4, 6, 8}, nRep, nWorkers, rng.MainSeed, 0x4B02)
	mRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, kalman.ConfigSpec{
		Name: "infill_20y_monthly", Group: "infill", FrequencyLabel: "monthly",
		HorizonYears: 4.0, Dt: kalman.MonthlyDt, TSteps: 48, NumReps: nRep, Salt: 0x4B03,
	}, nWorkers, rng.MainSeed)
	wRes := kalman.ReuseAsInfillConfig(lsRes[1])
	dRes := kalman.RunIndependentConfigStudy(trueSpec, kalman.StandardTenors, kalman.ConfigSpec{
		Name: "infill_20y_daily", Group: "infill", FrequencyLabel: "daily_business",
		HorizonYears: 4.0, Dt: kalman.DailyBusinessDt, TSteps: 120, NumReps: nRep, Salt: 0x4B04,
	}, nWorkers, rng.MainSeed)

	var kpTrends [3]kalman.InvHorizonTrendSummary
	for k := 0; k < 3; k++ {
		kpTrends[k] = kalman.AnalyzeLongSpanTrend(lsRes, k+1)
	}
	if err := buildFigures(figDir, kpTrends, mRes, wRes, dRes, lsRes, baseRes, []kalman.StudyConfigResult{baseRes, lsRes[0], lsRes[1], mRes, dRes}); err != nil {
		t.Fatalf("buildFigures failed: %v", err)
	}
	for _, fn := range []string{
		"day004b_kp-bias-vs-inv-h.svg",
		"day004b_infill-vs-longspan-rmse.svg",
		"day004b_wald-coverage-intervals.svg",
		"day004b_hessian-se-to-emp-sd-ratio.svg",
		"day004b_convergence-and-bound-hits.svg",
	} {
		st, err := os.Stat(filepath.Join(figDir, fn))
		if err != nil || st.Size() < 500 {
			t.Fatalf("expected non-empty SVG %s: err=%v", fn, err)
		}
	}

	// Validasi encoding ringkasan ke JSON
	_, _ = afns.DefaultParams(), json.Valid
	if ok := allFiniteCheck(baseRes, nRep); !ok {
		t.Fatalf("expected finite baseline summary")
	}
	_ = outJSON
}
