package kalman

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
)

func TestStudyWorkerDeterminismAndTrendLinearity(t *testing.T) {
	spec := DefaultSpec()
	cfg := ConfigSpec{
		Name:           "fast_determinism_fixture",
		Group:          "test",
		FrequencyLabel: "weekly",
		HorizonYears:   100.0 * WeeklyDt,
		Dt:             WeeklyDt,
		TSteps:         100,
		NumReps:        4,
		Salt:           0xDEAD,
	}

	res1 := RunIndependentConfigStudy(spec, StandardTenors, cfg, 1, rng.MainSeed)
	res4 := RunIndependentConfigStudy(spec, StandardTenors, cfg, 4, rng.MainSeed)

	for r := 0; r < cfg.NumReps; r++ {
		if math.Float64bits(res1.Replications[r].EstLogLik) != math.Float64bits(res4.Replications[r].EstLogLik) {
			t.Fatalf("rep %d EstLogLik mismatch 1 vs 4 workers: %v vs %v", r, res1.Replications[r].EstLogLik, res4.Replications[r].EstLogLik)
		}
		for i := 0; i < 11; i++ {
			if math.Float64bits(res1.Replications[r].EstParams[i]) != math.Float64bits(res4.Replications[r].EstParams[i]) {
				t.Fatalf("rep %d param %d mismatch 1 vs 4 workers", r, i)
			}
			if math.Float64bits(res1.Replications[r].HessianSE[i]) != math.Float64bits(res4.Replications[r].HessianSE[i]) {
				t.Fatalf("rep %d SE %d mismatch 1 vs 4 workers", r, i)
			}
		}
	}
	for i := 0; i < 11; i++ {
		if math.Float64bits(res1.Params[i].Bias) != math.Float64bits(res4.Params[i].Bias) {
			t.Fatalf("summary param %d bias mismatch 1 vs 4 workers", i)
		}
	}

	// Verifikasi studi long-span prefiks-bersama + konsistensi aljabar MCSE prefiks
	ls1 := RunLongSpanSharedPrefixStudy(spec, StandardTenors, []int{2, 4, 8}, 4, 1, rng.MainSeed, 0xCAFE)
	ls4 := RunLongSpanSharedPrefixStudy(spec, StandardTenors, []int{2, 4, 8}, 4, 4, rng.MainSeed, 0xCAFE)
	for h := range ls1 {
		for r := 0; r < 4; r++ {
			if math.Float64bits(ls1[h].Replications[r].EstLogLik) != math.Float64bits(ls4[h].Replications[r].EstLogLik) {
				t.Fatalf("longspan horizon %d rep %d diverged 1 vs 4 workers", h, r)
			}
		}
	}
	// Verifikasi bahwa dua replikasi gagal/batas pun tidak pernah dibuang dari ringkasan
	if len(ls1[0].Replications) != 4 || ls1[0].Params[1].TotalReps != 4 {
		t.Fatalf("unexpected total reps retained: %d", ls1[0].Params[1].TotalReps)
	}
	tr := AnalyzeLongSpanTrend(ls1, 1)
	if len(tr.Horizons) != 3 || tr.OLSSlopePrefixMCSE <= 0 || tr.WLSSlopePrefixMCSE <= 0 {
		t.Fatalf("invalid trend summary: %+v", tr)
	}
	infill := AnalyzeInfillScaling(ls1[0], ReuseAsInfillConfig(ls1[1]), ls1[2])
	if infill[1].SqrtNObsRatio <= 0 || infill[1].Name != "kp11" {
		t.Fatalf("unexpected infill summary: %+v", infill[1])
	}
}
