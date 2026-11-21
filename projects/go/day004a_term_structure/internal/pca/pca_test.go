package pca

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/kalman"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
)

func TestPCAEigenvaluesCosineSignPatternsAndWorkers(t *testing.T) {
	spec := kalman.DefaultSpec()
	res1, err := RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 2000, 20, 4, 1, rng.MainSeed)
	if err != nil {
		t.Fatalf("RunPCAStudy failed: %v", err)
	}
	res4, err := RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 2000, 20, 4, 4, rng.MainSeed)
	if err != nil {
		t.Fatalf("RunPCAStudy 4 workers failed: %v", err)
	}
	for k := range res1.EmpMeanEigenvalues {
		if math.Float64bits(res1.EmpMeanEigenvalues[k]) != math.Float64bits(res4.EmpMeanEigenvalues[k]) {
			t.Fatalf("PCA worker invariance violated at k=%d", k)
		}
	}

	if res1.TheoSignChanges != [3]int{0, 1, 2} {
		t.Fatalf("unexpected theoretical sign pattern: %v (want [0 1 2])", res1.TheoSignChanges)
	}
	if res1.EmpSignChanges != [3]int{0, 1, 2} {
		t.Fatalf("unexpected empirical sign pattern: %v (want [0 1 2])", res1.EmpSignChanges)
	}
	for k := 0; k < 3; k++ {
		if res1.EmpMeanTop3Cos[k] <= 0.99 {
			t.Fatalf("factor %d cosine similarity %.6f <= 0.99", k+1, res1.EmpMeanTop3Cos[k])
		}
		if res1.ZScoreSampleMean[k] > 3.0 {
			t.Fatalf("factor %d eigenvalue mean z-score %.2f > 3.0", k+1, res1.ZScoreSampleMean[k])
		}
	}
	if res1.NoiseBulkZScore > 3.0 {
		t.Fatalf("noise bulk mean z-score %.2f > 3.0 (emp=%g vs theo=%g)", res1.NoiseBulkZScore, res1.NoiseBulkEmpMean, res1.NoiseBulkTheoMean)
	}
	t.Logf("PCA top3 var=%.2f%% | signChanges=%v | cos=[%.5f, %.5f, %.5f] | bulkZ=%.2f",
		res1.TheoTop3CumVar*100.0, res1.TheoSignChanges,
		res1.EmpMeanTop3Cos[0], res1.EmpMeanTop3Cos[1], res1.EmpMeanTop3Cos[2], res1.NoiseBulkZScore)
}
