package parallel

import (
	"math"
	"testing"
)

func TestRunChunksWorkerInvariance(t *testing.T) {
	fn := func(idx int) float64 {
		x := float64(idx+1) * 0.41
		s := 0.0
		for k := 1; k <= 400; k++ {
			s += math.Cos(x*float64(k)) / float64(k)
		}
		return s
	}

	res1 := RunChunks(43, 1, fn)
	res4 := RunChunks(43, 4, fn)
	res8 := RunChunks(43, 8, fn)
	for i := range res1 {
		b1 := math.Float64bits(res1[i])
		if b1 != math.Float64bits(res4[i]) || b1 != math.Float64bits(res8[i]) {
			t.Fatalf("chunk %d diverged across workers: %v vs %v vs %v", i, res1[i], res4[i], res8[i])
		}
	}
}
