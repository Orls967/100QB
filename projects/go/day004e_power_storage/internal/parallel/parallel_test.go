package parallel

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

func TestRunChunksWorkerInvariance(t *testing.T) {
	fn := func(idx int) uint64 {
		rg := rng.New(rng.DeriveSeed(rng.MainSeed, 0x99, uint64(idx)))
		x := 0.0
		for k := 0; k < 20; k++ {
			x += rg.NormFloat64()
		}
		return math.Float64bits(x)
	}
	r1 := RunChunks(64, 1, fn)
	r8 := RunChunks(64, 8, fn)
	for i := range r1 {
		if r1[i] != r8[i] {
			t.Fatalf("chunk %d diverged across 1 vs 8 workers", i)
		}
	}
}
