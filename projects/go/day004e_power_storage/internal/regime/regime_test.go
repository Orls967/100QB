package regime

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

func TestTransitionMatrixAndStationary(t *testing.T) {
	P := DefaultTransitionMatrix()
	if err := ValidateTransitionMatrix(P); err != nil {
		t.Fatalf("ValidateTransitionMatrix failed: %v", err)
	}
	pi, err := StationaryDistribution(P)
	if err != nil {
		t.Fatalf("StationaryDistribution failed: %v", err)
	}
	sum := pi[0] + pi[1] + pi[2]
	if math.Abs(sum-1.0) > 1e-14 {
		t.Fatalf("stationary sum != 1: %.17g", sum)
	}
	res := StationaryResidual(P, pi)
	if res > 1e-14 {
		t.Fatalf("stationary residual too large: %.3e", res)
	}

	// Cek silang independen dengan pangkat matriks P^500 (power iteration).
	cur := P
	for step := 0; step < 9; step++ {
		var nxt [NumRegimes][NumRegimes]float64
		for i := 0; i < NumRegimes; i++ {
			for j := 0; j < NumRegimes; j++ {
				for k := 0; k < NumRegimes; k++ {
					nxt[i][j] += cur[i][k] * cur[k][j]
				}
			}
		}
		cur = nxt
	}
	for j := 0; j < NumRegimes; j++ {
		if math.Abs(cur[0][j]-pi[j]) > 1e-12 {
			t.Fatalf("power iteration mismatch for regime %d: %.15f vs %.15f", j, cur[0][j], pi[j])
		}
	}

	// Cek frekuensi empiris rantai Markov panjang.
	const nSteps = 200000
	rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x1101, 0))
	r := SampleStationary(pi, rg)
	var counts [NumRegimes]int
	for s := 0; s < nSteps; s++ {
		r = SampleNext(P, r, rg)
		counts[r]++
	}
	for j := 0; j < NumRegimes; j++ {
		emp := float64(counts[j]) / float64(nSteps)
		if math.Abs(emp-pi[j]) > 0.01 {
			t.Fatalf("empirical stationary frequency %s = %.4f far from pi=%.4f", Names[j], emp, pi[j])
		}
	}
}
