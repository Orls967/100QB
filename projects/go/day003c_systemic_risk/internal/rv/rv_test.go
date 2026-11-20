package rv_test

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/clearing"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rv"
)

func TestRVAlphaBetaOneMatchesENAndChainClosedForm(t *testing.T) {
	r := rng.New(20261121)
	net := graphs.BuildDirectedPoissonNetwork(100, 6.0, r)

	resEN := clearing.Picard(net, 1e-13, 2000)
	resRV1 := rv.Clear(net, 1.0, 1.0, 1e-13, 2000)

	for i := 0; i < net.N; i++ {
		if d := math.Abs(resEN.P[i] - resRV1.P[i]); d > 1e-12 {
			t.Fatalf("alpha=beta=1 mismatch idx %d: %e", i, d)
		}
	}

	// Rantai tertutup dengan biaya kebangkrutan alpha=0.85, beta=0.75
	n := 8
	e := []float64{0.10, 0.08, 0.12, 0.05, 0.09, 0.11, 0.07, 0.06}
	eFull := make([]float64, n+1)
	copy(eFull, e)
	l := make([][]float64, n+1)
	for i := range l {
		l[i] = make([]float64, n+1)
	}
	for i := 0; i < n; i++ {
		l[i][i+1] = 1.0
	}
	chainNet := graphs.NewFromDense(eFull, l)
	resChain := rv.Clear(chainNet, 0.85, 0.75, 1e-14, 1000)
	want := rv.ChainClosedForm(e, 0.85, 0.75)
	for i := 0; i < n; i++ {
		if d := math.Abs(resChain.P[i] - want[i]); d > 1e-12 {
			t.Fatalf("RV chain idx %d: got %.15f want %.15f diff %e", i, resChain.P[i], want[i], d)
		}
	}

	// Monotonisitas dalam (alpha, beta): (0.7, 0.7) <= (0.9, 0.9) <= (1.0, 1.0)
	resLow := rv.Clear(net, 0.7, 0.7, 1e-13, 2000)
	resMid := rv.Clear(net, 0.9, 0.9, 1e-13, 2000)
	for i := 0; i < net.N; i++ {
		if resLow.P[i] > resMid.P[i]+1e-11 || resMid.P[i] > resRV1.P[i]+1e-11 {
			t.Fatalf("Pelanggaran monotonisitas RV idx %d", i)
		}
	}
}

func TestStressDeterminismAndWorkerInvariance(t *testing.T) {
	r := rng.New(20261121)
	net := graphs.BuildDirectedPoissonNetwork(80, 5.0, r)

	s1 := rv.RunMonteCarloStress(net, 0.5, 0.2, 0.9, 0.9, 200, 1, 20261121)
	s4 := rv.RunMonteCarloStress(net, 0.5, 0.2, 0.9, 0.9, 200, 4, 20261121)

	if s1.ENMeanLoss != s4.ENMeanLoss || s1.RVMeanLoss != s4.RVMeanLoss {
		t.Fatalf("Hasil berbeda antara 1 dan 4 worker: %v vs %v", s1.ENMeanLoss, s4.ENMeanLoss)
	}
	if !s1.AllRVGeEN {
		t.Fatalf("Ada skenario dengan kerugian RV < kerugian EN")
	}
}
