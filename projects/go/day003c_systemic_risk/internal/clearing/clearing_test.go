package clearing_test

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/clearing"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
)

func TestPicardFictitiousBruteForce200Networks(t *testing.T) {
	r := rng.New(20261121)
	maxDiffPF := 0.0
	maxDiffPB := 0.0

	for trial := 0; trial < 210; trial++ {
		n := 4 + (trial % 6) // n in 4..9
		mode := trial % 3
		net := graphs.BuildSmallRandomNetwork(n, r, mode)

		resP := clearing.Picard(net, 1e-13, 5000)
		resF := clearing.FictitiousDefault(net)
		resB := clearing.BruteForce(net)

		for i := 0; i < n; i++ {
			if d := math.Abs(resP.P[i] - resF.P[i]); d > maxDiffPF {
				maxDiffPF = d
			}
			if d := math.Abs(resP.P[i] - resB.P[i]); d > maxDiffPB {
				maxDiffPB = d
			}
			if resP.P[i] > net.PBar[i]+1e-12 || resP.P[i] < -1e-12 {
				t.Fatalf("trial %d: p_i di luar [0, p_bar]", trial)
			}
		}
		if resP.TotalLoss < -1e-12 {
			t.Fatalf("trial %d: total loss negatif", trial)
		}
	}
	if maxDiffPF > 1e-10 {
		t.Fatalf("Picard vs FictitiousDefault maxDiff = %e > 1e-10", maxDiffPF)
	}
	if maxDiffPB > 1e-10 {
		t.Fatalf("Picard vs BruteForce maxDiff = %e > 1e-10", maxDiffPB)
	}
}

func TestExactChainAndTwoBankCycleAnchors(t *testing.T) {
	// Anchor 1: Rantai 0->1->...->n-1->sink (L_{i,i+1} = 1)
	// Tambahkan simpul sink n dengan e_n = 0, L_{n,*} = 0 agar p_bar_0..p_bar_{n-1} = 1.
	n := 8
	buildChain := func(e []float64) *graphs.InterbankNet {
		eFull := make([]float64, n+1)
		copy(eFull, e)
		l := make([][]float64, n+1)
		for i := range l {
			l[i] = make([]float64, n+1)
		}
		for i := 0; i < n; i++ {
			l[i][i+1] = 1.0
		}
		return graphs.NewFromDense(eFull, l)
	}

	// Pola A: e_i = 0.08 konstan -> p_i* = (i+1)*0.08
	eA := make([]float64, n)
	for i := 0; i < n; i++ {
		eA[i] = 0.08
	}
	netA := buildChain(eA)
	resA := clearing.FictitiousDefault(netA)
	for i := 0; i < n; i++ {
		want := float64(i+1) * 0.08
		if d := math.Abs(resA.P[i] - want); d > 1e-14 {
			t.Fatalf("Chain pattern A idx %d: got %.16f want %.16f diff %e", i, resA.P[i], want, d)
		}
	}

	// Pola B: e_0 = 0.25, e_i = 0.20 (i >= 1) -> p_i* = min(1.0, 0.25 + 0.20*i)
	eB := make([]float64, n)
	eB[0] = 0.25
	for i := 1; i < n; i++ {
		eB[i] = 0.20
	}
	netB := buildChain(eB)
	resB := clearing.FictitiousDefault(netB)
	for i := 0; i < n; i++ {
		want := math.Min(1.0, 0.25+0.20*float64(i))
		if d := math.Abs(resB.P[i] - want); d > 1e-14 {
			t.Fatalf("Chain pattern B idx %d: got %.16f want %.16f diff %e", i, resB.P[i], want, d)
		}
	}

	// Anchor 2: Dua bank saling berutang (0 <-> 1) dengan sink 2:
	// p_bar_0 = 1, Pi_{01} = g0 = 0.6, Pi_{02} = 0.4
	// p_bar_1 = 1, Pi_{10} = g1 = 0.5, Pi_{12} = 0.5
	// e_0 = 0.2, e_1 = 0.3 -> keduanya gagal, solusi tertutup:
	// p_0* = (e_0 + g1*e_1)/(1 - g0*g1), p_1* = (e_1 + g0*e_0)/(1 - g0*g1)
	g0, g1 := 0.6, 0.5
	e0, e1 := 0.2, 0.3
	l2 := [][]float64{
		{0, g0, 1 - g0},
		{g1, 0, 1 - g1},
		{0, 0, 0},
	}
	net2 := graphs.NewFromDense([]float64{e0, e1, 0}, l2)
	res2 := clearing.FictitiousDefault(net2)
	want0 := (e0 + g1*e1) / (1.0 - g0*g1)
	want1 := (e1 + g0*e0) / (1.0 - g0*g1)
	if math.Abs(res2.P[0]-want0) > 1e-14 || math.Abs(res2.P[1]-want1) > 1e-14 {
		t.Fatalf("Two-bank cycle mismatch: got (%.16f, %.16f) want (%.16f, %.16f)", res2.P[0], res2.P[1], want0, want1)
	}
}

func TestMonotonicityAndRegularUniqueness(t *testing.T) {
	r := rng.New(20261122)
	net1 := graphs.BuildSmallRandomNetwork(10, r, 0)
	for i := range net1.E {
		net1.E[i] += 0.05 // pastikan setiap simpul punya e_i > 0 (jaringan teratur)
	}
	eHigher := make([]float64, net1.N)
	for i := range eHigher {
		eHigher[i] = net1.E[i] + 0.2
	}
	net2 := net1.CloneWithE(eHigher)

	pTop := clearing.Picard(net1, 1e-13, 5000)
	pBot := clearing.PicardFromZero(net1, 1e-13, 5000)
	pHigh := clearing.Picard(net2, 1e-13, 5000)

	for i := 0; i < net1.N; i++ {
		if math.Abs(pTop.P[i]-pBot.P[i]) > 1e-10 {
			t.Fatalf("Ketidakunikan pada jaringan teratur idx %d: top=%f bot=%f", i, pTop.P[i], pBot.P[i])
		}
		if pHigh.P[i] < pTop.P[i]-1e-12 {
			t.Fatalf("Pelanggaran monotonisitas idx %d: high=%f low=%f", i, pHigh.P[i], pTop.P[i])
		}
	}
}
