package dp

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/battery"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/prices"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

func bruteForceBestProfit(pricePath []float64, bat battery.Spec) (bestNet, bestGross float64, bestSeq []battery.ActionID) {
	H := len(pricePath)
	bestGross = math.Inf(-1)
	seq := make([]battery.ActionID, H)
	bestSeq = make([]battery.ActionID, H)

	var dfs func(h, socIdx int, cashSum float64)
	dfs = func(h, socIdx int, cashSum float64) {
		if h == H {
			tot := cashSum + bat.TerminalValue(socIdx)
			if tot > bestGross {
				bestGross = tot
				copy(bestSeq, seq)
			}
			return
		}
		for aIdx := 0; aIdx < battery.NumActions; aIdx++ {
			act := battery.Actions[aIdx]
			sNext, ok := bat.NextSoCIdx(socIdx, act)
			if !ok {
				continue
			}
			seq[h] = act.ID
			dfs(h+1, sNext, cashSum+bat.HourlyCashflow(pricePath[h], act))
		}
	}
	dfs(0, bat.InitialSoCIdx, 0.0)
	bestNet = bestGross - bat.TerminalValue(bat.InitialSoCIdx)
	return bestNet, bestGross, bestSeq
}

func TestBellmanMatchesBruteForceSmallDeterministic(t *testing.T) {
	bat := battery.DefaultSpec()
	priceCases := [][]float64{
		{-40.0, 35.0, 240.0, 10.0},
		{150.0, -80.0, 30.0, 420.0, 5.0},
	}
	for idx, path := range priceCases {
		bfNet, bfGross, bfActs := bruteForceBestProfit(path, bat)
		dpNet, dpGross, dpActs, _ := SolveDeterministicDP(path, bat)
		if math.Abs(dpNet-bfNet) > 1e-12 || math.Abs(dpGross-bfGross) > 1e-12 {
			t.Fatalf("case %d profit mismatch: dpNet=%.15f bfNet=%.15f", idx, dpNet, bfNet)
		}
		for h := range path {
			if dpActs[h] != bfActs[h] {
				t.Fatalf("case %d hour %d action mismatch: dp=%v bf=%v", idx, h, dpActs[h], bfActs[h])
			}
		}
	}
}

func TestStochasticDPMatchesBruteForceOnDeterministicKernel(t *testing.T) {
	bat := battery.DefaultSpec()
	grid := NewCoarseGrid17()
	// Uji kasus kecil 4 jam (h=20..23) menggunakan kernel deterministik probabilitas 1 pada bin harga tetap.
	// Harga representatif bin: bin 3 (-75), bin 7 (30), bin 12 (275), bin 6 (10).
	testBins := [4]int{3, 7, 12, 6}
	path := make([]float64, 4)
	for i, b := range testBins {
		path[i] = grid.RepPrices[b]
	}
	bfNet, bfGross, bfActs := bruteForceBestProfit(path, bat)

	nBins := grid.NumBins
	nStates := prices.HoursPerDay * regime.NumRegimes * nBins
	stride := regime.NumRegimes * nBins
	probs := make([]float64, nStates*stride)
	for h := 0; h < prices.HoursPerDay; h++ {
		nextB := testBins[3]
		if h >= 20 && h < 23 {
			nextB = testBins[h-20+1]
		}
		for r := 0; r < regime.NumRegimes; r++ {
			for b := 0; b < nBins; b++ {
				stIdx := (h*regime.NumRegimes+r)*nBins + b
				probs[stIdx*stride+regime.Normal*nBins+nextB] = 1.0
			}
		}
	}
	kern := TransitionKernel{
		Grid:         grid,
		NumStates:    nStates,
		TargetStride: stride,
		Probs:        probs,
	}
	pol := SolveStochasticDP(kern, bat)
	if pol.MaxBellmanResidual > 1e-11 {
		t.Fatalf("Bellman residual too high: %.3e", pol.MaxBellmanResidual)
	}
	dpGross := pol.ValueAt(20, regime.Normal, testBins[0], bat.InitialSoCIdx)
	if math.Abs(dpGross-bfGross) > 1e-11 {
		t.Fatalf("stochastic DP on 4h deterministic kernel mismatch: dp=%.12f bf=%.12f (net=%.12f)", dpGross, bfGross, bfNet)
	}
	soc := bat.InitialSoCIdx
	for step := 0; step < 4; step++ {
		act := pol.ActionAt(20+step, regime.Normal, testBins[step], soc)
		if act != bfActs[step] {
			t.Fatalf("step %d action mismatch: dp=%v bf=%v", step, act, bfActs[step])
		}
		soc, _ = bat.NextSoCIdx(soc, battery.Actions[act])
	}
}

func TestGridRefinementAndTailBins(t *testing.T) {
	g17 := NewCoarseGrid17()
	g33 := NewFineGrid33()
	if g17.NumBins != 17 || g33.NumBins != 33 {
		t.Fatalf("unexpected bin counts: %d, %d", g17.NumBins, g33.NumBins)
	}
	if g17.BinIndex(-950.0) != 0 || g17.BinIndex(2500.0) != 16 {
		t.Fatalf("coarse tail bin mapping failed")
	}
	if g33.BinIndex(-950.0) != 0 || g33.BinIndex(2500.0) != 32 {
		t.Fatalf("fine tail bin mapping failed")
	}
	// Semua edge grid 17 bin harus menjadi subset grid 33 bin.
	edgeSet := make(map[float64]bool, len(g33.Edges))
	for _, e := range g33.Edges {
		edgeSet[e] = true
	}
	for _, e := range g17.Edges {
		if !edgeSet[e] {
			t.Fatalf("coarse edge %.1f missing in fine grid", e)
		}
	}
}

func TestPerfectForesightDominanceAndWorkerDeterminism(t *testing.T) {
	spec := prices.DefaultSpec()
	pi, err := regime.StationaryDistribution(spec.TransitionMatrix)
	if err != nil {
		t.Fatalf("stationary error: %v", err)
	}
	bat := battery.DefaultSpec()
	g17 := NewCoarseGrid17()

	k1 := BuildKernel(spec, g17, 250, 1, rng.KernelSeed, 0x7101)
	k4 := BuildKernel(spec, g17, 250, 4, rng.KernelSeed, 0x7101)
	if len(k1.Probs) != len(k4.Probs) {
		t.Fatalf("kernel length mismatch")
	}
	for i := range k1.Probs {
		if math.Float64bits(k1.Probs[i]) != math.Float64bits(k4.Probs[i]) {
			t.Fatalf("kernel prob %d diverged across 1 vs 4 workers", i)
		}
	}

	pol1 := SolveStochasticDP(k1, bat)
	pol4 := SolveStochasticDP(k4, bat)
	if pol1.MaxBellmanResidual > 1e-11 {
		t.Fatalf("Bellman residual too large: %.3e", pol1.MaxBellmanResidual)
	}

	days := prices.SimulateEvalDays(60, 4, spec, pi, rng.EvalSeed, 0x7102)
	res1 := EvalStochasticDP(StrategyDP17, days, pol1, bat)
	res4 := EvalStochasticDP(StrategyDP17, days, pol4, bat)
	for i := range res1.DailyProfits {
		if math.Float64bits(res1.DailyProfits[i]) != math.Float64bits(res4.DailyProfits[i]) {
			t.Fatalf("day %d profit diverged across 1 vs 4 workers", i)
		}
	}

	thresh := EvalThresholdPolicy(days, bat)
	pf := EvalPerfectForesightDP(days, bat)
	nobat := EvalNoBattery(days)

	for _, sr := range []*StrategyResult{&res1, &thresh, &nobat} {
		ComputeRegretVsPF(sr, pf, 1e-9)
		if sr.PerfectForesightDomViolations != 0 {
			t.Fatalf("strategy %s exceeded perfect foresight %d times (max shortfall %.3e)",
				sr.Name, sr.PerfectForesightDomViolations, sr.MaxDominationShortfall)
		}
		if sr.SoCViolationCount != 0 || sr.PowerViolationCount != 0 {
			t.Fatalf("strategy %s violated battery bounds: soc=%d pow=%d", sr.Name, sr.SoCViolationCount, sr.PowerViolationCount)
		}
	}
}

func TestExpandedTailGridOverflowAndRepresentation(t *testing.T) {
	g33 := NewFineGrid33()
	gExp := NewExpandedTailGrid()
	if gExp.LowerTailAt > -1000.0 || gExp.UpperTailAt < 5000.0 {
		t.Fatalf("expanded grid does not cover [-1000, 5000]: [%.1f, %.1f]", gExp.LowerTailAt, gExp.UpperTailAt)
	}
	edgeSet := make(map[float64]bool, len(gExp.Edges))
	for _, e := range gExp.Edges {
		edgeSet[e] = true
	}
	for _, e := range g33.Edges {
		if !edgeSet[e] {
			t.Fatalf("fine_33 edge %.1f missing from expanded_tail_38", e)
		}
	}
	// Overflow dan batas ekstrem
	cases := []struct {
		p   float64
		bin int
	}{
		{-1e6, 0},
		{-1000.01, 0},
		{-1000.0, 1},
		{-500.0, 3},
		{1999.99, 33},
		{2000.0, 34},
		{4999.99, 36},
		{5000.0, 37},
		{1e6, 37},
	}
	for _, tc := range cases {
		if got := gExp.BinIndex(tc.p); got != tc.bin {
			t.Fatalf("BinIndex(%.4f) = %d, want %d", tc.p, got, tc.bin)
		}
	}
}

func TestConditionalSeedCRNAndWorkerInvarianceExpanded(t *testing.T) {
	spec := prices.DefaultSpec()
	g33 := NewFineGrid33()
	gExp := NewExpandedTailGrid()

	crnMapper := func(h, r, b int) (uint64, int) {
		if b >= 3 && b <= 33 {
			return 0x4E33, b - 2
		}
		return 0x4E38, b
	}
	k33 := BuildKernel(spec, g33, 200, 2, rng.KernelSeed, 0x4E33)
	kExp1 := BuildKernelWithSeedMapper(spec, gExp, 200, 1, rng.KernelSeed, 0x4E38, crnMapper)
	kExp4 := BuildKernelWithSeedMapper(spec, gExp, 200, 4, rng.KernelSeed, 0x4E38, crnMapper)

	for i := range kExp1.Probs {
		if math.Float64bits(kExp1.Probs[i]) != math.Float64bits(kExp4.Probs[i]) {
			t.Fatalf("expanded kernel prob %d diverged across 1 vs 4 workers", i)
		}
	}
	// Verifikasi bahwa jumlah probabilitas setiap state tepat 1.0 (massa tail tidak dibuang)
	for st := 0; st < kExp1.NumStates; st++ {
		sum := 0.0
		row := kExp1.Probs[st*kExp1.TargetStride : (st+1)*kExp1.TargetStride]
		for _, p := range row {
			sum += p
		}
		if math.Abs(sum-1.0) > 1e-12 {
			t.Fatalf("state %d kernel prob sum = %.15f != 1.0", st, sum)
		}
	}
	// Pada bin interior bersama (mis. b33=15 <-> bExp=17, harga 45), cek bahwa jumlah probabilitas masuk ke interior g33 cocok
	if math.Abs(g33.RepPrices[15]-gExp.RepPrices[17]) > 1e-12 {
		t.Fatalf("rep price mismatch: %.3f vs %.3f", g33.RepPrices[15], gExp.RepPrices[17])
	}
	for r := 0; r < regime.NumRegimes; r++ {
		st33 := k33.StateIndex(12, r, 15)
		stExp := kExp1.StateIndex(12, r, 17)
		for b33 := 1; b33 < g33.NumBins-1; b33++ {
			p33 := k33.Probs[st33*k33.TargetStride+r*g33.NumBins+b33]
			pExp := kExp1.Probs[stExp*kExp1.TargetStride+r*gExp.NumBins+(b33+2)]
			if math.Abs(p33-pExp) > 1e-12 {
				t.Fatalf("CRN interior prob mismatch at r=%d b33=%d: %.12f vs %.12f", r, b33, p33, pExp)
			}
		}
	}
}
