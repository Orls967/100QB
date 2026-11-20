package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"time"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/cascade"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/clearing"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/jsonout"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/parallel"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rv"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/stats"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/svg"
)

const BaseSeed uint64 = 20261121

func findRepoRoot() string {
	wd, _ := os.Getwd()
	for cur := wd; cur != "/" && cur != "."; cur = filepath.Dir(cur) {
		if _, err := os.Stat(filepath.Join(cur, "outputs")); err == nil {
			return cur
		}
	}
	return wd
}

func maxSliceDiff(a, b []float64) float64 {
	mx := 0.0
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > mx {
			mx = d
		}
	}
	return mx
}

func logLogSlope(xs, ys []float64) float64 {
	n := float64(len(xs))
	if n < 2 {
		return 0
	}
	var sx, sy float64
	lx := make([]float64, len(xs))
	ly := make([]float64, len(ys))
	for i := range xs {
		lx[i] = math.Log(xs[i])
		ly[i] = math.Log(ys[i])
		sx += lx[i]
		sy += ly[i]
	}
	mx, my := sx/n, sy/n
	var num, den float64
	for i := range xs {
		num += (lx[i] - mx) * (ly[i] - my)
		den += (lx[i] - mx) * (lx[i] - mx)
	}
	if den <= 0 {
		return 0
	}
	return num / den
}

func main() {
	root := findRepoRoot()
	outDir := filepath.Join(root, "outputs", "day003c")
	figDir := filepath.Join(root, "outputs", "figures")
	nWorkers := parallel.DefaultWorkers()
	tStart := time.Now()

	// =========================================================================
	// BAGIAN C1: Kliring Eisenberg-Noe
	// =========================================================================
	tC1 := time.Now()
	rC1 := rng.New(BaseSeed)

	// 1. 210 jaringan acak kecil (siklik, berlapis, tidak-teratur): Picard vs Fictitious Default vs Brute-Force
	maxDiffPFSmall := 0.0
	maxDiffPBSmall := 0.0
	maxDiffFBSmall := 0.0
	allBounded := true
	allLossNonNeg := true
	for trial := 0; trial < 210; trial++ {
		n := 4 + (trial % 6)
		net := graphs.BuildSmallRandomNetwork(n, rC1, trial%3)
		resP := clearing.Picard(net, 1e-13, 5000)
		resF := clearing.FictitiousDefault(net)
		resB := clearing.BruteForce(net)
		if d := maxSliceDiff(resP.P, resF.P); d > maxDiffPFSmall {
			maxDiffPFSmall = d
		}
		if d := maxSliceDiff(resP.P, resB.P); d > maxDiffPBSmall {
			maxDiffPBSmall = d
		}
		if d := maxSliceDiff(resF.P, resB.P); d > maxDiffFBSmall {
			maxDiffFBSmall = d
		}
		for i := 0; i < n; i++ {
			if resP.P[i] < -1e-12 || resP.P[i] > net.PBar[i]+1e-12 {
				allBounded = false
			}
		}
		if resP.TotalLoss < -1e-12 {
			allLossNonNeg = false
		}
	}

	// 2. Anchor eksak: rantai n=10 (2 pola e) dan siklus 2 bank
	nChain := 10
	buildChain := func(e []float64) *graphs.InterbankNet {
		eFull := make([]float64, nChain+1)
		copy(eFull, e)
		l := make([][]float64, nChain+1)
		for i := range l {
			l[i] = make([]float64, nChain+1)
		}
		for i := 0; i < nChain; i++ {
			l[i][i+1] = 1.0
		}
		return graphs.NewFromDense(eFull, l)
	}
	ePatA := make([]float64, nChain)
	wantPatA := make([]float64, nChain)
	for i := 0; i < nChain; i++ {
		ePatA[i] = 0.07
		wantPatA[i] = float64(i+1) * 0.07
	}
	resPatA := clearing.FictitiousDefault(buildChain(ePatA))
	chainErrA := maxSliceDiff(resPatA.P[:nChain], wantPatA)

	ePatB := make([]float64, nChain)
	wantPatB := make([]float64, nChain)
	ePatB[0] = 0.25
	wantPatB[0] = 0.25
	for i := 1; i < nChain; i++ {
		ePatB[i] = 0.18
		wantPatB[i] = math.Min(1.0, 0.25+0.18*float64(i))
	}
	resPatB := clearing.FictitiousDefault(buildChain(ePatB))
	chainErrB := maxSliceDiff(resPatB.P[:nChain], wantPatB)

	g0, g1 := 0.6, 0.5
	e0, e1 := 0.2, 0.3
	netCycle2 := graphs.NewFromDense([]float64{e0, e1, 0}, [][]float64{
		{0, g0, 1 - g0},
		{g1, 0, 1 - g1},
		{0, 0, 0},
	})
	resCycle2 := clearing.FictitiousDefault(netCycle2)
	wantC0 := (e0 + g1*e1) / (1.0 - g0*g1)
	wantC1 := (e1 + g0*e0) / (1.0 - g0*g1)
	cycle2Err := math.Max(math.Abs(resCycle2.P[0]-wantC0), math.Abs(resCycle2.P[1]-wantC1))

	// 3. Monotonisitas & keunikan pada jaringan teratur
	netReg := graphs.BuildDirectedPoissonNetwork(50, 5.0, rng.New(BaseSeed+10))
	eUp := make([]float64, netReg.N)
	for i := range eUp {
		eUp[i] = netReg.E[i] * 1.25
	}
	resRegTop := clearing.Picard(netReg, 1e-13, 5000)
	resRegBot := clearing.PicardFromZero(netReg, 1e-13, 5000)
	resRegUp := clearing.Picard(netReg.CloneWithE(eUp), 1e-13, 5000)
	regUniquenessDiff := maxSliceDiff(resRegTop.P, resRegBot.P)
	monotoneOK := true
	for i := 0; i < netReg.N; i++ {
		if resRegUp.P[i]+1e-12 < resRegTop.P[i] {
			monotoneOK = false
		}
	}

	// 4. Skala n = 2000 bank: Poisson (mean deg 8) & Power-Law (gamma = 2.5)
	netPois2000 := graphs.BuildDirectedPoissonNetwork(2000, 8.0, rng.New(BaseSeed+100))
	tP0 := time.Now()
	resPoisPicard := clearing.Picard(netPois2000, 1e-12, 2000)
	poisPicardSec := time.Since(tP0).Seconds()
	tF0 := time.Now()
	resPoisFict := clearing.FictitiousDefault(netPois2000)
	poisFictSec := time.Since(tF0).Seconds()
	poisAlgoDiff := maxSliceDiff(resPoisPicard.P, resPoisFict.P)

	netPL2000 := graphs.BuildDirectedPowerLawNetwork(2000, 2.5, 8.0, rng.New(BaseSeed+200))
	tP1 := time.Now()
	resPLPicard := clearing.Picard(netPL2000, 1e-12, 2000)
	plPicardSec := time.Since(tP1).Seconds()
	tF1 := time.Now()
	resPLFict := clearing.FictitiousDefault(netPL2000)
	plFictSec := time.Since(tF1).Seconds()
	plAlgoDiff := maxSliceDiff(resPLPicard.P, resPLFict.P)
	timeC1 := time.Since(tC1).Seconds()

	// =========================================================================
	// BAGIAN C2: Kaskade Watts di Graf Acak (n = 100.000, phi = 0.18, rho0 = 0.001)
	// =========================================================================
	tC2 := time.Now()
	const nWatts = 100000
	const phiBase = 0.18
	const rho0Base = 0.001
	const nSeedsBase = 100 // rho0 * nWatts = 100

	// 1. Bukti invariansi urutan penyebaran (FIFO vs LIFO) pada n = 100.000
	rOrd := rng.New(BaseSeed + 300)
	gOrd := graphs.BuildUndirectedER(nWatts, 4.0, rOrd)
	seedsOrd := rOrd.Perm(nWatts)[:nSeedsBase]
	cntFIFO, maskFIFO := cascade.Simulate(gOrd, seedsOrd, phiBase, false)
	cntLIFO, maskLIFO := cascade.Simulate(gOrd, seedsOrd, phiBase, true)
	orderIdentical := (cntFIFO == cntLIFO)
	for i := 0; i < nWatts; i++ {
		if maskFIFO[i] != maskLIFO[i] {
			orderIdentical = false
			break
		}
	}

	// 2 & 3. S_simulasi (R = 200) vs S_teori untuk c in {0.5, 1, 2, 3, 4, 6, 8, 10, 12}
	cList := []float64{0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0}
	const rReps = 200
	cRows := make([]jsonout.Obj, 0, len(cList))
	simSMeans := make([]float64, len(cList))
	theoSList := make([]float64, len(cList))

	for idxC, cVal := range cList {
		qTheo, sTheo := cascade.MeanField(cVal, rho0Base, phiBase)
		theoSList[idxC] = sTheo

		chunks := parallel.RunChunks(rReps, 40, nWorkers, BaseSeed+1000+uint64(idxC*100), func(_ int, start, end int, cSeed uint64) []float64 {
			r := rng.New(cSeed)
			out := make([]float64, 0, end-start)
			for rep := start; rep < end; rep++ {
				g := graphs.BuildUndirectedER(nWatts, cVal, r)
				seeds := make([]int, nSeedsBase)
				used := make(map[int]struct{}, nSeedsBase)
				for k := 0; k < nSeedsBase; {
					cand := r.Intn(nWatts)
					if _, ok := used[cand]; !ok {
						used[cand] = struct{}{}
						seeds[k] = cand
						k++
					}
				}
				tot, _ := cascade.Simulate(g, seeds, phiBase, false)
				out = append(out, float64(tot)/float64(nWatts))
			}
			return out
		})
		sSamples := make([]float64, 0, rReps)
		for _, ch := range chunks {
			sSamples = append(sSamples, ch...)
		}
		sMean := stats.Mean(sSamples)
		sSE := stats.SE(sSamples)
		z := math.Abs(sMean-sTheo) / math.Max(sSE, 1e-12)
		simSMeans[idxC] = sMean
		// Dekat transisi fase (c1 ~ 1.07, c2 ~ 6.95 -> c in {1, 6, 8})
		nearBoundary := (cVal == 1.0 || cVal == 6.0 || cVal == 8.0)

		cRows = append(cRows, jsonout.Obj{
			{Key: "c", Val: cVal},
			{Key: "q_inf_theory", Val: qTheo},
			{Key: "s_theory", Val: sTheo},
			{Key: "s_sim_mean", Val: sMean},
			{Key: "s_sim_se", Val: sSE},
			{Key: "z_score", Val: z},
			{Key: "near_phase_transition", Val: nearBoundary},
		})
	}

	// 4. Batas fase c1, c2 untuk phi in {0.10, 0.18, 0.30} & scan simulasi resolusi dc = 0.25
	c1Base, c2Base, _ := cascade.PhaseBoundaries(phiBase)
	c1Phi01, c2Phi01, _ := cascade.PhaseBoundaries(0.10)
	c1Phi03, c2Phi03, okPhi03 := cascade.PhaseBoundaries(0.30)

	// Scan c dengan resolusi pra-deklarasi dc = 0.25 pada [0.5, 8.0] (R=20, n=25000 untuk memetakan lonjakan S)
	var scanC, scanSimS, scanTheoS []float64
	simJumpC1, simJumpC2 := 0.0, 0.0
	for cScan := 0.5; cScan <= 8.0+1e-9; cScan += cascade.PredeclaredCGridStep {
		_, st := cascade.MeanField(cScan, rho0Base, phiBase)
		chunks := parallel.RunChunks(20, 20, nWorkers, BaseSeed+5000+uint64(cScan*100), func(_ int, start, end int, cSeed uint64) float64 {
			r := rng.New(cSeed)
			sum := 0.0
			for rep := start; rep < end; rep++ {
				g := graphs.BuildUndirectedER(25000, cScan, r)
				seeds := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24}
				tot, _ := cascade.Simulate(g, seeds, phiBase, false)
				sum += float64(tot) / 25000.0
			}
			return sum
		})
		sm := 0.0
		for _, v := range chunks {
			sm += v
		}
		sm /= 20.0
		scanC = append(scanC, cScan)
		scanSimS = append(scanSimS, sm)
		scanTheoS = append(scanTheoS, st)
		if sm > 0.05 && simJumpC1 == 0.0 {
			simJumpC1 = cScan
		}
		if sm > 0.05 {
			simJumpC2 = cScan
		}
	}

	// 5. Distribusi ukuran kaskade satu benih tunggal di dekat c kritis (c = c1Base, n = 100.000, R = 2000)
	const rSingle = 2000
	singleChunks := parallel.RunChunks(rSingle, 100, nWorkers, BaseSeed+7000, func(_ int, start, end int, cSeed uint64) []int {
		r := rng.New(cSeed)
		out := make([]int, 0, end-start)
		for rep := start; rep < end; rep++ {
			g := graphs.BuildUndirectedER(nWatts, c1Base, r)
			seedNode := r.Intn(nWatts)
			tot, _ := cascade.Simulate(g, []int{seedNode}, phiBase, false)
			out = append(out, tot)
		}
		return out
	})
	singleSizes := make([]int, 0, rSingle)
	globalCnt := 0
	for _, ch := range singleChunks {
		for _, sz := range ch {
			singleSizes = append(singleSizes, sz)
			if float64(sz) >= 0.05*float64(nWatts) {
				globalCnt++
			}
		}
	}
	probGlobalSingle := float64(globalCnt) / float64(rSingle)

	// Hitung CCDF P(Size >= s) pada kaskade lokal (2 <= s < 5000) untuk kemiringan ekor log-log empiris
	sort.Ints(singleSizes)
	var tailX, tailCCDF []float64
	for _, thresh := range []int{2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50} {
		cnt := 0
		for _, sz := range singleSizes {
			if sz >= thresh && sz < 5000 {
				cnt++
			}
		}
		if cnt >= 5 {
			tailX = append(tailX, float64(thresh))
			tailCCDF = append(tailCCDF, float64(cnt)/float64(rSingle))
		}
	}
	ccdfSlope := logLogSlope(tailX, tailCCDF)
	timeC2 := time.Since(tC2).Seconds()

	// =========================================================================
	// BAGIAN C3: Biaya Kebangkrutan Rogers-Veraart & Uji Stres Monte Carlo
	// =========================================================================
	tC3 := time.Now()

	// 1. Verifikasi alpha = beta = 1 => EN, rantai tertutup, dan monotonisitas
	resRV1Pois := rv.Clear(netPois2000, 1.0, 1.0, 1e-12, 2000)
	rv1VsENDiff := maxSliceDiff(resPoisPicard.P, resRV1Pois.P)

	eChainRV := []float64{0.10, 0.08, 0.12, 0.05, 0.09, 0.11, 0.07, 0.06, 0.04, 0.08}
	netChainRV := buildChain(eChainRV)
	resChainRV := rv.Clear(netChainRV, 0.9, 0.9, 1e-14, 1000)
	wantChainRV := rv.ChainClosedForm(eChainRV, 0.9, 0.9)
	rvChainErr := maxSliceDiff(resChainRV.P[:nChain], wantChainRV)

	resRV08 := rv.Clear(netPois2000, 0.8, 0.8, 1e-12, 2000)
	resRV09 := rv.Clear(netPois2000, 0.9, 0.9, 1e-12, 2000)
	rvMonotoneOK := true
	for i := 0; i < netPois2000.N; i++ {
		if resRV08.P[i] > resRV09.P[i]+1e-10 || resRV09.P[i] > resRV1Pois.P[i]+1e-10 {
			rvMonotoneOK = false
		}
	}

	// 2. Uji stres Monte Carlo: 20.000 skenario per s in {0.2, 0.4, 0.6, 0.8}
	shockList := []float64{0.2, 0.4, 0.6, 0.8}
	stressRows := make([]jsonout.Obj, 0, len(shockList))
	stressSummaries := make([]rv.StressSummary, 0, len(shockList))
	for idxS, sVal := range shockList {
		sumS := rv.RunMonteCarloStress(netPois2000, sVal, 0.2, 0.9, 0.9, 20000, nWorkers, BaseSeed+9000+uint64(idxS*100))
		stressSummaries = append(stressSummaries, sumS)
		stressRows = append(stressRows, jsonout.Obj{
			{Key: "shock_s", Val: sVal},
			{Key: "en_mean_default_frac", Val: sumS.ENMeanDefaultFrac},
			{Key: "rv_mean_default_frac", Val: sumS.RVMeanDefaultFrac},
			{Key: "en_mean_loss", Val: sumS.ENMeanLoss},
			{Key: "rv_mean_loss", Val: sumS.RVMeanLoss},
			{Key: "en_prob_default_ge_10pct", Val: sumS.ENProbDefaultGe10Pct},
			{Key: "rv_prob_default_ge_10pct", Val: sumS.RVProbDefaultGe10Pct},
			{Key: "en_var_99", Val: sumS.ENVaR99},
			{Key: "en_es_975", Val: sumS.ENES975},
			{Key: "rv_var_99", Val: sumS.RVVaR99},
			{Key: "rv_es_975", Val: sumS.RVES975},
			{Key: "amplification_rv_over_en", Val: sumS.AmplificationRatio},
			{Key: "all_rv_loss_ge_en_loss", Val: sumS.AllRVGeEN},
		})
	}

	// 3. Perbandingan jaringan Poisson vs Pangkat-Hukum (hub): guncang top-20 hub vs 20 bank acak (s = 0.8)
	evalTargeted := func(net *graphs.InterbankNet, targets []int, sVal float64) (enDef, rvDef, enL, rvL float64) {
		eSh := make([]float64, net.N)
		copy(eSh, net.E)
		for _, idx := range targets {
			eSh[idx] *= (1.0 - sVal)
		}
		ns := net.CloneWithE(eSh)
		rE := clearing.Picard(ns, 1e-10, 500)
		rR := rv.Clear(ns, 0.9, 0.9, 1e-10, 500)
		return float64(rE.DefaultCount) / float64(net.N), float64(rR.DefaultCount) / float64(net.N), rE.TotalLoss, rR.TotalLoss
	}

	evalRandom20 := func(net *graphs.InterbankNet, sVal float64, seed uint64) (enDef, rvDef, enL, rvL float64) {
		r := rng.New(seed)
		const nTrials = 200
		for tr := 0; tr < nTrials; tr++ {
			perm := r.Perm(net.N)[:20]
			ed, rd, el, rl := evalTargeted(net, perm, sVal)
			enDef += ed
			rvDef += rd
			enL += el
			rvL += rl
		}
		return enDef / nTrials, rvDef / nTrials, enL / nTrials, rvL / nTrials
	}

	poisHubs := netPois2000.TopHubs(20)
	plHubs := netPL2000.TopHubs(20)
	poisHubENDef, poisHubRVDef, poisHubENLoss, poisHubRVLoss := evalTargeted(netPois2000, poisHubs, 0.8)
	poisRndENDef, poisRndRVDef, poisRndENLoss, poisRndRVLoss := evalRandom20(netPois2000, 0.8, BaseSeed+12000)
	plHubENDef, plHubRVDef, plHubENLoss, plHubRVLoss := evalTargeted(netPL2000, plHubs, 0.8)
	plRndENDef, plRndRVDef, plRndENLoss, plRndRVLoss := evalRandom20(netPL2000, 0.8, BaseSeed+13000)

	timeC3 := time.Since(tC3).Seconds()
	timeTotal := time.Since(tStart).Seconds()

	// =========================================================================
	// GRAFIK SVG (5 gambar)
	// =========================================================================
	// 1. S vs c teori vs simulasi + jendela kaskade
	_ = svg.WriteLineChart(
		filepath.Join(figDir, "day003c_watts-cascade-window.svg"),
		fmt.Sprintf("Kaskade Watts: Teori Pohon vs Simulasi (n=100k, c1=%.2f, c2=%.2f)", c1Base, c2Base),
		"Derajat Rata-Rata c",
		"Fraksi Gagal Akhir S",
		[]svg.Series{
			{Label: "Teori Pohon (dc=0.25)", X: scanC, Y: scanTheoS, Color: svg.Palette[0], Dashed: false},
			{Label: "Simulasi G(n, c/n) (R=200)", X: cList, Y: simSMeans, Color: svg.Palette[1], Dashed: true},
		},
		false, false,
	)

	// 2. Distribusi ukuran kaskade satu benih di dekat c1 (CCDF log-log)
	_ = svg.WriteLineChart(
		filepath.Join(figDir, "day003c_cascade-size-distribution.svg"),
		fmt.Sprintf("CCDF Ukuran Kaskade Lokal pada c1=%.3f (slope=%.3f)", c1Base, ccdfSlope),
		"Ukuran Kaskade s",
		"P(Size >= s)",
		[]svg.Series{
			{Label: "CCDF Empiris (R=2000)", X: tailX, Y: tailCCDF, Color: svg.Palette[2], Dashed: false},
		},
		true, true,
	)

	// 3. Konvergensi iterasi Picard (n = 2000)
	itPois := make([]float64, len(resPoisPicard.HistoryDiff))
	for i := range itPois {
		itPois[i] = float64(i + 1)
	}
	itPL := make([]float64, len(resPLPicard.HistoryDiff))
	for i := range itPL {
		itPL[i] = float64(i + 1)
	}
	_ = svg.WriteLineChart(
		filepath.Join(figDir, "day003c_picard-convergence.svg"),
		"Konvergensi Iterasi Picard Kliring Eisenberg-Noe (n=2000)",
		"Iterasi k",
		"max |p^(k) - p^(k-1)|",
		[]svg.Series{
			{Label: "Jaringan Poisson (deg=8)", X: itPois, Y: resPoisPicard.HistoryDiff, Color: svg.Palette[0], Dashed: false},
			{Label: "Jaringan Power-Law (gamma=2.5)", X: itPL, Y: resPLPicard.HistoryDiff, Color: svg.Palette[1], Dashed: true},
		},
		false, true,
	)

	// 4. Distribusi kerugian EN vs RV pada uji stres
	var enMeanPlot, rvMeanPlot, enVaRPlot, rvVaRPlot []float64
	for _, sm := range stressSummaries {
		enMeanPlot = append(enMeanPlot, sm.ENMeanLoss)
		rvMeanPlot = append(rvMeanPlot, sm.RVMeanLoss)
		enVaRPlot = append(enVaRPlot, sm.ENVaR99)
		rvVaRPlot = append(rvVaRPlot, sm.RVVaR99)
	}
	_ = svg.WriteLineChart(
		filepath.Join(figDir, "day003c_stress-loss-en-vs-rv.svg"),
		"Kerugian Nilai Sistemik Uji Stres: Eisenberg-Noe vs Rogers-Veraart",
		"Magnitudo Guncangan s",
		"Total Kehilangan Nilai",
		[]svg.Series{
			{Label: "EN Mean Loss", X: shockList, Y: enMeanPlot, Color: svg.Palette[0], Dashed: false},
			{Label: "RV Mean Loss (alpha=beta=0.9)", X: shockList, Y: rvMeanPlot, Color: svg.Palette[1], Dashed: false},
			{Label: "EN VaR 99%", X: shockList, Y: enVaRPlot, Color: svg.Palette[0], Dashed: true},
			{Label: "RV VaR 99%", X: shockList, Y: rvVaRPlot, Color: svg.Palette[1], Dashed: true},
		},
		false, false,
	)

	// 5. Poisson vs Hub (guncangan 20 acak vs top-20 hub)
	_ = svg.WriteBarChart(
		filepath.Join(figDir, "day003c_poisson-vs-hub-shock.svg"),
		"Kehilangan Nilai RV (s=0.8): 20 Bank Acak vs Top-20 Hub",
		"Total Kehilangan Nilai RV",
		[]string{"Poisson Acak-20", "Poisson Hub-20", "PowerLaw Acak-20", "PowerLaw Hub-20"},
		[]float64{poisRndRVLoss, poisHubRVLoss, plRndRVLoss, plHubRVLoss},
	)

	// =========================================================================
	// RESULTS.JSON
	// =========================================================================
	rootObj := jsonout.Obj{
		{Key: "n_workers", Val: nWorkers},
		{Key: "seed", Val: int64(BaseSeed)},
		{Key: "c1_eisenberg_noe", Val: jsonout.Obj{
			{Key: "small_networks_count", Val: 210},
			{Key: "picard_vs_fictitious_max_diff", Val: maxDiffPFSmall},
			{Key: "picard_vs_bruteforce_max_diff", Val: maxDiffPBSmall},
			{Key: "fictitious_vs_bruteforce_max_diff", Val: maxDiffFBSmall},
			{Key: "chain_anchor_pattern_a_err", Val: chainErrA},
			{Key: "chain_anchor_pattern_b_err", Val: chainErrB},
			{Key: "two_bank_cycle_anchor_err", Val: cycle2Err},
			{Key: "monotonicity_verified", Val: monotoneOK},
			{Key: "bounded_p_le_pbar_verified", Val: allBounded},
			{Key: "non_negative_loss_verified", Val: allLossNonNeg},
			{Key: "regular_uniqueness_top_vs_zero_diff", Val: regUniquenessDiff},
			{Key: "scale_n2000_poisson", Val: jsonout.Obj{
				{Key: "picard_iterations", Val: resPoisPicard.Iterations},
				{Key: "picard_sec", Val: poisPicardSec},
				{Key: "fictitious_rounds", Val: resPoisFict.Iterations},
				{Key: "fictitious_sec", Val: poisFictSec},
				{Key: "inter_algo_max_diff", Val: poisAlgoDiff},
				{Key: "default_count", Val: resPoisPicard.DefaultCount},
				{Key: "total_loss", Val: resPoisPicard.TotalLoss},
			}},
			{Key: "scale_n2000_powerlaw", Val: jsonout.Obj{
				{Key: "picard_iterations", Val: resPLPicard.Iterations},
				{Key: "picard_sec", Val: plPicardSec},
				{Key: "fictitious_rounds", Val: resPLFict.Iterations},
				{Key: "fictitious_sec", Val: plFictSec},
				{Key: "inter_algo_max_diff", Val: plAlgoDiff},
				{Key: "default_count", Val: resPLPicard.DefaultCount},
				{Key: "total_loss", Val: resPLPicard.TotalLoss},
			}},
		}},
		{Key: "c2_watts_cascade", Val: jsonout.Obj{
			{Key: "n_nodes", Val: nWatts},
			{Key: "phi", Val: phiBase},
			{Key: "rho0", Val: rho0Base},
			{Key: "order_invariance_fifo_vs_lifo", Val: orderIdentical},
			{Key: "fifo_failed_count", Val: cntFIFO},
			{Key: "lifo_failed_count", Val: cntLIFO},
			{Key: "formula_correction_note", Val: "Rumus prompt memiliki indikator (m+1)/k >= phi pada persamaan q yang sirkular dan gagal di m=0; dikoreksi sesuai Watts (2002) Eq. 3-4 menjadi m/k >= phi."},
			{Key: "c_comparison_r200", Val: cRows},
			{Key: "phase_boundaries", Val: jsonout.Obj{
				{Key: "phi_018_c1", Val: c1Base},
				{Key: "phi_018_c2", Val: c2Base},
				{Key: "predeclared_scan_dc", Val: cascade.PredeclaredCGridStep},
				{Key: "sim_jump_c1", Val: simJumpC1},
				{Key: "sim_jump_c2", Val: simJumpC2},
				{Key: "phi_010_c1", Val: c1Phi01},
				{Key: "phi_010_c2", Val: c2Phi01},
				{Key: "phi_030_has_window", Val: okPhi03},
				{Key: "phi_030_c1", Val: c1Phi03},
				{Key: "phi_030_c2", Val: c2Phi03},
			}},
			{Key: "single_seed_near_c1_r2000", Val: jsonout.Obj{
				{Key: "c_tested", Val: c1Base},
				{Key: "global_cascade_prob", Val: probGlobalSingle},
				{Key: "empirical_ccdf_tail_loglog_slope", Val: ccdfSlope},
			}},
		}},
		{Key: "c3_rogers_veraart_and_stress", Val: jsonout.Obj{
			{Key: "alpha_beta_1_vs_en_max_diff", Val: rv1VsENDiff},
			{Key: "rv_chain_closed_form_max_err", Val: rvChainErr},
			{Key: "monotonicity_in_alpha_beta_verified", Val: rvMonotoneOK},
			{Key: "mc_stress_20k_per_s", Val: stressRows},
			{Key: "poisson_vs_powerlaw_hub_s08", Val: jsonout.Obj{
				{Key: "poisson_random20_en_default_frac", Val: poisRndENDef},
				{Key: "poisson_random20_rv_default_frac", Val: poisRndRVDef},
				{Key: "poisson_random20_en_loss", Val: poisRndENLoss},
				{Key: "poisson_random20_rv_loss", Val: poisRndRVLoss},
				{Key: "poisson_top20_hub_en_default_frac", Val: poisHubENDef},
				{Key: "poisson_top20_hub_rv_default_frac", Val: poisHubRVDef},
				{Key: "poisson_top20_hub_en_loss", Val: poisHubENLoss},
				{Key: "poisson_top20_hub_rv_loss", Val: poisHubRVLoss},
				{Key: "powerlaw_random20_en_default_frac", Val: plRndENDef},
				{Key: "powerlaw_random20_rv_default_frac", Val: plRndRVDef},
				{Key: "powerlaw_random20_en_loss", Val: plRndENLoss},
				{Key: "powerlaw_random20_rv_loss", Val: plRndRVLoss},
				{Key: "powerlaw_top20_hub_en_default_frac", Val: plHubENDef},
				{Key: "powerlaw_top20_hub_rv_default_frac", Val: plHubRVDef},
				{Key: "powerlaw_top20_hub_en_loss", Val: plHubENLoss},
				{Key: "powerlaw_top20_hub_rv_loss", Val: plHubRVLoss},
			}},
		}},
		{Key: "timings_seconds", Val: jsonout.Obj{
			{Key: "c1_clearing", Val: timeC1},
			{Key: "c2_watts_cascade", Val: timeC2},
			{Key: "c3_rv_and_stress", Val: timeC3},
			{Key: "total", Val: timeTotal},
		}},
	}

	if err := jsonout.WriteFile(filepath.Join(outDir, "results.json"), rootObj); err != nil {
		panic(err)
	}

	fmt.Printf("day003c selesai dalam %.3fs (workers=%d): C1=%.3fs, C2=%.3fs, C3=%.3fs\n",
		timeTotal, nWorkers, timeC1, timeC2, timeC3)
	fmt.Printf("  C1: 210 small nets maxDiff PF=%.2e PB=%.2e; chainA=%.2e chainB=%.2e cycle2=%.2e\n",
		maxDiffPFSmall, maxDiffPBSmall, chainErrA, chainErrB, cycle2Err)
	fmt.Printf("  C1 n=2000 Poisson: Picard=%d it (%.4fs) vs Fict=%d rnd (%.4fs) diff=%.2e\n",
		resPoisPicard.Iterations, poisPicardSec, resPoisFict.Iterations, poisFictSec, poisAlgoDiff)
	fmt.Printf("  C2: FIFO==LIFO: %v (%d); phase phi=0.18: c1=%.4f c2=%.4f (sim jump [%.2f, %.2f]); single-seed globalP=%.4f slope=%.3f\n",
		orderIdentical, cntFIFO, c1Base, c2Base, simJumpC1, simJumpC2, probGlobalSingle, ccdfSlope)
	fmt.Printf("  C3: RV(1,1)-EN=%.2e, RVchain=%.2e; Poisson vs PL hub RV loss: Pois(rnd=%.2f, hub=%.2f) PL(rnd=%.2f, hub=%.2f)\n",
		rv1VsENDiff, rvChainErr, poisRndRVLoss, poisHubRVLoss, plRndRVLoss, plHubRVLoss)
}
