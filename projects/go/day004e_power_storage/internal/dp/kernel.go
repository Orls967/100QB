package dp

import (
	"time"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/parallel"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/prices"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

// TransitionKernel menyimpan peluang transisi satu jam P(r_next, bin_next | h, r, bin)
// beserta statistik massa probabilitas di tail bin bawah dan atas.
type TransitionKernel struct {
	Grid           PriceGrid
	DrawsPerState  int
	NumStates      int
	TargetStride   int // NumRegimes * Grid.NumBins
	TotalDraws     int64
	Probs          []float64 // panjang NumStates * TargetStride
	LowerTailMass  float64
	UpperTailMass  float64
	LowerTailCount int64
	UpperTailCount int64
	LowerTailMeanP float64
	UpperTailMeanP float64
	BuildSeconds   float64
	WorkersUsed    int
	KernelSeedUsed uint64
}

type stateKernelChunk struct {
	probs          []float64
	lowerTailCount int64
	upperTailCount int64
	lowerTailSumP  float64
	upperTailSumP  float64
}

// StateIndex mengembalikan indeks datar state (h, r, b).
func (k TransitionKernel) StateIndex(h, r, b int) int {
	return (h*regime.NumRegimes+r)*k.Grid.NumBins + b
}

// TargetIndex mengembalikan indeks datar target (rNext, bNext).
func (k TransitionKernel) TargetIndex(rNext, bNext int) int {
	return rNext*k.Grid.NumBins + bNext
}

// BuildKernel mengestimasi kernel transisi satu jam via simulasi kondisional Monte Carlo independen
// per state (h, r, b) dengan seed turunan deterministik tetap.
func BuildKernel(spec prices.ModelSpec, grid PriceGrid, drawsPerState, nWorkers int, kernelSeed, salt uint64) TransitionKernel {
	return BuildKernelWithSeedMapper(spec, grid, drawsPerState, nWorkers, kernelSeed, salt, nil)
}

// BuildKernelWithSeedMapper membangun kernel transisi dengan pemetaan seed opsional per bin
// agar draw kondisional common random numbers (CRN) dapat disejajarkan pada bin interior yang sama.
func BuildKernelWithSeedMapper(
	spec prices.ModelSpec,
	grid PriceGrid,
	drawsPerState, nWorkers int,
	kernelSeed, defaultSalt uint64,
	seedMapper func(h, r, b int) (salt uint64, origBin int),
) TransitionKernel {
	t0 := time.Now()
	nBins := grid.NumBins
	nStates := prices.HoursPerDay * regime.NumRegimes * nBins
	stride := regime.NumRegimes * nBins

	chunks := parallel.RunChunks(nStates, nWorkers, func(stIdx int) stateKernelChunk {
		h := stIdx / (regime.NumRegimes * nBins)
		rem := stIdx % (regime.NumRegimes * nBins)
		r := rem / nBins
		b := rem % nBins

		useSalt := defaultSalt
		useBin := b
		if seedMapper != nil {
			useSalt, useBin = seedMapper(h, r, b)
		}
		rg := rng.New(rng.DeriveStateSeed(kernelSeed, useSalt, uint64(h), uint64(r), uint64(useBin)))
		pNow := grid.RepPrices[b]
		counts := make([]int, stride)
		var loCnt, hiCnt int64
		var loSum, hiSum float64

		for d := 0; d < drawsPerState; d++ {
			rNext, pNext, _, _ := spec.Step(r, h, pNow, rg)
			bNext := grid.BinIndex(pNext)
			counts[rNext*nBins+bNext]++
			if bNext == 0 {
				loCnt++
				loSum += pNext
			} else if bNext == nBins-1 {
				hiCnt++
				hiSum += pNext
			}
		}

		invN := 1.0 / float64(drawsPerState)
		probs := make([]float64, stride)
		for i, c := range counts {
			probs[i] = float64(c) * invN
		}
		return stateKernelChunk{
			probs:          probs,
			lowerTailCount: loCnt,
			upperTailCount: hiCnt,
			lowerTailSumP:  loSum,
			upperTailSumP:  hiSum,
		}
	})

	allProbs := make([]float64, nStates*stride)
	var totLo, totHi int64
	var sumLoP, sumHiP float64
	for stIdx, ch := range chunks {
		copy(allProbs[stIdx*stride:(stIdx+1)*stride], ch.probs)
		totLo += ch.lowerTailCount
		totHi += ch.upperTailCount
		sumLoP += ch.lowerTailSumP
		sumHiP += ch.upperTailSumP
	}
	totDraws := int64(nStates) * int64(drawsPerState)
	loMean := grid.RepPrices[0]
	if totLo > 0 {
		loMean = sumLoP / float64(totLo)
	}
	hiMean := grid.RepPrices[nBins-1]
	if totHi > 0 {
		hiMean = sumHiP / float64(totHi)
	}

	return TransitionKernel{
		Grid:           grid,
		DrawsPerState:  drawsPerState,
		NumStates:      nStates,
		TargetStride:   stride,
		TotalDraws:     totDraws,
		Probs:          allProbs,
		LowerTailMass:  float64(totLo) / float64(totDraws),
		UpperTailMass:  float64(totHi) / float64(totDraws),
		LowerTailCount: totLo,
		UpperTailCount: totHi,
		LowerTailMeanP: loMean,
		UpperTailMeanP: hiMean,
		BuildSeconds:   time.Since(t0).Seconds(),
		WorkersUsed:    nWorkers,
		KernelSeedUsed: kernelSeed,
	}
}
