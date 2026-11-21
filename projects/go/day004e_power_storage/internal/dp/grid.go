package dp

import "sort"

// PriceGrid mendefinisikan diskretisasi harga menjadi bin interior dan dua tail bin terbuka:
// bin 0 = (-inf, Edges[0]), bin NumBins-1 = [Edges[len-1], +inf).
type PriceGrid struct {
	Name        string
	NumBins     int
	Edges       []float64 // panjang = NumBins - 1
	RepPrices   []float64 // nilai representatif tiap bin (panjang = NumBins)
	LowerTailAt float64   // batas bawah tail bin bawah (-500)
	UpperTailAt float64   // batas atas tail bin atas (2000)
}

func buildGrid(name string, edges []float64, lowerTailRep, upperTailRep float64) PriceGrid {
	numBins := len(edges) + 1
	reps := make([]float64, numBins)
	reps[0] = lowerTailRep
	for b := 1; b < numBins-1; b++ {
		reps[b] = 0.5 * (edges[b-1] + edges[b])
	}
	reps[numBins-1] = upperTailRep
	cpEdges := make([]float64, len(edges))
	copy(cpEdges, edges)
	return PriceGrid{
		Name:        name,
		NumBins:     numBins,
		Edges:       cpEdges,
		RepPrices:   reps,
		LowerTailAt: edges[0],
		UpperTailAt: edges[len(edges)-1],
	}
}

// NewCoarseGrid17 membangun grid harga kasar 17 bin termasuk tail bin < -500 dan >= 2000 EUR/MWh.
func NewCoarseGrid17() PriceGrid {
	edges := []float64{
		-500, -200, -100, -50, -20, 0, 20, 40, 60, 80, 120, 200, 350, 600, 1000, 2000,
	}
	return buildGrid("coarse_17", edges, -600.0, 2500.0)
}

// NewFineGrid33 membangun grid harga halus 33 bin (perhalusan ketat atas grid 17 bin)
// termasuk tail bin < -500 dan >= 2000 EUR/MWh.
func NewFineGrid33() PriceGrid {
	edges := []float64{
		-500, -350, -200, -150, -100, -75, -50, -35, -20, -10,
		0, 10, 20, 30, 40, 50, 60, 70, 80, 100,
		120, 140, 160, 200, 275, 350, 475, 600, 800, 1000,
		1500, 2000,
	}
	return buildGrid("fine_33", edges, -600.0, 2500.0)
}

// NewExpandedTailGrid membangun grid perluasan 38 bin yang mempertahankan seluruh batas
// interior fine_33 dan memperluas domain menjadi [-1000, 5000] EUR/MWh dengan tail bin eksplisit
// < -1000 dan >= 5000 EUR/MWh.
func NewExpandedTailGrid() PriceGrid {
	edges := []float64{
		-1000, -750,
		-500, -350, -200, -150, -100, -75, -50, -35, -20, -10,
		0, 10, 20, 30, 40, 50, 60, 70, 80, 100,
		120, 140, 160, 200, 275, 350, 475, 600, 800, 1000,
		1500, 2000,
		2500, 3500, 5000,
	}
	return buildGrid("expanded_tail_38", edges, -1150.0, 5500.0)
}

// BinIndex memetakan harga kontinu ke indeks bin [0..NumBins-1] tanpa membuang ekor ekstrem.
func (g PriceGrid) BinIndex(price float64) int {
	if price < g.Edges[0] {
		return 0
	}
	last := len(g.Edges) - 1
	if price >= g.Edges[last] {
		return g.NumBins - 1
	}
	// Cari indeks i terkecil sehingga Edges[i] > price; maka Edges[i-1] <= price < Edges[i].
	idx := sort.Search(len(g.Edges), func(i int) bool {
		return g.Edges[i] > price
	})
	return idx
}
