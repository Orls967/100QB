package graphs

import (
	"math"
	"sort"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
)

// InEdge menyatakan kewajiban proporsional Pi_{j->i} dari debitur From=j kepada kreditur i.
type InEdge struct {
	From int
	Pi   float64
}

// InterbankNet menyimpan jaringan kewajiban antarbank berarah dengan struktur sparse dan dense opsional.
type InterbankNet struct {
	N      int
	E      []float64 // aset eksternal e_i
	PBar   []float64 // total kewajiban p_bar_i = sum_j L_ij
	In     [][]InEdge
	OutDeg []int
}

// NewFromDense membangun InterbankNet dari aset eksternal e dan matriks kewajiban L (n x n).
func NewFromDense(e []float64, l [][]float64) *InterbankNet {
	n := len(e)
	eCopy := make([]float64, n)
	copy(eCopy, e)
	pBar := make([]float64, n)
	in := make([][]InEdge, n)
	outDeg := make([]int, n)

	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			if i != j && l[i][j] > 0 {
				sum += l[i][j]
				outDeg[i]++
			}
		}
		pBar[i] = sum
	}
	for i := 0; i < n; i++ {
		if pBar[i] <= 0 {
			continue
		}
		for j := 0; j < n; j++ {
			if i != j && l[i][j] > 0 {
				in[j] = append(in[j], InEdge{From: i, Pi: l[i][j] / pBar[i]})
			}
		}
	}
	return &InterbankNet{N: n, E: eCopy, PBar: pBar, In: in, OutDeg: outDeg}
}

// CloneWithE membuat salinan ringan jaringan dengan vektor aset eksternal baru.
func (net *InterbankNet) CloneWithE(eNew []float64) *InterbankNet {
	ec := make([]float64, net.N)
	copy(ec, eNew)
	return &InterbankNet{
		N:      net.N,
		E:      ec,
		PBar:   net.PBar,
		In:     net.In,
		OutDeg: net.OutDeg,
	}
}

// BuildSmallRandomNetwork membangkitkan jaringan acak kecil (siklik/berlapis/tidak-teratur) untuk uji brute-force.
func BuildSmallRandomNetwork(n int, r *rng.RNG, mode int) *InterbankNet {
	e := make([]float64, n)
	l := make([][]float64, n)
	for i := 0; i < n; i++ {
		l[i] = make([]float64, n)
		e[i] = r.Float64() * 1.5
		if mode == 2 && i%3 == 0 {
			e[i] = 0.0 // beberapa bank tanpa aset eksternal
		}
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			var prob float64
			switch mode {
			case 0: // umum / siklik
				prob = 0.45
			case 1: // berlapis (DAG + sedikit umpan balik)
				if j > i {
					prob = 0.65
				} else {
					prob = 0.1
				}
			default: // jarang / tidak-teratur
				prob = 0.25
			}
			if r.Float64() < prob {
				l[i][j] = math.Exp(0.5 * r.NormFloat64())
			}
		}
	}
	return NewFromDense(e, l)
}

// BuildDirectedPoissonNetwork membangkitkan jaringan berarah n bank dengan derajat keluar Poisson(meanDeg)
// dan bobot kewajiban lognormal exp(0.5*Z).
func BuildDirectedPoissonNetwork(n int, meanDeg float64, r *rng.RNG) *InterbankNet {
	pBar := make([]float64, n)
	in := make([][]InEdge, n)
	outDeg := make([]int, n)

	for i := 0; i < n; i++ {
		deg := r.Poisson(meanDeg)
		if deg < 1 {
			deg = 1
		}
		if deg >= n {
			deg = n - 1
		}
		outDeg[i] = deg
		targets := sampleDistinctTargets(n, i, deg, r)
		weights := make([]float64, deg)
		sumW := 0.0
		for k := 0; k < deg; k++ {
			w := math.Exp(0.5 * r.NormFloat64())
			weights[k] = w
			sumW += w
		}
		pBar[i] = sumW
		for k, t := range targets {
			in[t] = append(in[t], InEdge{From: i, Pi: weights[k] / sumW})
		}
	}

	e := make([]float64, n)
	for i := 0; i < n; i++ {
		z := r.NormFloat64()
		e[i] = 0.1 * pBar[i] * math.Exp(0.5*z-0.125)
	}
	return &InterbankNet{N: n, E: e, PBar: pBar, In: in, OutDeg: outDeg}
}

// BuildDirectedPowerLawNetwork membangkitkan jaringan berarah pangkat-hukum P(k) ~ k^{-gamma}
// lewat pencocokan kuantil derajat (configuration-style) dengan rata-rata derajat ~ meanDeg.
func BuildDirectedPowerLawNetwork(n int, gamma, meanDeg float64, r *rng.RNG) *InterbankNet {
	pBar := make([]float64, n)
	in := make([][]InEdge, n)
	outDeg := make([]int, n)

	kMin := 2.0
	kMax := 120.0
	p1 := math.Pow(kMin, 1.0-gamma)
	p2 := math.Pow(kMax, 1.0-gamma)

	rawDeg := make([]float64, n)
	sumRaw := 0.0
	for i := 0; i < n; i++ {
		u := r.Float64()
		kVal := math.Pow(p1+u*(p2-p1), 1.0/(1.0-gamma))
		rawDeg[i] = kVal
		sumRaw += kVal
	}
	scale := (meanDeg * float64(n)) / sumRaw

	for i := 0; i < n; i++ {
		deg := int(math.Round(rawDeg[i] * scale))
		if deg < 1 {
			deg = 1
		}
		if deg >= n {
			deg = n - 1
		}
		outDeg[i] = deg
		targets := sampleDistinctTargets(n, i, deg, r)
		weights := make([]float64, deg)
		sumW := 0.0
		for k := 0; k < deg; k++ {
			w := math.Exp(0.5 * r.NormFloat64())
			weights[k] = w
			sumW += w
		}
		pBar[i] = sumW
		for k, t := range targets {
			in[t] = append(in[t], InEdge{From: i, Pi: weights[k] / sumW})
		}
	}

	e := make([]float64, n)
	for i := 0; i < n; i++ {
		z := r.NormFloat64()
		e[i] = 0.1 * pBar[i] * math.Exp(0.5*z-0.125)
	}
	return &InterbankNet{N: n, E: e, PBar: pBar, In: in, OutDeg: outDeg}
}

// TopHubs mengembalikan indeks k bank dengan derajat total (keluar + masuk) tertinggi.
func (net *InterbankNet) TopHubs(k int) []int {
	type pair struct {
		idx int
		deg int
	}
	arr := make([]pair, net.N)
	for i := 0; i < net.N; i++ {
		arr[i] = pair{idx: i, deg: net.OutDeg[i] + len(net.In[i])}
	}
	sort.Slice(arr, func(a, b int) bool {
		if arr[a].deg != arr[b].deg {
			return arr[a].deg > arr[b].deg
		}
		return arr[a].idx < arr[b].idx
	})
	if k > net.N {
		k = net.N
	}
	out := make([]int, k)
	for i := 0; i < k; i++ {
		out[i] = arr[i].idx
	}
	return out
}

func sampleDistinctTargets(n, self, k int, r *rng.RNG) []int {
	out := make([]int, 0, k)
	used := make(map[int]struct{}, k)
	for len(out) < k {
		cand := r.Intn(n)
		if cand == self {
			continue
		}
		if _, ok := used[cand]; !ok {
			used[cand] = struct{}{}
			out = append(out, cand)
		}
	}
	return out
}

// UndirectedCSR menyimpan graf tak berarah jarang dalam format Compressed Sparse Row datar (ramah cache, nol alokasi per simpul).
type UndirectedCSR struct {
	N      int
	Offset []int32
	Edges  []int32
}

// BuildUndirectedER membangkitkan graf Erdos-Renyi tak berarah G(n, c/n) dalam waktu O(n + m)
// menggunakan metode Batagelj-Brandes (lompatan geometrik).
func BuildUndirectedER(n int, c float64, r *rng.RNG) *UndirectedCSR {
	p := c / float64(n)
	if p <= 0 {
		return &UndirectedCSR{N: n, Offset: make([]int32, n+1)}
	}
	log1mp := math.Log(1.0 - p)
	estM := int(float64(n)*c*0.6) + 64
	uList := make([]int32, 0, estM)
	vList := make([]int32, 0, estM)
	deg := make([]int32, n)

	v := 1
	w := -1
	for v < n {
		u := r.Float64()
		if u <= 1e-15 {
			u = 1e-15
		}
		skip := int(math.Floor(math.Log(u) / log1mp))
		w += 1 + skip
		for w >= v && v < n {
			w -= v
			v++
		}
		if v < n {
			uList = append(uList, int32(v))
			vList = append(vList, int32(w))
			deg[v]++
			deg[w]++
		}
	}

	offset := make([]int32, n+1)
	for i := 0; i < n; i++ {
		offset[i+1] = offset[i] + deg[i]
	}
	edges := make([]int32, offset[n])
	cur := make([]int32, n)
	copy(cur, offset[:n])
	for idx := range uList {
		a, b := uList[idx], vList[idx]
		edges[cur[a]] = b
		cur[a]++
		edges[cur[b]] = a
		cur[b]++
	}
	return &UndirectedCSR{N: n, Offset: offset, Edges: edges}
}
