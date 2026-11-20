package rv

import (
	"math"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/clearing"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/parallel"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/stats"
)

// Clear menjalankan iterasi titik tetap terbesar Rogers-Veraart (2013) dari p^(0) = p_bar.
// Untuk bank sehat (e_i + sum_j Pi_{ji} p_j >= p_bar_i): Phi_RV_i(p) = p_bar_i.
// Untuk bank gagal (e_i + sum_j Pi_{ji} p_j < p_bar_i):  Phi_RV_i(p) = alpha * e_i + beta * sum_j Pi_{ji} p_j.
func Clear(net *graphs.InterbankNet, alpha, beta, tol float64, maxIter int) clearing.Result {
	n := net.N
	p := make([]float64, n)
	copy(p, net.PBar)
	next := make([]float64, n)

	it := 0
	for it < maxIter {
		it++
		maxDiff := 0.0
		for i := 0; i < n; i++ {
			interbank := 0.0
			for _, ed := range net.In[i] {
				interbank += ed.Pi * p[ed.From]
			}
			totalAsset := net.E[i] + interbank
			var val float64
			if totalAsset >= net.PBar[i]-1e-14 {
				val = net.PBar[i]
			} else {
				val = alpha*net.E[i] + beta*interbank
				if val > net.PBar[i] {
					val = net.PBar[i]
				}
			}
			next[i] = val
			if d := math.Abs(val - p[i]); d > maxDiff {
				maxDiff = d
			}
		}
		copy(p, next)
		if maxDiff <= tol {
			break
		}
	}

	defCnt := 0
	loss := 0.0
	for i := 0; i < n; i++ {
		diff := net.PBar[i] - p[i]
		if diff > 1e-9 {
			defCnt++
		}
		if diff > 0 {
			loss += diff
		}
	}
	return clearing.Result{
		P:            p,
		DefaultCount: defCnt,
		TotalLoss:    loss,
		Iterations:   it,
	}
}

// ChainClosedForm menghitung solusi tertutup rantai 0->1->...->n-1 (L_{i,i+1}=1, p_bar=1)
// ketika setiap bank memiliki e_i kecil sehingga semua gagal: p_0 = alpha*e_0, p_i = alpha*e_i + beta*p_{i-1}.
func ChainClosedForm(e []float64, alpha, beta float64) []float64 {
	n := len(e)
	p := make([]float64, n)
	for i := 0; i < n; i++ {
		inflow := 0.0
		if i > 0 {
			inflow = p[i-1]
		}
		if e[i]+inflow >= 1.0 {
			p[i] = 1.0
		} else {
			p[i] = alpha*e[i] + beta*inflow
		}
	}
	return p
}

// StressSummary merangkum hasil uji stres Monte Carlo pada satu nilai shock s.
type StressSummary struct {
	ShockS               float64
	ENMeanDefaultFrac    float64
	RVMeanDefaultFrac    float64
	ENMeanLoss           float64
	RVMeanLoss           float64
	ENProbDefaultGe10Pct float64
	RVProbDefaultGe10Pct float64
	ENVaR99              float64
	ENES975              float64
	RVVaR99              float64
	RVES975              float64
	AmplificationRatio   float64
	AllRVGeEN            bool
	SampleENLosses       []float64
	SampleRVLosses       []float64
}

type scenarioOut struct {
	enDefFrac float64
	rvDefFrac float64
	enLoss    float64
	rvLoss    float64
}

// RunMonteCarloStress menjalankan nScenarios skenario guncangan e_i <- e_i * (1 - s * u_i), u_i ~ Bernoulli(pShock)
// secara paralel deterministik untuk kliring EN dan RV(alpha, beta).
func RunMonteCarloStress(
	baseNet *graphs.InterbankNet,
	shockS, pShock, alpha, beta float64,
	nScenarios, nWorkers int,
	seed uint64,
) StressSummary {
	n := baseNet.N
	chunks := parallel.RunChunks(nScenarios, 100, nWorkers, seed, func(_ int, start, end int, cSeed uint64) []scenarioOut {
		r := rng.New(cSeed)
		out := make([]scenarioOut, 0, end-start)
		eShock := make([]float64, n)
		for sc := start; sc < end; sc++ {
			for i := 0; i < n; i++ {
				if r.Float64() < pShock {
					eShock[i] = baseNet.E[i] * (1.0 - shockS)
				} else {
					eShock[i] = baseNet.E[i]
				}
			}
			netSc := baseNet.CloneWithE(eShock)
			resEN := clearing.Picard(netSc, 1e-9, 250)
			resRV := Clear(netSc, alpha, beta, 1e-9, 250)
			out = append(out, scenarioOut{
				enDefFrac: float64(resEN.DefaultCount) / float64(n),
				rvDefFrac: float64(resRV.DefaultCount) / float64(n),
				enLoss:    resEN.TotalLoss,
				rvLoss:    resRV.TotalLoss,
			})
		}
		return out
	})

	enDef := make([]float64, 0, nScenarios)
	rvDef := make([]float64, 0, nScenarios)
	enLoss := make([]float64, 0, nScenarios)
	rvLoss := make([]float64, 0, nScenarios)
	enGe10 := 0
	rvGe10 := 0
	allGe := true

	for _, ch := range chunks {
		for _, item := range ch {
			enDef = append(enDef, item.enDefFrac)
			rvDef = append(rvDef, item.rvDefFrac)
			enLoss = append(enLoss, item.enLoss)
			rvLoss = append(rvLoss, item.rvLoss)
			if item.enDefFrac >= 0.10 {
				enGe10++
			}
			if item.rvDefFrac >= 0.10 {
				rvGe10++
			}
			if item.rvLoss+1e-7 < item.enLoss {
				allGe = false
			}
		}
	}

	enVaR, enES := stats.VaRAndES(enLoss, 0.99, 0.975)
	rvVaR, rvES := stats.VaRAndES(rvLoss, 0.99, 0.975)
	mEN := stats.Mean(enLoss)
	mRV := stats.Mean(rvLoss)
	amp := 1.0
	if mEN > 1e-12 {
		amp = mRV / mEN
	}
	sampleCnt := 500
	if len(enLoss) < sampleCnt {
		sampleCnt = len(enLoss)
	}

	return StressSummary{
		ShockS:               shockS,
		ENMeanDefaultFrac:    stats.Mean(enDef),
		RVMeanDefaultFrac:    stats.Mean(rvDef),
		ENMeanLoss:           mEN,
		RVMeanLoss:           mRV,
		ENProbDefaultGe10Pct: float64(enGe10) / float64(nScenarios),
		RVProbDefaultGe10Pct: float64(rvGe10) / float64(nScenarios),
		ENVaR99:              enVaR,
		ENES975:              enES,
		RVVaR99:              rvVaR,
		RVES975:              rvES,
		AmplificationRatio:   amp,
		AllRVGeEN:            allGe,
		SampleENLosses:       enLoss[:sampleCnt],
		SampleRVLosses:       rvLoss[:sampleCnt],
	}
}
