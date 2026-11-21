package ak

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/panjer"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/parallel"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/stats"
)

// Row merangkum estimasi Asmussen-Kroese bersyarat dibandingkan batas Panjer dan asimtotik subeksponensial.
type Row struct {
	DistCode   string
	DistName   string
	U          float64
	PsiAK      float64
	SEAK       float64
	RelSEAK    float64
	PanjerLow  float64
	PanjerMid  float64
	PanjerHigh float64
	WithinBand bool
	ZVsMid     float64
	Asymptote  float64
	AKOverAsym float64
	SD         float64
	VRFTheory  float64
	ElapsedS   float64
}

func sampleGeom0(rho, logRho float64, rg *rng.RNG) int {
	return int(math.Floor(math.Log(rg.Float64Open()) / logRho))
}

// Simulate menjalankan estimator Asmussen-Kroese berbobot-ukuran (size-biased Pollaczek-Khinchine):
//
//	N = 1 + G1 + G2 dengan G1, G2 iid ~ Geom0(rho) sehingga P(N=n) = n (1-rho)^2 rho^{n-1} (n >= 1),
//	psi_hat_AK(u) = (rho / (1 - rho)) * F_bar_L(max(M_{N-1}, u - S_{N-1})).
//
// Terjamin terbatas dalam [0, rho/(1-rho)] dengan nol sampel sia-sia bernilai K=0.
func Simulate(d claims.Dist, u float64, pSer panjer.Series, nReps, nChunks, nWorkers int, seed uint64) Row {
	t0 := time.Now()
	rho := claims.Rho
	logRho := math.Log(rho)
	scaleAK := rho / (1.0 - rho) // 5.0

	type chunkAK struct {
		sum, sumSq float64
	}
	perChunk := nReps / nChunks
	rem := nReps % nChunks
	salt := uint64(0x5100) ^ uint64(math.Round(u*10.0)) ^ uint64(d.Kind*1337)

	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkAK {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, salt, uint64(c)))
		var ca chunkAK
		for k := 0; k < cnt; k++ {
			nMinus1 := sampleGeom0(rho, logRho, rg) + sampleGeom0(rho, logRho, rg)
			sPrev := 0.0
			mPrev := 0.0
			for j := 0; j < nMinus1; j++ {
				y := d.SampleLadder(rg)
				sPrev += y
				if y > mPrev {
					mPrev = y
				}
			}
			arg := u - sPrev
			if mPrev > arg {
				arg = mPrev
			}
			est := scaleAK * d.LadderSurvival(arg)
			ca.sum += est
			ca.sumSq += est * est
		}
		return ca
	})

	sum, sumSq := 0.0, 0.0
	for _, r := range res {
		sum += r.sum
		sumSq += r.sumSq
	}
	meanAK, seAK := stats.SummaryMoments(sum, sumSq, nReps)
	relSE := 0.0
	if meanAK > 0 {
		relSE = seAK / meanAK
	}
	lo, mid, hi, _, _, _ := pSer.At(u)
	tolBand := (hi-lo)*0.5 + 3.0*seAK
	within := math.Abs(meanAK-mid) <= tolBand
	zMid := 0.0
	if seAK > 0 {
		zMid = math.Abs(meanAK-mid) / seAK
	}
	asym := d.SubexpAsymptotic(u)
	ratioAs := 0.0
	if asym > 0 {
		ratioAs = meanAK / asym
	}
	sd := seAK * math.Sqrt(float64(nReps))
	varCrude := mid * (1.0 - mid)
	vrf := 0.0
	if sd > 0 {
		vrf = varCrude / (sd * sd)
	}

	return Row{
		DistCode:   d.Code,
		DistName:   d.Name,
		U:          u,
		PsiAK:      meanAK,
		SEAK:       seAK,
		RelSEAK:    relSE,
		PanjerLow:  lo,
		PanjerMid:  mid,
		PanjerHigh: hi,
		WithinBand: within,
		ZVsMid:     zMid,
		Asymptote:  asym,
		AKOverAsym: ratioAs,
		SD:         sd,
		VRFTheory:  vrf,
		ElapsedS:   time.Since(t0).Seconds(),
	}
}
