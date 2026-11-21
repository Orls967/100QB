package is

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/closedform"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/parallel"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/stats"
)

// Row merangkum hasil estimator Importance Sampling eksponensial dengan eksponen tilt rTilt.
type Row struct {
	DistCode  string
	DistName  string
	U         float64
	RTilt     float64
	TiltRatio float64
	PsiExact  float64
	PsiIS     float64
	SEIS      float64
	RelSEIS   float64
	ZScore    float64
	SDIS      float64
	VRFTheory float64
	ElapsedS  float64
}

// Simulate menjalankan estimator Importance Sampling eksponensial tak-bias
// pada ukuran tilted: laju kedatangan lambda * M_X(rTilt), klaim ~ e^{rTilt x} f(x)/M_X(rTilt),
// premi c tetap sama. Apabila drift c - lambda * M_X'(rTilt) < 0, ruin terjadi hampir pasti (tau < inf).
func Simulate(d claims.Dist, u, rTilt, tiltRatio float64, nPaths, nChunks, nWorkers int, seed uint64) Row {
	t0 := time.Now()
	mgfR, _, ok := d.MGFAndDeriv(rTilt)
	if !ok {
		return Row{DistCode: d.Code, DistName: d.Name, U: u, RTilt: rTilt, TiltRatio: tiltRatio}
	}
	lamTilt := claims.Lambda * mgfR
	driftFactor := claims.Lambda*(mgfR-1.0) - claims.Premium*rTilt

	cf, _ := closedform.Derive(d)
	psiExact := cf.Eval(u)

	type chunkRes struct {
		sum, sumSq float64
		count      int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	salt := uint64(0x4100) ^ uint64(math.Round(u*10.0)) ^ uint64(math.Round(tiltRatio*1000.0)) ^ uint64(d.Kind*997)

	res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkRes {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, salt, uint64(c)))
		var cr chunkRes
		cr.count = cnt
		for k := 0; k < cnt; k++ {
			t := 0.0
			surplus := u
			for step := 0; step < 5_000_000; step++ {
				interT := rg.ExpRate(lamTilt)
				t += interT
				surplus += claims.Premium*interT - d.SampleTilted(rTilt, rg)
				if surplus < 0 {
					w := math.Exp(-rTilt*(u-surplus) + driftFactor*t)
					cr.sum += w
					cr.sumSq += w * w
					break
				}
			}
		}
		return cr
	})

	sum, sumSq := 0.0, 0.0
	for _, r := range res {
		sum += r.sum
		sumSq += r.sumSq
	}
	meanIS, seIS := stats.SummaryMoments(sum, sumSq, nPaths)
	relSE := 0.0
	if meanIS > 0 {
		relSE = seIS / meanIS
	}
	z := 0.0
	if seIS > 0 {
		z = math.Abs(meanIS-psiExact) / seIS
	}
	sdIS := seIS * math.Sqrt(float64(nPaths))
	varCrude := psiExact * (1.0 - psiExact)
	vrf := 0.0
	if sdIS > 0 {
		vrf = varCrude / (sdIS * sdIS)
	}

	return Row{
		DistCode:  d.Code,
		DistName:  d.Name,
		U:         u,
		RTilt:     rTilt,
		TiltRatio: tiltRatio,
		PsiExact:  psiExact,
		PsiIS:     meanIS,
		SEIS:      seIS,
		RelSEIS:   relSE,
		ZScore:    z,
		SDIS:      sdIS,
		VRFTheory: vrf,
		ElapsedS:  time.Since(t0).Seconds(),
	}
}

// HeavyTailDivergenceRow menunjukkan divergensi integral parsial MGF e^{r x} f(x) hingga batas xMax
// untuk Pareto dan Weibull sebagai bukti kegagalan Importance Sampling eksponensial pada ekor subeksponensial.
type HeavyTailDivergenceRow struct {
	DistCode   string
	DistName   string
	RTest      float64
	XMax       float64
	PartialMGF float64
}

// VerifyHeavyTailMGFFailure menghitung integral terpotong integral_0^{xMax} e^{r x} f(x) dx
// via aturan Simpson majemuk (log-transformasi) untuk memperlihatkan ledakan M_X(r) -> inf pada r > 0.
func VerifyHeavyTailMGFFailure(dists []claims.Dist, rVal float64, xMaxList []float64) []HeavyTailDivergenceRow {
	rows := make([]HeavyTailDivergenceRow, 0, len(dists)*len(xMaxList))
	for _, d := range dists {
		for _, xMax := range xMaxList {
			// Transformasi x = z^2 menyederhanakan singularitas x^{-1/2} Weibull di x = 0.
			const nPanels = 4000
			zMax := math.Sqrt(xMax)
			h := zMax / float64(nPanels)
			acc := 0.0
			for i := 0; i <= nPanels; i++ {
				z := float64(i) * h
				x := z * z
				// dx = 2*z dz; pada Weibull f(x) * 2z = sqrt(2) * e^{-sqrt(2)*z} terbatas baik di z=0!
				integrand := math.Exp(rVal*x) * 2.0 * z * d.PDF(x)
				if i == 0 && d.Kind == claims.KindWeibull {
					integrand = math.Sqrt(1.0 / d.Scale)
				}
				w := 1.0
				if i > 0 && i < nPanels {
					if i%2 == 1 {
						w = 4.0
					} else {
						w = 2.0
					}
				}
				acc += w * integrand
			}
			partial := (h / 3.0) * acc
			rows = append(rows, HeavyTailDivergenceRow{
				DistCode:   d.Code,
				DistName:   d.Name,
				RTest:      rVal,
				XMax:       xMax,
				PartialMGF: partial,
			})
		}
	}
	return rows
}
