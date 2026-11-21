package process

import (
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/closedform"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/parallel"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/stats"
)

// MomentsRow merangkum momen teoritis vs empiris proses surplus U_t di waktu t.
type MomentsRow struct {
	DistCode  string
	DistName  string
	U0        float64
	TimeT     float64
	TheoEN    float64
	SimEN     float64
	SimENSE   float64
	ZN        float64
	TheoVarN  float64
	SimVarN   float64
	TheoEUt   float64
	SimEUt    float64
	SimEUtSE  float64
	ZUt       float64
	TheoVarSt float64
	SimVarSt  float64
	VarStDiff float64
}

// VerifyMoments menyimulasikan N lintasan sampai waktu t untuk memverifikasi momen N_t, S_t, dan U_t.
func VerifyMoments(dists []claims.Dist, u0, tHorizon float64, nPaths, nChunks, nWorkers int, seed uint64) []MomentsRow {
	rows := make([]MomentsRow, len(dists))
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks

	for di, d := range dists {
		_, ex2, _ := d.Moments()
		lamT := claims.Lambda * tHorizon
		theoEUt := u0 + claims.Theta*claims.Lambda*claims.Mu*tHorizon
		theoVarSt := claims.Lambda * ex2 * tHorizon

		type chunkRes struct {
			sN, sqN float64
			sS, sqS float64
			sU, sqU float64
			count   int
		}
		res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkRes {
			cnt := perChunk
			if c < rem {
				cnt++
			}
			rg := rng.New(rng.DeriveChunkSeed(seed, 0x3100^uint64(di), uint64(c)))
			var cr chunkRes
			cr.count = cnt
			for k := 0; k < cnt; k++ {
				nClaims := rg.Poisson(lamT)
				sVal := 0.0
				for j := 0; j < nClaims; j++ {
					sVal += d.Sample(rg)
				}
				uVal := u0 + claims.Premium*tHorizon - sVal
				fn := float64(nClaims)
				cr.sN += fn
				cr.sqN += fn * fn
				cr.sS += sVal
				cr.sqS += sVal * sVal
				cr.sU += uVal
				cr.sqU += uVal * uVal
			}
			return cr
		})

		sN, sqN, sS, sqS, sU, sqU := 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
		for _, r := range res {
			sN += r.sN
			sqN += r.sqN
			sS += r.sS
			sqS += r.sqS
			sU += r.sU
			sqU += r.sqU
		}
		meanN, seN := stats.SummaryMoments(sN, sqN, nPaths)
		varN := (sqN - sN*sN/float64(nPaths)) / float64(nPaths-1)
		meanU, seU := stats.SummaryMoments(sU, sqU, nPaths)
		varSt := (sqS - sS*sS/float64(nPaths)) / float64(nPaths-1)

		rows[di] = MomentsRow{
			DistCode:  d.Code,
			DistName:  d.Name,
			U0:        u0,
			TimeT:     tHorizon,
			TheoEN:    lamT,
			SimEN:     meanN,
			SimENSE:   seN,
			ZN:        math.Abs(meanN-lamT) / seN,
			TheoVarN:  lamT,
			SimVarN:   varN,
			TheoEUt:   theoEUt,
			SimEUt:    meanU,
			SimEUtSE:  seU,
			ZUt:       math.Abs(meanU-theoEUt) / seU,
			TheoVarSt: theoVarSt,
			SimVarSt:  varSt,
			VarStDiff: (varSt - theoVarSt) / theoVarSt,
		}
	}
	return rows
}

// DeficitRow merangkum uji sifat tanpa-memori defisit saat ruin (Exp(mean=1)) beserta uji KS nol-simulasi.
type DeficitRow struct {
	U             float64
	PsiExact      float64
	PsiMC         float64
	PsiMCSE       float64
	PsiZScore     float64
	NRuined       int
	DeficitMean   float64
	DeficitMeanSE float64
	DeficitSD     float64
	KSStat        float64
	KSCrit95      float64
	KSPValue      float64
	KSPass        bool
}

// VerifyExponentialDeficit menguji distribusi bersyarat |U_tau| | tau < inf ~ Exp(1) untuk klaim Exp(1)
// lewat MC dan uji KS satu-sampel menggunakan nilai kritis simulasi nol sendiri.
func VerifyExponentialDeficit(us []float64, nPaths, ksNullReps, nChunks, nWorkers int, seed uint64) []DeficitRow {
	expDist := claims.Exponential()
	cf, _ := closedform.Derive(expDist)
	const stopTime = 180.0
	const uCeiling = 110.0 // sisa e^{-(1/6)*110} ~ 1.1e-8

	rows := make([]DeficitRow, len(us))
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks

	for ui, u := range us {
		type chunkDef struct {
			deficits []float64
			ruined   int
		}
		res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkDef {
			cnt := perChunk
			if c < rem {
				cnt++
			}
			rg := rng.New(rng.DeriveChunkSeed(seed, 0x3200^uint64(ui), uint64(c)))
			var cd chunkDef
			for k := 0; k < cnt; k++ {
				t := 0.0
				surplus := u
				for {
					interT := rg.ExpRate(claims.Lambda)
					t += interT
					surplus += claims.Premium*interT - expDist.Sample(rg)
					if surplus < 0 {
						cd.ruined++
						cd.deficits = append(cd.deficits, -surplus)
						break
					}
					if t > stopTime && surplus > uCeiling {
						break
					}
				}
			}
			return cd
		})

		allDef := make([]float64, 0, nPaths/2)
		totRuined := 0
		for _, r := range res {
			totRuined += r.ruined
			allDef = append(allDef, r.deficits...)
		}
		psiMC := float64(totRuined) / float64(nPaths)
		psiSE := math.Sqrt(psiMC * (1.0 - psiMC) / float64(nPaths))
		psiEx := cf.Eval(u)

		meanDef := stats.Mean(allDef)
		sdDef := stats.SD(allDef)
		seDef := stats.SE(allDef)

		ksStat := stats.KSOneSampleStat(allDef, func(y float64) float64 {
			if y <= 0 {
				return 0
			}
			return 1.0 - math.Exp(-y)
		})
		// Simulasikan distribusi nol KS persis pada jumlah sampel n = len(allDef).
		crit, pVal := stats.KSEvalAgainstSimulatedNull(ksStat, len(allDef), ksNullReps, 0.05, nChunks, nWorkers, seed^uint64(ui*97+1))

		rows[ui] = DeficitRow{
			U:             u,
			PsiExact:      psiEx,
			PsiMC:         psiMC,
			PsiMCSE:       psiSE,
			PsiZScore:     math.Abs(psiMC-psiEx) / psiSE,
			NRuined:       totRuined,
			DeficitMean:   meanDef,
			DeficitMeanSE: seDef,
			DeficitSD:     sdDef,
			KSStat:        ksStat,
			KSCrit95:      crit,
			KSPValue:      pVal,
			KSPass:        ksStat < crit,
		}
	}
	return rows
}

// CrudeMCRow merangkum probabilitas ruin Monte Carlo biasa di T=50, T=200, dan T=inf.
type CrudeMCRow struct {
	DistCode string
	DistName string
	U        float64
	PsiExact float64
	PsiT50   float64
	PsiT200  float64
	PsiInfMC float64
	PsiInfSE float64
	ZScore   float64
	RelSE    float64
	ElapsedS float64
}

// SimulateCrudeMC menyimulasikan N = nPaths lintasan proses surplus klaim L_k = S_{T_k} - c T_k
// secara simultan untuk daftar level modal uList pada horizon T=50, T=200, dan T=inf.
// Aturan henti eksak untuk T=inf menghentikan lintasan setelah T > 200 saat L_max yang dibutuhkan > 150
// di atas L_k sekarang (peluang residual batas atas Lundberg e^{-R * 150} < 1.4e-11 untuk R >= 0.16).
func SimulateCrudeMC(dists []claims.Dist, uList []float64, nPaths, nChunks, nWorkers int, seed uint64) []CrudeMCRow {
	rows := make([]CrudeMCRow, 0, len(dists)*len(uList))
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks

	uMax := 0.0
	for _, u := range uList {
		if u > uMax {
			uMax = u
		}
	}
	for di, d := range dists {
		cf, errCF := closedform.Derive(d)
		hasCF := errCF == nil
		safetyBuffer := 150.0
		if d.HeavyTail {
			safetyBuffer = 2500.0
		}
		type chunkMC struct {
			cntT50  []int
			cntT200 []int
			cntInf  []int
		}
		res := parallel.RunChunks(nChunks, nWorkers, func(c int) chunkMC {
			cnt := perChunk
			if c < rem {
				cnt++
			}
			cm := chunkMC{
				cntT50:  make([]int, len(uList)),
				cntT200: make([]int, len(uList)),
				cntInf:  make([]int, len(uList)),
			}
			rg := rng.New(rng.DeriveChunkSeed(seed, 0x3300^uint64(di), uint64(c)))
			for k := 0; k < cnt; k++ {
				t := 0.0
				loss := 0.0
				max50 := 0.0
				max200 := 0.0
				maxInf := 0.0
				recorded50 := false
				recorded200 := false
				for {
					interT := rg.ExpRate(claims.Lambda)
					tNext := t + interT
					if !recorded50 && tNext > 50.0 {
						max50 = maxInf
						recorded50 = true
					}
					if !recorded200 && tNext > 200.0 {
						max200 = maxInf
						recorded200 = true
					}
					t = tNext
					loss += d.Sample(rg) - claims.Premium*interT
					if loss > maxInf {
						maxInf = loss
					}
					if !recorded50 && t > 50.0 {
						max50 = maxInf
						recorded50 = true
					}
					if !recorded200 && t > 200.0 {
						max200 = maxInf
						recorded200 = true
					}
					// Kondisi penghentian aman setelah horizon 200 terlewati:
					if recorded200 && (loss < -safetyBuffer) && (maxInf > uMax || uMax-loss > safetyBuffer) {
						break
					}
					// Bila sudah ruin untuk semua u in uList dan melewati T=200:
					if recorded200 && maxInf > uMax {
						break
					}
				}
				for ui, u := range uList {
					if max50 > u {
						cm.cntT50[ui]++
					}
					if max200 > u {
						cm.cntT200[ui]++
					}
					if maxInf > u {
						cm.cntInf[ui]++
					}
				}
			}
			return cm
		})

		for ui, u := range uList {
			c50, c200, cInf := 0, 0, 0
			for _, r := range res {
				c50 += r.cntT50[ui]
				c200 += r.cntT200[ui]
				cInf += r.cntInf[ui]
			}
			p50 := float64(c50) / float64(nPaths)
			p200 := float64(c200) / float64(nPaths)
			pInf := float64(cInf) / float64(nPaths)
			seInf := math.Sqrt(pInf * (1.0 - pInf) / float64(nPaths))
			psiEx := math.NaN()
			z := math.NaN()
			if hasCF {
				psiEx = cf.Eval(u)
				if seInf > 0 {
					z = math.Abs(pInf-psiEx) / seInf
				}
			}
			relSE := 0.0
			if pInf > 0 {
				relSE = seInf / pInf
			}
			rows = append(rows, CrudeMCRow{
				DistCode: d.Code,
				DistName: d.Name,
				U:        u,
				PsiExact: psiEx,
				PsiT50:   p50,
				PsiT200:  p200,
				PsiInfMC: pInf,
				PsiInfSE: seInf,
				ZScore:   z,
				RelSE:    relSE,
			})
		}
	}
	return rows
}
