package kalman

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/afns"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/parallel"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/stats"
)

// ConfigSpec mendefinisikan satu konfigurasi eksperimen MLE AFNS.
type ConfigSpec struct {
	Name              string
	Group             string
	FrequencyLabel    string
	HorizonYears      float64
	Dt                float64
	TSteps            int
	NumReps           int
	Salt              uint64
	SharedPrefixGroup string
	ReusedFrom        string
}

// ParamSummary menyimpan statistik ringkasan lengkap seluruh R replikasi untuk satu parameter.
// Seluruh replikasi dipertahankan tanpa membuang kasus gagal konvergen atau solusi batas.
type ParamSummary struct {
	Name                        string
	TrueValue                   float64
	MeanEstimate                float64
	Bias                        float64
	BiasMCSE                    float64
	BiasInMCSE                  float64
	EmpiricalSD                 float64
	RMSE                        float64
	MeanHessianSE               float64
	MedianHessianSE             float64
	MeanHessianSEToEmpSDRatio   float64
	MedianHessianSEToEmpSDRatio float64
	Coverage95                  float64
	CoverageWilsonLower         float64
	CoverageWilsonUpper         float64
	CoveredCount                int
	TotalReps                   int
	BoundHitCount               int
	BoundHitPct                 float64
}

// StudyConfigResult menyimpan ringkasan agregat dan seluruh hasil per replikasi satu konfigurasi.
type StudyConfigResult struct {
	Spec            ConfigSpec
	NumWorkers      int
	ElapsedSeconds  float64
	ConvergedCount  int
	ConvergedPct    float64
	LLGeTrueCount   int
	LLGeTruePct     float64
	BoundHitCount   int
	BoundHitPct     float64
	HessianRegCount int
	Params          [11]ParamSummary
	Replications    []MLEResult
}

// summarizeReplications menghitung statistik ringkasan seluruh replikasi tanpa eksklusi apa pun.
func summarizeReplications(trueSpec ModelSpec, cfg ConfigSpec, nWorkers int, elapsed float64, reps []MLEResult) StudyConfigResult {
	trueV := SpecToVec(trueSpec)
	nRep := len(reps)
	convCnt, llGeCnt, hitCnt, regCnt := 0, 0, 0, 0
	var estCols [11][]float64
	var hessCols [11][]float64
	var covCnts [11]int
	var hitParamCnts [11]int
	for i := 0; i < 11; i++ {
		estCols[i] = make([]float64, nRep)
		hessCols[i] = make([]float64, nRep)
	}

	for r, res := range reps {
		if res.Converged {
			convCnt++
		}
		if res.EstLogLik >= res.TrueLogLik-1e-8 {
			llGeCnt++
		}
		if res.HitBound {
			hitCnt++
		}
		if res.HessianRegularized {
			regCnt++
		}
		for i := 0; i < 11; i++ {
			est := res.EstParams[i]
			se := res.HessianSE[i]
			estCols[i][r] = est
			hessCols[i][r] = se
			if res.HitBoundByParam[i] {
				hitParamCnts[i]++
			}
			if se > 0 && math.Abs(est-trueV[i]) <= stats.Z95*se {
				covCnts[i]++
			}
		}
	}

	var pSum [11]ParamSummary
	fRep := float64(nRep)
	for i := 0; i < 11; i++ {
		meanEst := stats.Mean(estCols[i])
		bias := meanEst - trueV[i]
		empSD := stats.SD(estCols[i])
		mcse := stats.SE(estCols[i])
		rmse := stats.RMSE(estCols[i], trueV[i])
		meanHSE := stats.Mean(hessCols[i])
		medHSE := stats.Quantile(hessCols[i], 0.5)
		biasZ := 0.0
		if mcse > 0 {
			biasZ = math.Abs(bias) / mcse
		}
		meanRatio, medRatio := 0.0, 0.0
		if empSD > 0 {
			meanRatio = meanHSE / empSD
			medRatio = medHSE / empSD
		}
		cov, wLo, wHi := stats.WilsonCI(covCnts[i], nRep, stats.Z95)
		pSum[i] = ParamSummary{
			Name:                        ParamNames[i],
			TrueValue:                   trueV[i],
			MeanEstimate:                meanEst,
			Bias:                        bias,
			BiasMCSE:                    mcse,
			BiasInMCSE:                  biasZ,
			EmpiricalSD:                 empSD,
			RMSE:                        rmse,
			MeanHessianSE:               meanHSE,
			MedianHessianSE:             medHSE,
			MeanHessianSEToEmpSDRatio:   meanRatio,
			MedianHessianSEToEmpSDRatio: medRatio,
			Coverage95:                  cov,
			CoverageWilsonLower:         wLo,
			CoverageWilsonUpper:         wHi,
			CoveredCount:                covCnts[i],
			TotalReps:                   nRep,
			BoundHitCount:               hitParamCnts[i],
			BoundHitPct:                 100.0 * float64(hitParamCnts[i]) / math.Max(fRep, 1.0),
		}
	}

	return StudyConfigResult{
		Spec:            cfg,
		NumWorkers:      nWorkers,
		ElapsedSeconds:  elapsed,
		ConvergedCount:  convCnt,
		ConvergedPct:    100.0 * float64(convCnt) / math.Max(fRep, 1.0),
		LLGeTrueCount:   llGeCnt,
		LLGeTruePct:     100.0 * float64(llGeCnt) / math.Max(fRep, 1.0),
		BoundHitCount:   hitCnt,
		BoundHitPct:     100.0 * float64(hitCnt) / math.Max(fRep, 1.0),
		HessianRegCount: regCnt,
		Params:          pSum,
		Replications:    reps,
	}
}

// RunIndependentConfigStudy menjalankan satu konfigurasi studi MLE dengan benih independen per replikasi
// sehingga keluaran identik bit-demi-bit untuk 1 atau banyak worker.
func RunIndependentConfigStudy(trueSpec ModelSpec, taus []float64, cfg ConfigSpec, nWorkers int, mainSeed uint64) StudyConfigResult {
	t0 := time.Now()
	reps := parallel.RunChunks(cfg.NumReps, nWorkers, func(repIdx int) MLEResult {
		rg := rng.New(rng.DeriveSeed(mainSeed, cfg.Salt, uint64(repIdx)))
		_, yields := afns.SimulateStationaryPath(trueSpec.AFNS, taus, cfg.TSteps, cfg.Dt, trueSpec.SigmaEps, rg)
		initSpec := PerturbInitialGuess(trueSpec, rg)
		return EstimateMLE(initSpec, trueSpec, taus, cfg.Dt, yields)
	})
	return summarizeReplications(trueSpec, cfg, nWorkers, time.Since(t0).Seconds(), reps)
}

// SharedPrefixPath menyimpan satu lintasan penuh panjang maksimum dan tebakan awal bersama satu replikasi.
type SharedPrefixPath struct {
	Yields   [][]float64
	InitSpec ModelSpec
}

// RunLongSpanSharedPrefixStudy menjalankan studi long-span mingguan untuk daftar horizon kalender (tahun),
// memakai prefiks lintasan dan tebakan awal yang sama per replikasi agar perbandingan antar horizon terisolasi.
func RunLongSpanSharedPrefixStudy(trueSpec ModelSpec, taus []float64, horizonYears []int, nRep, nWorkers int, mainSeed uint64, baseSalt uint64) []StudyConfigResult {
	if len(horizonYears) == 0 {
		return nil
	}
	maxH := 0
	for _, h := range horizonYears {
		if h > maxH {
			maxH = h
		}
	}
	maxSteps := 52 * maxH

	// Bangkitkan lintasan panjang maksimum per replikasi secara deterministik
	paths := parallel.RunChunks(nRep, nWorkers, func(repIdx int) SharedPrefixPath {
		rg := rng.New(rng.DeriveSeed(mainSeed, baseSalt, uint64(repIdx)))
		_, y := afns.SimulateStationaryPath(trueSpec.AFNS, taus, maxSteps, WeeklyDt, trueSpec.SigmaEps, rg)
		initSpec := PerturbInitialGuess(trueSpec, rg)
		return SharedPrefixPath{Yields: y, InitSpec: initSpec}
	})

	results := make([]StudyConfigResult, len(horizonYears))
	for hIdx, hYr := range horizonYears {
		steps := 52 * hYr
		cfg := ConfigSpec{
			Name:              formatLongSpanName(hYr),
			Group:             "longspan",
			FrequencyLabel:    "weekly",
			HorizonYears:      float64(hYr),
			Dt:                WeeklyDt,
			TSteps:            steps,
			NumReps:           nRep,
			Salt:              baseSalt,
			SharedPrefixGroup: "longspan_weekly_80y_prefix",
		}
		t0 := time.Now()
		reps := parallel.RunChunks(nRep, nWorkers, func(repIdx int) MLEResult {
			p := paths[repIdx]
			return EstimateMLE(p.InitSpec, trueSpec, taus, WeeklyDt, p.Yields[:steps])
		})
		results[hIdx] = summarizeReplications(trueSpec, cfg, nWorkers, time.Since(t0).Seconds(), reps)
	}
	return results
}

func formatLongSpanName(hYr int) string {
	return "longspan_" + itoa(hYr) + "y_weekly"
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [16]byte
	pos := len(buf)
	for n > 0 {
		pos--
		buf[pos] = byte('0' + n%10)
		n /= 10
	}
	return string(buf[pos:])
}

// ReuseAsInfillConfig membungkus hasil horizon 20 tahun mingguan dari studi long-span untuk studi infill
// tanpa menjalankan ulang optimisasi identik, sambil mendokumentasikan sumber pemakaian ulangnya.
func ReuseAsInfillConfig(res20y StudyConfigResult) StudyConfigResult {
	cp := res20y
	cp.Spec.Name = "infill_20y_weekly"
	cp.Spec.Group = "infill"
	cp.Spec.ReusedFrom = res20y.Spec.Name
	cp.ElapsedSeconds = 0.0
	return cp
}

// Day004aBaselineRef menyimpan statistik acuan R=50, T=1000 mingguan dari outputs/day004a/results.json.
type Day004aBaselineRef struct {
	Name          string
	TrueValue     float64
	MeanEstimate  float64
	Bias          float64
	RepSE         float64
	EmpiricalSD   float64
	RMSE          float64
	MeanHessianSE float64
	Coverage95    float64
	CoveredCount  int
}

// Day004aReferenceBaseline adalah nilai eksak terukur dari outputs/day004a/results.json (R=50, T=1000 mingguan).
var Day004aReferenceBaseline = [11]Day004aBaselineRef{
	{Name: "lambda", TrueValue: 0.6, MeanEstimate: 0.60004078725797383, Bias: 4.0787257973851077e-05, RepSE: 0.00042673588988184716, EmpiricalSD: 0.0030174784151112996, RMSE: 0.0029874296755509383, MeanHessianSE: 0.002895217250646133, Coverage95: 0.96, CoveredCount: 48},
	{Name: "kp11", TrueValue: 0.05, MeanEstimate: 0.26355171427733082, Bias: 0.21355171427733083, RepSE: 0.028112136192094257, EmpiricalSD: 0.19878282135069616, RMSE: 0.29039396091111586, MeanHessianSE: 0.43862652291416704, Coverage95: 0.80, CoveredCount: 40},
	{Name: "kp22", TrueValue: 0.5, MeanEstimate: 0.7061276499407555, Bias: 0.2061276499407555, RepSE: 0.042351788907655868, EmpiricalSD: 0.29947237131984672, RMSE: 0.36107981836479325, MeanHessianSE: 0.26483748807649637, Coverage95: 0.94, CoveredCount: 47},
	{Name: "kp33", TrueValue: 0.4, MeanEstimate: 0.67069829555246641, Bias: 0.27069829555246638, RepSE: 0.041844592161070907, EmpiricalSD: 0.29588594873078689, RMSE: 0.39884218931577703, MeanHessianSE: 0.25645451769317618, Coverage95: 0.88, CoveredCount: 44},
	{Name: "theta1", TrueValue: 0.045, MeanEstimate: 0.044509533535344675, Bias: -0.00049046646465532306, RepSE: 0.0024644747090346963, EmpiricalSD: 0.017426467788211775, RMSE: 0.017258293697091334, MeanHessianSE: 0.25608631384207753, Coverage95: 0.68, CoveredCount: 34},
	{Name: "theta2", TrueValue: -0.02, MeanEstimate: -0.018985360778742414, Bias: 0.0010146392212575864, RepSE: 0.00057036728995014723, EmpiricalSD: 0.0040331057849074285, RMSE: 0.0041194800856567916, MeanHessianSE: 0.003540069519436364, Coverage95: 0.82, CoveredCount: 41},
	{Name: "theta3", TrueValue: 0.01, MeanEstimate: 0.0083410409474636675, Bias: -0.0016589590525363327, RepSE: 0.0012578566341397467, EmpiricalSD: 0.00889438955760701, RMSE: 0.0089599167087876172, MeanHessianSE: 0.0074504225301608388, Coverage95: 0.80, CoveredCount: 40},
	{Name: "sigma1", TrueValue: 0.006, MeanEstimate: 0.0060030573084533965, Bias: 3.0573084533963915e-06, RepSE: 2.8569452098983222e-06, EmpiricalSD: 2.0201653313975281e-05, RMSE: 2.0230961313415059e-05, MeanHessianSE: 1.8467452132643473e-05, Coverage95: 0.92, CoveredCount: 46},
	{Name: "sigma2", TrueValue: 0.01, MeanEstimate: 0.010020699817819889, Bias: 2.0699817819888958e-05, RepSE: 3.9714476596956914e-05, EmpiricalSD: 0.00028082375712982675, RMSE: 0.00027877091917005671, MeanHessianSE: 0.0002737129910850619, Coverage95: 0.92, CoveredCount: 46},
	{Name: "sigma3", TrueValue: 0.02, MeanEstimate: 0.020184098018196438, Bias: 0.00018409801819643801, RepSE: 0.00010602853538012562, EmpiricalSD: 0.00074973496366564603, RMSE: 0.00076469114401288101, MeanHessianSE: 0.00067317049581555415, Coverage95: 0.90, CoveredCount: 45},
	{Name: "sigma_eps", TrueValue: 0.0005, MeanEstimate: 0.00050084053188757837, Bias: 8.4053188757835916e-07, RepSE: 6.4462256695265759e-07, EmpiricalSD: 4.5581698839810344e-06, RMSE: 4.5899747593470443e-06, MeanHessianSE: 4.149886664089327e-06, Coverage95: 0.96, CoveredCount: 48},
}

// BaselineParamDiff membandingkan bias estimasi baseline day004b terhadap day004a memakai MCSE gabungan.
type BaselineParamDiff struct {
	Name                string
	Bias4b              float64
	MCSE4b              float64
	Bias4a              float64
	MCSE4a              float64
	BiasDiff            float64
	CombinedMCSE        float64
	DiffZScore          float64
	Within2CombinedMCSE bool
	Within3CombinedMCSE bool
}

// BaselineComparisonSummary merangkum perbandingan baseline day004b (1000 minggu) vs day004a.
type BaselineComparisonSummary struct {
	Params         [11]BaselineParamDiff
	MaxAbsZ        float64
	CountWithin2SE int
	CountWithin3SE int
	AllWithin3SE   bool
}

// CompareBaselineToDay004a membandingkan bias replikasi baseline day004b terhadap day004a.
func CompareBaselineToDay004a(baseline StudyConfigResult) BaselineComparisonSummary {
	var out BaselineComparisonSummary
	out.AllWithin3SE = true
	for i := 0; i < 11; i++ {
		p4b := baseline.Params[i]
		p4a := Day004aReferenceBaseline[i]
		diff := p4b.Bias - p4a.Bias
		combSE := math.Sqrt(p4b.BiasMCSE*p4b.BiasMCSE + p4a.RepSE*p4a.RepSE)
		z := 0.0
		if combSE > 0 {
			z = diff / combSE
		}
		w2 := math.Abs(z) <= 2.0
		w3 := math.Abs(z) <= 3.0
		if w2 {
			out.CountWithin2SE++
		}
		if w3 {
			out.CountWithin3SE++
		} else {
			out.AllWithin3SE = false
		}
		if math.Abs(z) > out.MaxAbsZ {
			out.MaxAbsZ = math.Abs(z)
		}
		out.Params[i] = BaselineParamDiff{
			Name:                p4b.Name,
			Bias4b:              p4b.Bias,
			MCSE4b:              p4b.BiasMCSE,
			Bias4a:              p4a.Bias,
			MCSE4a:              p4a.RepSE,
			BiasDiff:            diff,
			CombinedMCSE:        combSE,
			DiffZScore:          z,
			Within2CombinedMCSE: w2,
			Within3CombinedMCSE: w3,
		}
	}
	return out
}

// InvHorizonTrendSummary menganalisis hubungan bias terhadap 1/H beserta ketidakpastian prefiks-bersama.
type InvHorizonTrendSummary struct {
	ParamIndex             int
	ParamName              string
	Horizons               []float64
	InvHorizons            []float64
	Biases                 []float64
	BiasMCSEs              []float64
	RMSEs                  []float64
	BiasTimesH             []float64
	BiasTimesHMCSE         []float64
	Trend                  stats.LinearTrend
	OLSInterceptPrefixMCSE float64
	OLSSlopePrefixMCSE     float64
	WLSInterceptPrefixMCSE float64
	WLSSlopePrefixMCSE     float64
}

// AnalyzeLongSpanTrend mengestimasi tren linier bias terhadap 1/H dan menguji bias * H vs 4.
// Memakai varians antar-replikasi dari kombinasi linier terobservasi untuk menangkap korelasi prefiks-bersama.
func AnalyzeLongSpanTrend(longSpanResults []StudyConfigResult, paramIndex int) InvHorizonTrendSummary {
	nH := len(longSpanResults)
	horizons := make([]float64, nH)
	invH := make([]float64, nH)
	biases := make([]float64, nH)
	mcses := make([]float64, nH)
	rmses := make([]float64, nH)
	bTimesH := make([]float64, nH)
	bTimesHSE := make([]float64, nH)

	for h, res := range longSpanResults {
		H := res.Spec.HorizonYears
		p := res.Params[paramIndex]
		horizons[h] = H
		invH[h] = 1.0 / H
		biases[h] = p.Bias
		mcses[h] = p.BiasMCSE
		rmses[h] = p.RMSE
		bTimesH[h] = p.Bias * H
		bTimesHSE[h] = p.BiasMCSE * H
	}

	trend := stats.FitLinearTrend(invH, biases, mcses)

	// Hitung MCSE tepat dengan prefiks-bersama: evaluasi slope_r dan intercept_r per replikasi r.
	var olsIntMCSE, olsSlopeMCSE, wlsIntMCSE, wlsSlopeMCSE float64
	if nH >= 2 {
		nRep := len(longSpanResults[0].Replications)
		canPrefix := nRep > 1
		for h := 1; h < nH; h++ {
			if len(longSpanResults[h].Replications) != nRep {
				canPrefix = false
				break
			}
		}
		if canPrefix {
			mx := stats.Mean(invH)
			olsIntW := make([]float64, nH)
			for h := 0; h < nH; h++ {
				olsIntW[h] = (1.0 / float64(nH)) - mx*trend.OLSWeightsX[h]
			}
			// Bobot intercept WLS
			wlsIntW := make([]float64, nH)
			if trend.WLSWeightsAvailable {
				tau := make([]float64, nH)
				sumTau, sumTauX := 0.0, 0.0
				for h := 0; h < nH; h++ {
					tau[h] = 1.0 / (mcses[h] * mcses[h])
					sumTau += tau[h]
					sumTauX += tau[h] * invH[h]
				}
				xBarW := sumTauX / sumTau
				for h := 0; h < nH; h++ {
					wlsIntW[h] = (tau[h] / sumTau) - xBarW*trend.WLSWeightsX[h]
				}
			} else {
				copy(wlsIntW, olsIntW)
			}

			olsSlopes := make([]float64, nRep)
			olsInts := make([]float64, nRep)
			wlsSlopes := make([]float64, nRep)
			wlsInts := make([]float64, nRep)
			trueVal := longSpanResults[0].Params[paramIndex].TrueValue

			for r := 0; r < nRep; r++ {
				var os, oi, ws, wi float64
				for h := 0; h < nH; h++ {
					b := longSpanResults[h].Replications[r].EstParams[paramIndex] - trueVal
					os += trend.OLSWeightsX[h] * b
					oi += olsIntW[h] * b
					ws += trend.WLSWeightsX[h] * b
					wi += wlsIntW[h] * b
				}
				olsSlopes[r] = os
				olsInts[r] = oi
				wlsSlopes[r] = ws
				wlsInts[r] = wi
			}
			olsSlopeMCSE = stats.SE(olsSlopes)
			olsIntMCSE = stats.SE(olsInts)
			wlsSlopeMCSE = stats.SE(wlsSlopes)
			wlsIntMCSE = stats.SE(wlsInts)
		}
	}

	return InvHorizonTrendSummary{
		ParamIndex:             paramIndex,
		ParamName:              ParamNames[paramIndex],
		Horizons:               horizons,
		InvHorizons:            invH,
		Biases:                 biases,
		BiasMCSEs:              mcses,
		RMSEs:                  rmses,
		BiasTimesH:             bTimesH,
		BiasTimesHMCSE:         bTimesHSE,
		Trend:                  trend,
		OLSInterceptPrefixMCSE: olsIntMCSE,
		OLSSlopePrefixMCSE:     olsSlopeMCSE,
		WLSInterceptPrefixMCSE: wlsIntMCSE,
		WLSSlopePrefixMCSE:     wlsSlopeMCSE,
	}
}

// InfillParamScaling merangkum perbandingan metrik infill bulanan -> mingguan -> harian pada H=20y.
type InfillParamScaling struct {
	Name                   string
	MonthlyBias            float64
	WeeklyBias             float64
	DailyBias              float64
	MonthlyRMSE            float64
	WeeklyRMSE             float64
	DailyRMSE              float64
	MonthlyEmpSD           float64
	WeeklyEmpSD            float64
	DailyEmpSD             float64
	DailyOverMonthlyRMSE   float64
	DailyOverMonthlyEmpSD  float64
	DailyOverMonthlyHessSE float64
	SqrtNObsRatio          float64 // sqrt(T_monthly / T_daily) sebagai pembanding laju sqrt(N)
}

// AnalyzeInfillScaling menyusun tabel perbandingan perubahan frekuensi observasi pada horizon tetap 20 tahun.
func AnalyzeInfillScaling(monthly, weekly, daily StudyConfigResult) [11]InfillParamScaling {
	var out [11]InfillParamScaling
	sqrtRatio := math.Sqrt(float64(monthly.Spec.TSteps) / float64(daily.Spec.TSteps))
	for i := 0; i < 11; i++ {
		pm, pw, pd := monthly.Params[i], weekly.Params[i], daily.Params[i]
		safeDiv := func(a, b float64) float64 {
			if b == 0 {
				return 0
			}
			return a / b
		}
		out[i] = InfillParamScaling{
			Name:                   pm.Name,
			MonthlyBias:            pm.Bias,
			WeeklyBias:             pw.Bias,
			DailyBias:              pd.Bias,
			MonthlyRMSE:            pm.RMSE,
			WeeklyRMSE:             pw.RMSE,
			DailyRMSE:              pd.RMSE,
			MonthlyEmpSD:           pm.EmpiricalSD,
			WeeklyEmpSD:            pw.EmpiricalSD,
			DailyEmpSD:             pd.EmpiricalSD,
			DailyOverMonthlyRMSE:   safeDiv(pd.RMSE, pm.RMSE),
			DailyOverMonthlyEmpSD:  safeDiv(pd.EmpiricalSD, pm.EmpiricalSD),
			DailyOverMonthlyHessSE: safeDiv(pd.MeanHessianSE, pm.MeanHessianSE),
			SqrtNObsRatio:          sqrtRatio,
		}
	}
	return out
}
