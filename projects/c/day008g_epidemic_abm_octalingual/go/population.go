package main

import (
	"os"
	"strconv"
	"strings"
)

const (
	BP      int64 = 10000
	MaxDays       = 200
	StS     int32 = 0
	StE     int32 = 1
	StIA    int32 = 2
	StIS    int32 = 3
	StH     int32 = 4
	StICU   int32 = 5
	StR     int32 = 6
	StD     int32 = 7
)

type Config struct {
	masterSeed, parityN, parityDays, mainDays           int64
	betaBpInit, asympProbBp, asympRelBp                 int64
	latentDays, infDaysMin, infDaysMax                  int64
	hospBedsPer10k, icuBedsPer10k                       int64
	waningMinDays, waningBp, variantDay, variantMultBp  int64
	schoolCloseDay, lockdownDay, maskDay, maskEffBp     int64
	testDay, testProbBp, traceProbBp                    int64
	vaxDay, vaxDailyBp, vaxEff1Bp, vaxEff2Bp, vaxHospBp int64
	initInfected                                        int64
}

type Person struct {
	ageGrp, hhID, swID, status, variant         int32
	stateDay, durDays, infDay, vaxDoses, vaxDay int32
	compliance, quarUntil, superFlag, offspring int32
	everInf, hospRej                            int32
}

type RunParams struct {
	n, days, betaBp                              int64
	complianceOvr, testOvr, vaxDayOvr, vaxCovOvr int64
	schoolOvr, lockOvr, varDayOvr, maskOvr       int64
	contactScaleBp, vaxPerfect                   int64
	intervOn                                     int32
}

func defaultRunParams(n, days, betaBp int64, intervOn int32) RunParams {
	return RunParams{
		n: n, days: days, betaBp: betaBp,
		complianceOvr: -1, testOvr: -1, vaxDayOvr: -1, vaxCovOvr: -1,
		schoolOvr: -1, lockOvr: -1, varDayOvr: -1, maskOvr: -1,
		contactScaleBp: 10000, vaxPerfect: 0, intervOn: intervOn,
	}
}

type TransRec struct {
	day    int32
	src    int
	dst    int
	lay    int
	vari   int32
	genInt int32
}

type ProgRec struct {
	day    int32
	pid    int
	fromSt int32
	toSt   int32
}

type Sim struct {
	cfg            Config
	rp             RunParams
	n              int64
	hospCap        int64
	icuCap         int64
	hospOcc        int64
	icuOcc         int64
	pop            []Person
	nHH            int
	nSW            int
	hhOff          []int
	hhMem          []int
	swOff          []int
	swMem          []int
	seenDay        []int32
	rng            [5]Xoshiro256
	candTgt        []int
	candSrc        []int
	candLay        []int
	candVar        []int32
	dailyCnt       [MaxDays][8]int64
	dailyIncAge    [MaxDays][3]int64
	dailyDeathsAge [MaxDays][3]int64
	dailyHosp      [MaxDays][4]int64
	dailyTests     [MaxDays][2]int64
	dailyQuar      [MaxDays]int64
	dailyVax       [MaxDays][2]int64
	dailyVar       [MaxDays][2]int64
	layerInf       [3]int64
	ageInf         [3]int64
	ageTot         [3]int64
	branches       [21]int64
	invViolations  int64
	totalTrans     int64
	totalReinf     int64
	eventHash      uint64
	logEnabled     bool
	trans          []TransRec
	progs          []ProgRec
}

func loadConfig(path string) Config {
	data, err := os.ReadFile(path)
	if err != nil {
		panic(err)
	}
	var c Config
	for _, line := range strings.Split(string(data), "\n") {
		k, v, ok := strings.Cut(line, "=")
		if !ok {
			continue
		}
		val, _ := strconv.ParseInt(strings.TrimSpace(v), 10, 64)
		switch strings.TrimSpace(k) {
		case "master_seed":
			c.masterSeed = val
		case "parity_N":
			c.parityN = val
		case "parity_days":
			c.parityDays = val
		case "main_days":
			c.mainDays = val
		case "beta_bp_init":
			c.betaBpInit = val
		case "asymp_prob_bp":
			c.asympProbBp = val
		case "asymp_rel_bp":
			c.asympRelBp = val
		case "latent_days":
			c.latentDays = val
		case "inf_days_min":
			c.infDaysMin = val
		case "inf_days_max":
			c.infDaysMax = val
		case "hosp_beds_per_10k":
			c.hospBedsPer10k = val
		case "icu_beds_per_10k":
			c.icuBedsPer10k = val
		case "waning_min_days":
			c.waningMinDays = val
		case "waning_bp":
			c.waningBp = val
		case "variant_day":
			c.variantDay = val
		case "variant_mult_bp":
			c.variantMultBp = val
		case "school_close_day":
			c.schoolCloseDay = val
		case "lockdown_day":
			c.lockdownDay = val
		case "mask_day":
			c.maskDay = val
		case "mask_eff_bp":
			c.maskEffBp = val
		case "test_day":
			c.testDay = val
		case "test_prob_bp":
			c.testProbBp = val
		case "trace_prob_bp":
			c.traceProbBp = val
		case "vax_day":
			c.vaxDay = val
		case "vax_daily_bp":
			c.vaxDailyBp = val
		case "vax_eff1_bp":
			c.vaxEff1Bp = val
		case "vax_eff2_bp":
			c.vaxEff2Bp = val
		case "vax_hosp_bp":
			c.vaxHospBp = val
		case "init_infected":
			c.initInfected = val
		}
	}
	return c
}

func newSim(c Config, rp RunParams, logEn bool) *Sim {
	s := &Sim{
		cfg:        c,
		rp:         rp,
		n:          rp.n,
		hospCap:    max64((rp.n*c.hospBedsPer10k)/BP, 1),
		icuCap:     max64((rp.n*c.icuBedsPer10k)/BP, 1),
		eventHash:  0xcbf29ce484222325,
		logEnabled: logEn,
	}
	for i := uint32(0); i < 5; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(100, i)
		s.rng[i] = newXoshiro(splitmix64(&x))
	}
	n := int(rp.n)
	s.pop = make([]Person, n)
	s.seenDay = make([]int32, n)
	for i := 0; i < n; i++ {
		s.pop[i].swID = -1
		s.pop[i].infDay = -1
		s.pop[i].vaxDay = -1
		s.pop[i].quarUntil = -1
		s.seenDay[i] = -1
	}

	hhCnt := make([]int, 0, n/3+1)
	pid, hh := 0, 0
	for pid < n {
		sz := int(s.rng[0].randInt(1, 6))
		if sz > n-pid {
			sz = n - pid
		}
		hhCnt = append(hhCnt, sz)
		for k := 0; k < sz; k++ {
			s.pop[pid].hhID = int32(hh)
			u := s.rng[0].randInt(1, 100)
			if u <= 22 {
				s.pop[pid].ageGrp = 0
			} else if u <= 82 {
				s.pop[pid].ageGrp = 1
			} else {
				s.pop[pid].ageGrp = 2
			}
			s.ageTot[s.pop[pid].ageGrp]++
			s.pop[pid].compliance = int32(s.rng[0].randInt(4000, 9500))
			if s.rng[0].randInt(1, 10) == 1 {
				s.pop[pid].superFlag = 1
			}
			pid++
		}
		hh++
	}
	s.nHH = hh
	s.hhOff = make([]int, hh+1)
	for i := 0; i < hh; i++ {
		s.hhOff[i+1] = s.hhOff[i] + hhCnt[i]
	}
	s.hhMem = make([]int, n)
	cur := make([]int, hh)
	for i := 0; i < n; i++ {
		h := int(s.pop[i].hhID)
		s.hhMem[s.hhOff[h]+cur[h]] = i
		cur[h]++
	}

	nSchools := int(s.ageTot[0]/25) + 1
	nWorks := int(s.ageTot[1]/15) + 1
	s.nSW = nSchools + nWorks
	swCnt := make([]int, s.nSW)
	for i := 0; i < n; i++ {
		if s.pop[i].ageGrp == 0 {
			sw := int(s.rng[0].randInt(0, int64(nSchools)-1))
			s.pop[i].swID = int32(sw)
			swCnt[sw]++
		} else if s.pop[i].ageGrp == 1 {
			sw := nSchools + int(s.rng[0].randInt(0, int64(nWorks)-1))
			s.pop[i].swID = int32(sw)
			swCnt[sw]++
		}
	}
	s.swOff = make([]int, s.nSW+1)
	for i := 0; i < s.nSW; i++ {
		s.swOff[i+1] = s.swOff[i] + swCnt[i]
	}
	s.swMem = make([]int, n)
	swCur := make([]int, s.nSW)
	for i := 0; i < n; i++ {
		if s.pop[i].swID >= 0 {
			sw := int(s.pop[i].swID)
			s.swMem[s.swOff[sw]+swCur[sw]] = i
			swCur[sw]++
		}
	}
	return s
}

func max64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}
