package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

type Config struct {
	masterSeed, numDays, parityDays, parityEventsPerDay, ticksPerDay  int64
	ewmaAlphaMilli, storeSBase, storeUpperS, whSBase, whUpperS, whCap int64
	carrierCapDaily, kitFee, holdingCostPerDay, backorderPenalty      int64
	disruptProbPct, shockProbPct, creditMarginBps, kChk               int64
}
type Account struct {
	cash, creditLimit              int64
	onHand, inTransit, backorder   [4]int64
	ewmaMilli                      [4]int32
	status, plantState, disruptRem int32
}
type SaleRec struct {
	day, tick       int64
	store, sku, qty int32
	rev             int64
}
type ShipRec struct {
	day, tick                           int64
	src, dst, carrier, sku, qty, status int32
	goods, freight                      int64
}
type ProdRec struct {
	day, tick     int64
	fac, sku, qty int32
	cost          int64
}
type KitRec struct {
	day, tick    int64
	wh, qty, imb int32
	fee          int64
}
type FeeRec struct {
	day, tick          int64
	fromID, toID, kind int32
	amount             int64
}
type Sim struct {
	cfg                                                                 Config
	acc                                                                 [numAccounts]Account
	rng                                                                 [numAccounts + 2]Xoshiro256
	heap                                                                []SimEvent
	nextSeq                                                             uint64
	initCash, totalKits, allocTotalReq, allocTotalFill                  int64
	allocCases, allocNonzero, invViolations                             int64
	initStock, totalProduced, totalSold                                 [4]int64
	branches                                                            [16]int64
	carrierDayVol, carrierTierCnt                                       [3]int64
	logEnabled                                                          bool
	dailyDemand, dayInvSnap                                             [6][4]int64
	dailyCostEchelon                                                    [6][4]int64
	dailyFillStore                                                      [10][2]int64
	shockCnt, shockResp                                                 [3]int64
	kitImbSeries                                                        []int64
	leadtimePairs                                                       [][2]int64
	binDemand                                                           [4][100]int64
	binOrders                                                           [4][100]int64
	binStoreDem                                                         [10][100]int64
	shockEvs                                                            [][3]int
	nTransships, demEventsTot, demBoEvents, demUnitsTot, demUnitsFilled int64
	carrierCapTot, carrierVolTot, capScale                              int64
	sales                                                               []SaleRec
	ships                                                               []ShipRec
	prods                                                               []ProdRec
	kits                                                                []KitRec
	fees                                                                []FeeRec
	fOrders, fAlloc, fDisrupt, fCredit, fRisk, fTrans, fShocks, fAccSer *bufio.Writer
	files                                                               []*os.File
	eventHash                                                           uint64
}

func loadConfig(path string) Config {
	b, _ := os.ReadFile(path)
	var c Config
	for _, ln := range strings.Split(string(b), "\n") {
		k, v, ok := strings.Cut(ln, "=")
		if !ok {
			continue
		}
		val, _ := strconv.ParseInt(strings.TrimSpace(v), 10, 64)
		switch strings.TrimSpace(k) {
		case "master_seed":
			c.masterSeed = val
		case "num_days":
			c.numDays = val
		case "parity_days":
			c.parityDays = val
		case "parity_events_per_day":
			c.parityEventsPerDay = val
		case "ticks_per_day":
			c.ticksPerDay = val
		case "ewma_alpha_milli":
			c.ewmaAlphaMilli = val
		case "store_s_base":
			c.storeSBase = val
		case "store_S_base":
			c.storeUpperS = val
		case "wh_s_base":
			c.whSBase = val
		case "wh_S_base":
			c.whUpperS = val
		case "wh_cap":
			c.whCap = val
		case "carrier_cap_daily":
			c.carrierCapDaily = val
		case "kit_fee":
			c.kitFee = val
		case "holding_cost_per_day":
			c.holdingCostPerDay = val
		case "backorder_penalty":
			c.backorderPenalty = val
		case "disrupt_prob_pct":
			c.disruptProbPct = val
		case "shock_prob_pct":
			c.shockProbPct = val
		case "credit_margin_bps":
			c.creditMarginBps = val
		case "k_chk":
			c.kChk = val
		}
	}
	return c
}
func openBW(s *Sim, dir, name, hdr string) *bufio.Writer {
	f, _ := os.Create(dir + "/" + name)
	s.files = append(s.files, f)
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, hdr)
	return bw
}
func newSim(c Config, logEn bool, outDir string) *Sim {
	return newSimScaled(c, logEn, outDir, 2000)
}
func newSimScaled(c Config, logEn bool, outDir string, evPerDay int64) *Sim {
	cs := evPerDay / 2000
	if cs < 1 {
		cs = 1
	}
	s := &Sim{cfg: c, nextSeq: 1, logEnabled: logEn, eventHash: 0xCBF29CE484222325, capScale: cs}
	for i := 0; i < numAccounts+2; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(uint32(i/10), uint32(i))
		s.rng[i] = x256Seed(splitmix64(&x))
	}
	for i := 0; i < numStores; i++ {
		s.acc[i].cash, s.acc[i].creditLimit = 2000000000*cs, 500000000*cs
		s.acc[i].onHand = [4]int64{60 * cs, 60 * cs, 60 * cs, 40 * cs}
		s.acc[i].ewmaMilli = [4]int32{12000, 12000, 12000, 12000}
	}
	for i := 10; i < 13; i++ {
		s.acc[i].cash, s.acc[i].creditLimit = 5000000000*cs, 1000000000*cs
		s.acc[i].onHand = [4]int64{1600 * cs, 1400 * cs, 1200 * cs, 600 * cs}
	}
	s.acc[13].cash, s.acc[13].creditLimit = 5000000000*cs, 1000000000*cs
	s.acc[13].onHand = [4]int64{1800 * cs, 1300 * cs, 900 * cs, 500 * cs}
	for i := 14; i < 17; i++ {
		s.acc[i].cash, s.acc[i].creditLimit = 5000000000*cs, 2000000000*cs
		s.acc[i].onHand[i-14] = 4000 * cs
	}
	for i := 17; i < 20; i++ {
		s.acc[i].cash, s.acc[i].creditLimit = 500000000*cs, 200000000*cs
	}
	s.acc[accCust].cash, s.acc[accSupp].cash, s.acc[accFee].cash = 50000000000*cs, 50000000000*cs, 1000000000*cs
	for i := 0; i < numAccounts; i++ {
		s.initCash += s.acc[i].cash
		for k := 0; k < 4; k++ {
			s.initStock[k] += s.acc[i].onHand[k]
		}
	}
	if logEn && outDir != "" {
		s.fOrders = openBW(s, outDir, "orders_log.csv", "day,tick,seq,etype,src,dst,sku,qty")
		s.fAlloc = openBW(s, outDir, "allocations.csv", "day,wh,sku,store,req,fill")
		s.fDisrupt = openBW(s, outDir, "disruptions.csv", "day,tick,fac,state,dur")
		s.fCredit = openBW(s, outDir, "credit_calls.csv", "day,acc,cash,req_margin,status")
		s.fRisk = openBW(s, outDir, "risk_actions.csv", "day,tick,acc,action,amount")
		s.fTrans = openBW(s, outDir, "transshipments.csv", "day,tick,src,dst,sku,qty")
		s.fShocks = openBW(s, outDir, "shocks_log.csv", "day,tick,cls,store,mag")
		s.fAccSer = openBW(s, outDir, "accounts_series.csv", "day,acc,cash,oh_tot,it_tot,status")
	}
	return s
}
func (s *Sim) closeLogs() {
	for _, bw := range []*bufio.Writer{s.fOrders, s.fAlloc, s.fDisrupt, s.fCredit, s.fRisk, s.fTrans, s.fShocks, s.fAccSer} {
		if bw != nil {
			bw.Flush()
		}
	}
	for _, f := range s.files {
		f.Close()
	}
}
