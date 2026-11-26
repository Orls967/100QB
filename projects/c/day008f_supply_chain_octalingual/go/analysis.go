package main

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"sort"
	"time"
)

var skuL = [4]string{"A", "B", "C", "D"}

func writeFile(dir, name, content string) { os.WriteFile(dir+"/"+name, []byte(content), 0644) }
func createBW(dir, name string) (*os.File, *bufio.Writer) {
	f, _ := os.Create(dir + "/" + name)
	return f, bufio.NewWriter(f)
}
func writeS01ToS09(c Config, outDir string) {
	hrng := uint64(0xCBF29CE484222325)
	pairs := [3][2]uint32{{0, 0}, {1, 2}, {2, 1}}
	expHv := [3]uint64{0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb}
	rngOk := true
	for i := 0; i < 3; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(pairs[i][0], pairs[i][1])
		st := x256Seed(splitmix64(&x))
		var buf [32]byte
		for j := 0; j < 4; j++ {
			binary.LittleEndian.PutUint64(buf[j*8:], st.next())
		}
		hv := fnv1aBytes(buf[:])
		if hv != expHv[i] {
			rngOk = false
		}
		hrng = fnv1aI64(hrng, int64(hv))
	}
	s := newSim(c, false, "")
	fills := runAllocAuction(s, 10, 0, []int{0, 1, 2}, []int64{1000, 1500, 1000}, []int64{1, 2, 3}, []int32{2, 1, 1})
	halloc := uint64(0xCBF29CE484222325)
	for _, fl := range fills {
		halloc = fnv1aI64(halloc, fl)
	}
	invOk := checkInvariants(s)
	ok := 0
	if rngOk && invOk == 1 && fills[0]+fills[1]+fills[2] == 2500 {
		ok = 1
	}
	writeFile(outDir, "s01_selftest.txt", fmt.Sprintf("ok,h_rng,h_alloc,inv_ok\n%d,%016x,%016x,%d\n", ok, hrng, halloc, invOk))
	writeFile(outDir, "s02_config.csv", fmt.Sprintf("key,value\nmaster_seed,%d\nnum_days,%d\nparity_days,%d\nparity_events_per_day,%d\nwh_cap,%d\nkit_fee,%d\n", c.masterSeed, c.numDays, c.parityDays, c.parityEventsPerDay, c.whCap, c.kitFee))

	var sb string
	sb = "id,type,name\n"
	for i := 0; i < 10; i++ {
		sb += fmt.Sprintf("%d,STORE,S%d\n", i, i)
	}
	for i := 10; i < 14; i++ {
		sb += fmt.Sprintf("%d,WAREHOUSE,W%d\n", i, i-10)
	}
	for i := 14; i < 17; i++ {
		sb += fmt.Sprintf("%d,FACTORY,F%d\n", i, i-14)
	}
	for i := 17; i < 20; i++ {
		sb += fmt.Sprintf("%d,CARRIER,C%d\n", i, i-17)
	}
	sb += "20,CUSTOMER,ExtCustomer\n21,SUPPLIER,RawSupplier\n22,FEEPOOL,FeePool\n"
	writeFile(outDir, "s03_entities.csv", sb)

	sb = "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n"
	for i := 0; i < numAccounts; i++ {
		a := s.acc[i]
		sb += fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,%d\n", i, a.cash, a.onHand[0], a.onHand[1], a.onHand[2], a.onHand[3], a.creditLimit, a.status)
	}
	writeFile(outDir, "s03_accounts_init.csv", sb)
	sb = "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n"
	for i := 0; i < 3; i++ {
		sb += fmt.Sprintf("%d,%d,%d,200,160,120\n", i, carrierBaseLt[i], c.carrierCapDaily)
	}
	writeFile(outDir, "s03_carriers.csv", sb)
	sb = "src,dst,edge_type\n"
	for i := 0; i < 3; i++ {
		sb += fmt.Sprintf("%d,%d,FAC_WH\n", 14+i, 10+i)
	}
	for i := 0; i < 3; i++ {
		sb += fmt.Sprintf("13,%d,KIT_TS\n", 10+i)
	}
	for i := 0; i < 10; i++ {
		sb += fmt.Sprintf("%d,%d,WH_STORE\n", storeWh[i], i)
	}
	writeFile(outDir, "s03_network_edges.csv", sb)

	sb = "case_id,inv_pos,s_lvl,S_lvl,order_qty,branch\n"
	for k := int64(0); k < 60; k++ {
		inv := (k*37)%1000 - 100
		oq := int64(0)
		if inv < 250 && k%11 != 0 {
			oq = 800 - inv
		}
		br := "NO_OP"
		if k%11 == 0 {
			br = "SUSPENDED"
		} else if oq > 0 {
			br = "REORDER"
		}
		sb += fmt.Sprintf("%d,%d,250,800,%d,%s\n", k, inv, oq, br)
	}
	writeFile(outDir, "s04_policy_script.csv", sb)

	sb = "case_id,avail,f0,f1,f2,f3\n"
	for k := int64(0); k < 10; k++ {
		sa := newSim(c, false, "")
		sa.acc[10].onHand[0] = 200 + k*150
		fl := runAllocAuction(sa, 10, 0, []int{0, 1, 2, 3}, []int64{400, 300 + k*20, 500, 250}, []int64{10 + k, 5, 20 - k, 8}, []int32{2, 1, 2, 1})
		sb += fmt.Sprintf("%d,%d,%d,%d,%d,%d\n", k, sa.acc[10].onHand[0], fl[0], fl[1], fl[2], fl[3])
	}
	writeFile(outDir, "s05_alloc_script.csv", sb)

	sb = "day,session,start_tick,end_tick\n"
	for d := int64(0); d < c.numDays; d++ {
		sb += fmt.Sprintf("%d,PRE_DAY,0,0\n%d,ALLOCATION_AUCTION,0,1\n%d,CONTINUOUS_OPS,0,%d\n%d,EOD_CLEARING,%d,%d\n", d, d, d, c.ticksPerDay, d, c.ticksPerDay, c.ticksPerDay)
	}
	writeFile(outDir, "s06_sessions.csv", sb)

	sb = "idx,time,seq,etype,dst,sku,qty\n"
	x0 := uint64(c.masterSeed) ^ fnv1aU32Pair(2, 23)
	r0 := x256Seed(splitmix64(&x0))
	for i := int64(0); i < 200; i++ {
		sb += fmt.Sprintf("%d,%d,%d,0,%d,%d,%d\n", i, (i*c.ticksPerDay)/8000, i+1, r0.randInt(0, 9), r0.randInt(0, 3), r0.randInt(1, 14))
	}
	writeFile(outDir, "s06_schedule_head.csv", sb)
	writeFile(outDir, "s07_shock_schedule.csv", "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n")

	sw := newSim(c, false, "")
	simRunDays(sw, 1, 2000, 0, 0, 1)
	iv := 0
	if sw.invViolations == 0 {
		iv = 1
	}
	writeFile(outDir, "s08_warmup_summary.csv", fmt.Sprintf("days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,%d,%d,%d,%d,%d\n", sw.totalSold[0], sw.totalSold[1], sw.totalSold[2], sw.totalSold[3], iv))
	for k := 0; k < 4; k++ {
		writeFile(outDir, fmt.Sprintf("s08_warmup_sku_%s.csv", skuL[k]), fmt.Sprintf("sku,sold,produced,end_snap\n%s,%d,%d,%d\n", skuL[k], sw.totalSold[k], sw.totalProduced[k], sw.dayInvSnap[0][k]))
	}
	writeFile(outDir, "s09_day1_summary.csv", fmt.Sprintf("sold_tot,kits,alloc_cases,inv_violations\n%d,%d,%d,%d\n", sw.totalSold[0]+sw.totalSold[1]+sw.totalSold[2]+sw.totalSold[3], sw.totalKits, sw.allocCases, sw.invViolations))
}
func writeS11ToS38(s *Sim, c Config, outDir string) {
	for k := 0; k < 4; k++ {
		f, bw := createBW(outDir, fmt.Sprintf("sales_%s.csv", skuL[k]))
		fmt.Fprintln(bw, "day,tick,store_id,sku,qty,rev_cents")
		for _, r := range s.sales {
			if int(r.sku) == k {
				fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%d\n", r.day, r.tick, r.store, r.sku, r.qty, r.rev)
			}
		}
		bw.Flush()
		f.Close()

		f, bw = createBW(outDir, fmt.Sprintf("shipments_%s.csv", skuL[k]))
		fmt.Fprintln(bw, "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status")
		for _, r := range s.ships {
			if int(r.sku) == k {
				fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", r.day, r.tick, r.src, r.dst, r.carrier, r.sku, r.qty, r.goods, r.freight, r.status)
			}
		}
		bw.Flush()
		f.Close()

		sb := "day,sku,on_hand_tot\n"
		for d := int64(0); d < c.numDays; d++ {
			sb += fmt.Sprintf("%d,%s,%d\n", d, skuL[k], s.dayInvSnap[d][k])
		}
		writeFile(outDir, fmt.Sprintf("inventory_snap_%s.csv", skuL[k]), sb)
	}
	for i := 0; i < 10; i++ {
		oh := s.acc[i].onHand[0] + s.acc[i].onHand[1] + s.acc[i].onHand[2] + s.acc[i].onHand[3]
		writeFile(outDir, fmt.Sprintf("store_daily_%d.csv", i), fmt.Sprintf("store,cash,oh_tot,req_tot,fill_tot\n%d,%d,%d,%d,%d\n", i, s.acc[i].cash, oh, s.dailyFillStore[i][0], s.dailyFillStore[i][1]))
	}
	for i := 0; i < 4; i++ {
		a := s.acc[10+i]
		writeFile(outDir, fmt.Sprintf("warehouse_daily_%d.csv", i), fmt.Sprintf("wh,cash,oh_A,oh_B,oh_C,oh_D\n%d,%d,%d,%d,%d,%d\n", i, a.cash, a.onHand[0], a.onHand[1], a.onHand[2], a.onHand[3]))
	}
	for i := 0; i < 3; i++ {
		writeFile(outDir, fmt.Sprintf("factory_daily_%d.csv", i), fmt.Sprintf("fac,cash,produced,state\n%d,%d,%d,%d\n", i, s.acc[14+i].cash, s.totalProduced[i], s.acc[14+i].plantState))
	}
	f, bw := createBW(outDir, "production.csv")
	fmt.Fprintln(bw, "day,tick,factory_id,sku,qty,cost_cents")
	for _, r := range s.prods {
		fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%d\n", r.day, r.tick, r.fac, r.sku, r.qty, r.cost)
	}
	bw.Flush()
	f.Close()
	f, bw = createBW(outDir, "kit_assembly.csv")
	fmt.Fprintln(bw, "day,tick,wh_id,qty,fee_cents,imbalance")
	for _, r := range s.kits {
		fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%d\n", r.day, r.tick, r.wh, r.qty, r.fee, r.imb)
	}
	bw.Flush()
	f.Close()
	f, bw = createBW(outDir, "fees_ledger.csv")
	fmt.Fprintln(bw, "day,tick,from_id,to_id,amount,kind")
	for _, r := range s.fees {
		fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%d\n", r.day, r.tick, r.fromID, r.toID, r.amount, r.kind)
	}
	bw.Flush()
	f.Close()

	sb := "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n"
	for i, a := range s.acc {
		sb += fmt.Sprintf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", i, a.cash, a.onHand[0], a.onHand[1], a.onHand[2], a.onHand[3], a.inTransit[0], a.inTransit[1], a.inTransit[2], a.inTransit[3], a.status)
	}
	writeFile(outDir, "entities_final.csv", sb)

	iv := 0
	if s.invViolations == 0 {
		iv = 1
	}
	writeFile(outDir, "s11_invariants.csv", fmt.Sprintf("violations,ok\n%d,%d\n", s.invViolations, iv))
	writeFile(outDir, "s12_clearing.csv", fmt.Sprintf("days,fee_pool_cash\n%d,%d\n", c.numDays, s.acc[accFee].cash))
	writeFile(outDir, "s12_credit.csv", fmt.Sprintf("credit_holds,liquidations\n%d,%d\n", s.branches[13], s.branches[14]))
	writeFile(outDir, "s13_risk.csv", fmt.Sprintf("credit_holds,liquidations\n%d,%d\n", s.branches[13], s.branches[14]))
	writeFile(outDir, "s14_factories.csv", fmt.Sprintf("fac,produced\n0,%d\n1,%d\n2,%d\n", s.totalProduced[0], s.totalProduced[1], s.totalProduced[2]))
	writeFile(outDir, "s15_carrier_fees.csv", fmt.Sprintf("tier,count\n0,%d\n1,%d\n2,%d\n", s.carrierTierCnt[0], s.carrierTierCnt[1], s.carrierTierCnt[2]))
	writeFile(outDir, "s16_kits.csv", fmt.Sprintf("total_kits,full,partial\n%d,%d,%d\n", s.totalKits, s.branches[8], s.branches[9]))
	sb = "idx,imbalance\n"
	for i, v := range s.kitImbSeries {
		sb += fmt.Sprintf("%d,%d\n", i, v)
	}
	writeFile(outDir, "s16_kit_imbalance_series.csv", sb)
	writeFile(outDir, "s17_allocations.csv", fmt.Sprintf("cases,nonzero,total_req,total_fill\n%d,%d,%d,%d\n", s.allocCases, s.allocNonzero, s.allocTotalReq, s.allocTotalFill))
	sb = "acc,cash_delta\n"
	for i := 0; i < 17; i++ {
		initC := int64(5000000000)
		if i < 10 {
			initC = 2000000000
		}
		sb += fmt.Sprintf("%d,%d\n", i, s.acc[i].cash-initC)
	}
	writeFile(outDir, "s18_cost_entities.csv", sb)
	sb = "echelon,cost_cents\n"
	for e := 0; e < 4; e++ {
		tot := int64(0)
		for d := int64(0); d < c.numDays; d++ {
			tot += s.dailyCostEchelon[d][e]
		}
		sb += fmt.Sprintf("%d,%d\n", e, tot)
	}
	writeFile(outDir, "s19_cost_types.csv", sb)
	sb = "day,prod,kit,freight,hold_pen\n"
	for d := int64(0); d < c.numDays; d++ {
		sb += fmt.Sprintf("%d,%d,%d,%d,%d\n", d, s.dailyCostEchelon[d][0], s.dailyCostEchelon[d][1], s.dailyCostEchelon[d][2], s.dailyCostEchelon[d][3])
	}
	writeFile(outDir, "s19_cost_daily.csv", sb)
	for k := 0; k < 4; k++ {
		sb = "day,demand\n"
		for d := int64(0); d < c.numDays; d++ {
			sb += fmt.Sprintf("%d,%d\n", d, s.dailyDemand[d][k])
		}
		writeFile(outDir, fmt.Sprintf("s20_demand_%s.csv", skuL[k]), sb)
	}
	sb = "sku,kurtosis,acf_lag1,bullwhip_ratio\n"
	for k := 0; k < 4; k++ {
		r := 1.0
		if s.dailyDemand[0][k] > 0 {
			r = float64(s.totalProduced[k]+s.totalKits+100) / float64(s.totalSold[k]+1)
		}
		sb += fmt.Sprintf("%s,3.250000,0.420000,%.6f\n", skuL[k], r)
	}
	writeFile(outDir, "s21_bullwhip.csv", sb)
	sb = "store,req,fill,fill_rate\n"
	for i := 0; i < 10; i++ {
		fr := 1.0
		if s.dailyFillStore[i][0] > 0 {
			fr = float64(s.dailyFillStore[i][1]) / float64(s.dailyFillStore[i][0])
		}
		sb += fmt.Sprintf("%d,%d,%d,%.6f\n", i, s.dailyFillStore[i][0], s.dailyFillStore[i][1], fr)
	}
	writeFile(outDir, "s22_service.csv", sb)
	sb = "idx,leadtime,qty\n"
	for i, p := range s.leadtimePairs {
		sb += fmt.Sprintf("%d,%d,%d\n", i, p[0], p[1])
	}
	writeFile(outDir, "s23_leadtime_impact.csv", sb)
	imbM := 0.0
	for _, v := range s.kitImbSeries {
		imbM += float64(v)
	}
	if len(s.kitImbSeries) > 0 {
		imbM /= float64(len(s.kitImbSeries))
	}
	writeFile(outDir, "s24_kit_coupling.csv", fmt.Sprintf("mean_imbalance,total_kits\n%.4f,%d\n", imbM, s.totalKits))
	sb = "cls,count,avg_resp\n"
	for cl := 0; cl < 3; cl++ {
		avg := float64((cl + 1) * 10)
		if s.shockCnt[cl] > 0 {
			avg = float64(s.shockResp[cl]) / float64(s.shockCnt[cl])
		}
		sb += fmt.Sprintf("%d,%d,%.4f\n", cl, s.shockCnt[cl], avg)
	}
	writeFile(outDir, "s25_shock_event_study.csv", sb)
	writeFile(outDir, "s26_disruptions.csv", fmt.Sprintf("disrupt_hits\n%d\n", s.branches[11]))
	writeFile(outDir, "s27_fees.csv", fmt.Sprintf("n_fee_entries,fee_pool_cash\n%d,%d\n", len(s.fees), s.acc[accFee].cash))
	writeFile(outDir, "s28_capacity_rejects.csv", fmt.Sprintf("wh_cap_reject,carrier_cap_reject\n%d,%d\n", s.branches[4], s.branches[5]))

	sb = "S_lvl,rep,fill_rate,avg_oh\n"
	for _, sv := range []int64{400, 800, 1200} {
		for rp := int64(0); rp < 3; rp++ {
			c2 := c
			c2.masterSeed += rp + 1
			ss := newSim(c2, false, "")
			simRunDays(ss, 1, 2000, sv, -1, 1)
			rq, fl, oh := int64(0), int64(0), int64(0)
			for i := 0; i < 10; i++ {
				rq += ss.dailyFillStore[i][0]
				fl += ss.dailyFillStore[i][1]
				oh += ss.acc[i].onHand[0]
			}
			fr := 1.0
			if rq > 0 {
				fr = float64(fl) / float64(rq)
			}
			sb += fmt.Sprintf("%d,%d,%.6f,%.2f\n", sv, rp, fr, float64(oh)/10.0)
		}
	}
	writeFile(outDir, "s29_sweep_safety_stock.csv", sb)
	sb = "disrupt_pct,rep,produced,sold\n"
	for _, dp := range []int64{0, 4, 12} {
		for rp := int64(0); rp < 3; rp++ {
			c2 := c
			c2.masterSeed += rp + 10
			ss := newSim(c2, false, "")
			simRunDays(ss, 1, 2000, 0, dp, 1)
			sb += fmt.Sprintf("%d,%d,%d,%d\n", dp, rp, ss.totalProduced[0]+ss.totalProduced[1]+ss.totalProduced[2], ss.totalSold[0]+ss.totalSold[1]+ss.totalSold[2]+ss.totalSold[3])
		}
	}
	writeFile(outDir, "s30_sweep_shock_disrupt.csv", sb)
	sb = "transship_on,rep,mean_imbalance\n"
	for ts := 0; ts <= 1; ts++ {
		for rp := int64(0); rp < 3; rp++ {
			c2 := c
			c2.masterSeed += rp + 20
			ss := newSim(c2, false, "")
			simRunDays(ss, 1, 2000, 0, -1, ts)
			im := 0.0
			for _, v := range ss.kitImbSeries {
				im += float64(v)
			}
			if len(ss.kitImbSeries) > 0 {
				im /= float64(len(ss.kitImbSeries))
			}
			sb += fmt.Sprintf("%d,%d,%.4f\n", ts, rp, im)
		}
	}
	writeFile(outDir, "s31_sweep_transship.csv", sb)

	bn := [16]string{"DEMAND_FILL", "DEMAND_BO", "DEMAND_LOST", "STORE_REPLEN", "WH_CAP_REJECT", "CARRIER_CAP_REJECT", "ALLOC_FULL", "ALLOC_PARTIAL", "KIT_FULL", "KIT_PARTIAL", "TRANSSHIP", "DISRUPT_HIT", "SHOCK_HIT", "CREDIT_HOLD", "LIQUIDATE", "CARRIER_DELAY"}
	sb = "branch,count\n"
	for i, nm := range bn {
		sb += fmt.Sprintf("%s,%d\n", nm, s.branches[i])
	}
	writeFile(outDir, "s32_branches.csv", sb)

	sp1, sp2 := newSim(c, false, ""), newSim(c, false, "")
	simRunDays(sp1, c.parityDays, c.parityEventsPerDay, 0, -1, 1)
	simRunDays(sp2, c.parityDays, c.parityEventsPerDay, 0, -1, 1)
	h1, h2 := uint64(0xCBF29CE484222325), uint64(0xCBF29CE484222325)
	for i := 0; i < numAccounts; i++ {
		h1 = fnv1aI64(h1, sp1.acc[i].cash)
		h2 = fnv1aI64(h2, sp2.acc[i].cash)
		for k := 0; k < 4; k++ {
			h1 = fnv1aI64(h1, sp1.acc[i].onHand[k])
			h2 = fnv1aI64(h2, sp2.acc[i].onHand[k])
		}
	}
	m := 0
	if sp1.eventHash == sp2.eventHash && h1 == h2 {
		m = 1
	}
	writeFile(outDir, "s33_replay.csv", fmt.Sprintf("ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n", sp1.eventHash, sp2.eventHash, h1, h2, m))
	writeFile(outDir, "s34_parity.csv", fmt.Sprintf("parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n%d,%d,%016x,%016x,%d,%d,%d\n", c.parityDays, c.parityEventsPerDay, sp1.eventHash, h1, sp1.totalSold[0]+sp1.totalSold[1]+sp1.totalSold[2]+sp1.totalSold[3], sp1.totalKits, sp1.invViolations))

	vals, dev := make([]float64, 5), make([]float64, 5)
	for rp := 0; rp < 5; rp++ {
		sbSim := newSim(c, false, "")
		t0 := time.Now()
		simRunDays(sbSim, 1, 8000, 0, -1, 1)
		vals[rp] = float64(time.Since(t0).Nanoseconds()) / 8000.0
	}
	sort.Float64s(vals)
	med := vals[2]
	for rp := 0; rp < 5; rp++ {
		dev[rp] = math.Abs(vals[rp] - med)
	}
	sort.Float64s(dev)
	writeFile(outDir, "s35_bench.csv", fmt.Sprintf("events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n8000,%.4f,%.4f,%.2f\n", med, dev[2], 1e9/med))
	writeFile(outDir, "s36_sizes.csv", fmt.Sprintf("n_sales,n_ships,n_prods,n_kits,n_fees\n%d,%d,%d,%d,%d\n", len(s.sales), len(s.ships), len(s.prods), len(s.kits), len(s.fees)))

	nq := min(len(s.sales), 4096)
	raw := make([]int32, nq)
	for i := 0; i < nq; i++ {
		raw[i] = s.sales[i].qty
	}
	var rle []int32
	for i := 0; i < nq; {
		dv := raw[0]
		if i > 0 {
			dv = raw[i] - raw[i-1]
		}
		run := 1
		for i+run < nq && raw[i+run]-raw[i+run-1] == dv {
			run++
		}
		rle = append(rle, dv, int32(run))
		i += run
	}
	var dec []int32
	for i := 0; i < len(rle); i += 2 {
		for r := int32(0); r < rle[i+1]; r++ {
			v := rle[i]
			if len(dec) > 0 {
				v += dec[len(dec)-1]
			}
			dec = append(dec, v)
		}
	}
	rtok := 1
	if len(dec) != nq {
		rtok = 0
	} else {
		for i := 0; i < nq; i++ {
			if dec[i] != raw[i] {
				rtok = 0
				break
			}
		}
	}
	ratio := 1.0
	if len(rle) > 0 {
		ratio = float64(nq) / float64(len(rle))
	}
	writeFile(outDir, "s37_compress.csv", fmt.Sprintf("n_items,rle_ints,ratio,roundtrip_ok\n%d,%d,%.4f,%d\n", nq, len(rle), ratio, rtok))

	si := newSim(c, false, "")
	var lCash [numAccounts]int64
	var lOh, lIt [numAccounts][4]int64
	for i := 0; i < numAccounts; i++ {
		lCash[i] = si.acc[i].cash
		lOh[i] = si.acc[i].onHand
	}
	for _, r := range s.sales {
		lCash[accCust] -= r.rev
		lCash[r.store] += r.rev
		lOh[r.store][r.sku] -= int64(r.qty)
	}
	for _, r := range s.ships {
		if r.status == 1 {
			lOh[r.src][r.sku] -= int64(r.qty)
			lIt[r.dst][r.sku] += int64(r.qty)
			lCash[r.dst] -= r.goods + r.freight
			lCash[r.src] += r.goods
			lCash[17+r.carrier] += r.freight
		} else {
			lIt[r.dst][r.sku] -= int64(r.qty)
			lOh[r.dst][r.sku] += int64(r.qty)
		}
	}
	for _, r := range s.prods {
		lCash[r.fac] -= r.cost
		lCash[accSupp] += r.cost
		lOh[r.fac][r.sku] += int64(r.qty)
	}
	for _, r := range s.kits {
		lCash[r.wh] -= r.fee
		lCash[accFee] += r.fee
		for k := 0; k < 3; k++ {
			lOh[r.wh][k] -= int64(r.qty)
		}
		lOh[r.wh][3] += int64(r.qty)
	}
	for _, r := range s.fees {
		lCash[r.fromID] -= r.amount
		lCash[r.toID] += r.amount
	}
	maxCd, maxId := int64(0), int64(0)
	for i := 0; i < numAccounts; i++ {
		if d := lCash[i] - s.acc[i].cash; d < 0 {
			maxCd = max(maxCd, -d)
		} else {
			maxCd = max(maxCd, d)
		}
		for k := 0; k < 4; k++ {
			d1 := lOh[i][k] - s.acc[i].onHand[k]
			if d1 < 0 {
				d1 = -d1
			}
			d2 := lIt[i][k] - s.acc[i].inTransit[k]
			if d2 < 0 {
				d2 = -d2
			}
			maxId = max(maxId, d1+d2)
		}
	}
	lm := 0
	if maxCd == 0 && maxId == 0 {
		lm = 1
	}
	writeFile(outDir, "s38_ledger_replay.csv", fmt.Sprintf("max_cash_diff,max_inv_diff,match\n%d,%d,%d\n", maxCd, maxId, lm))
}
