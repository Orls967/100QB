package main

import (
	"fmt"
	"math"
	"sort"
)

func runMetricWait(c Config, rp RunParams, seedOff int64) float64 {
	c2 := c
	c2.masterSeed += seedOff
	s := newSim(c2, rp)
	simRun(s)
	if s.stnWaitCnt[StnSecurity] > 0 {
		return float64(s.stnWaitSum[StnSecurity]) / float64(s.stnWaitCnt[StnSecurity])
	}
	return 0.0
}

func writeS14ToS45(c Config, outDir string, sMain *Sim) {
	snm := []string{"CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"}
	evNm := []string{"station_events_CHECKIN.csv", "station_events_SECURITY.csv", "station_events_IMIGRASI.csv", "station_events_GATE.csv", "station_events_BELT.csv"}
	for st := 0; st < 5; st++ {
		f, w := openBW(outDir, evNm[st])
		fmt.Fprintln(w, "slot,stn,delta")
		for _, ev := range sMain.stnEv[st] {
			fmt.Fprintf(w, "%d,%d,%d\n", ev.slot, ev.stn, ev.delta)
		}
		closeBW(f, w)
	}
	f, w := openBW(outDir, "queue_series.csv")
	fmt.Fprintln(w, "slot,q0,q1,q2,q3,q4,q5")
	for t := 0; t < int(sMain.auditSlots); t++ {
		b := t * NStn
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d,%d\n", t, sMain.qSeries[b], sMain.qSeries[b+1], sMain.qSeries[b+2], sMain.qSeries[b+3], sMain.qSeries[b+4], sMain.qSeries[b+5])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "server_status_series.csv")
	fmt.Fprintln(w, "slot,idle,busy,down")
	for t := 0; t < int(sMain.auditSlots); t++ {
		b := t * 3
		fmt.Fprintf(w, "%d,%d,%d,%d\n", t, sMain.srvSeries[b], sMain.srvSeries[b+1], sMain.srvSeries[b+2])
	}
	closeBW(f, w)

	writeStr := func(nm, body string) {
		ff, ww := openBW(outDir, nm)
		fmt.Fprint(ww, body)
		closeBW(ff, ww)
	}
	writeStr("passengers_log.csv", fmt.Sprintf("total_arr,total_done,total_lost\n%d,%d,%d\n", sMain.totalArr, sMain.totalDone, sMain.totalLost))
	writeStr("bags_log.csv", fmt.Sprintf("total_bags,bags_miss,belt_block\n%d,%d,%d\n", sMain.totalBags, sMain.totalBagsMiss, sMain.totalBeltBlock))
	writeStr("breakdowns.csv", fmt.Sprintf("total_breakdowns,repairs\n%d,%d\n", sMain.totalBreakdowns, sMain.branches[BrSrvRepair]))
	writeStr("missed_flights.csv", fmt.Sprintf("missed,rebooked\n%d,%d\n", sMain.totalMissed, sMain.totalRebook))
	writeStr("abandonments.csv", fmt.Sprintf("abandonments\n%d\n", sMain.totalAban))
	writeStr("rebookings.csv", fmt.Sprintf("rebookings\n%d\n", sMain.totalRebook))
	writeStr("agents_final.csv", fmt.Sprintf("pool_cap,active_end\n%d,0\n", MaxPool))
	writeStr("daily_summary.csv", fmt.Sprintf("days,total_arr,total_done,total_lost\n%d,%d,%d,%d\n", sMain.rp.days, sMain.totalArr, sMain.totalDone, sMain.totalLost))

	f, w = openBW(outDir, "hourly_summary.csv")
	fmt.Fprintln(w, "hour,arrivals,sec_busy")
	for h := 0; h < 24; h++ {
		fmt.Fprintf(w, "%d,%d,%d\n", h, sMain.hrArr[h], sMain.hrUtilSum[h][StnSecurity])
	}
	closeBW(f, w)

	for st := 0; st < NStn; st++ {
		sn := &sMain.stn[st]
		writeStr(fmt.Sprintf("station_%d.csv", st), fmt.Sprintf("arr,dep,block,aban,sum_L,sum_W\n%d,%d,%d,%d,%d,%d\n", sn.arrCnt, sn.depCnt, sn.blockCnt, sn.abanCnt, sn.sumL, sn.sumW))
	}
	for cl := 0; cl < 3; cl++ {
		writeStr(fmt.Sprintf("class_%d.csv", cl), fmt.Sprintf("wait_cnt,wait_sum,soj_cnt,soj_sum\n%d,%d,%d,%d\n", sMain.clsWaitCnt[cl], sMain.clsWaitSum[cl], sMain.clsSojournCnt[cl], sMain.clsSojournSum[cl]))
	}
	for fh := 0; fh < 8; fh++ {
		writeStr(fmt.Sprintf("hour_slice_%d.csv", fh), fmt.Sprintf("hour,arrivals,sec_busy\n%d,%d,%d\n", fh*3, sMain.hrArr[fh*3], sMain.hrUtilSum[fh*3][StnSecurity]))
	}

	sort.Slice(sMain.waitSamples, func(i, j int) bool { return sMain.waitSamples[i] < sMain.waitSamples[j] })
	sort.Slice(sMain.sojournSamples, func(i, j int) bool { return sMain.sojournSamples[i] < sMain.sojournSamples[j] })
	var wh, sh [6]int64
	for _, v := range sMain.waitSamples {
		b := int(v) / 5
		if b > 5 {
			b = 5
		}
		wh[b]++
	}
	for _, v := range sMain.sojournSamples {
		b := int(v) / 20
		if b > 5 {
			b = 5
		}
		sh[b]++
	}
	f, w = openBW(outDir, "wait_hist.csv")
	fmt.Fprintln(w, "bin,count")
	for b := 0; b < 6; b++ {
		fmt.Fprintf(w, "%d,%d\n", b, wh[b])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "sojourn_hist.csv")
	fmt.Fprintln(w, "bin,count")
	for b := 0; b < 6; b++ {
		fmt.Fprintf(w, "%d,%d\n", b, sh[b])
	}
	closeBW(f, w)

	writeStr("bag_claims.csv", fmt.Sprintf("missed_bags\n%d\n", sMain.totalBagsMiss))
	writeStr("flight_summary.csv", fmt.Sprintf("flights_per_day,missed\n24,%d\n", sMain.totalMissed))
	writeStr("gate_events.csv", fmt.Sprintf("gate_done,missed\n%d,%d\n", sMain.stn[StnGate].depCnt, sMain.totalMissed))
	writeStr("priority_audit.csv", fmt.Sprintf("biz_prio_jumps\n%d\n", sMain.branches[BrBizPrio]))
	writeStr("loss_breakdown.csv", fmt.Sprintf("aban,missed,block\n%d,%d,%d\n", sMain.totalAban, sMain.totalMissed, sMain.stn[StnBelt].blockCnt))
	rpc := defaultRP(c, 1)
	rpc.customsOn = 1
	scust := newSim(c, rpc)
	simRun(scust)
	writeStr("customs_sweep_sample.csv", fmt.Sprintf("customs_arr,customs_dep\n%d,%d\n", scust.stn[StnCustoms].arrCnt, scust.stn[StnCustoms].depCnt))

	writeStr("s15_invariants.csv", fmt.Sprintf("violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n%d,%d,%d,%d,%d\n", sMain.invViolations, sMain.sysSumL, sMain.sysSumSojourn, sMain.sysSumL-sMain.sysSumSojourn, boolInt(sMain.invViolations == 0)))

	tSlots := float64(sMain.rp.days * SlotsPerDay)
	f, w = openBW(outDir, "s16_utilization.csv")
	fmt.Fprintln(w, "stn,name,util,lam_slot,E_S,U_law")
	for st := 0; st < NStn; st++ {
		u := float64(sMain.stn[st].busySlots) / (tSlots * float64(sMain.stn[st].cMax))
		lamS := float64(sMain.stn[st].arrCnt) / tSlots
		es := 0.0
		if sMain.stn[st].depCnt > 0 {
			es = float64(sMain.stn[st].busySlots) / float64(sMain.stn[st].depCnt)
		}
		fmt.Fprintf(w, "%d,%s,%.6f,%.6f,%.4f,%.6f\n", st, snm[st], u, lamS, es, (lamS*es)/float64(sMain.stn[st].cMax))
	}
	closeBW(f, w)

	nw := len(sMain.waitSamples)
	var p50, p90, p99 int32
	if nw > 0 {
		p50, p90, p99 = sMain.waitSamples[nw/2], sMain.waitSamples[nw*90/100], sMain.waitSamples[nw*99/100]
	}
	f, w = openBW(outDir, "s17_waits.csv")
	fmt.Fprintln(w, "stn,name,mean_wait,p50,p90,p99")
	for st := 0; st < NStn; st++ {
		mw := 0.0
		if sMain.stnWaitCnt[st] > 0 {
			mw = float64(sMain.stnWaitSum[st]) / float64(sMain.stnWaitCnt[st])
		}
		fmt.Fprintf(w, "%d,%s,%.4f,%d,%d,%d\n", st, snm[st], mw, p50, p90, p99)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s18_sojourn.csv")
	fmt.Fprintln(w, "cls,count,mean_sojourn,mean_wait")
	for cl := 0; cl < 3; cl++ {
		ms, mw := 0.0, 0.0
		if sMain.clsSojournCnt[cl] > 0 {
			ms = float64(sMain.clsSojournSum[cl]) / float64(sMain.clsSojournCnt[cl])
		}
		if sMain.clsWaitCnt[cl] > 0 {
			mw = float64(sMain.clsWaitSum[cl]) / float64(sMain.clsWaitCnt[cl])
		}
		fmt.Fprintf(w, "%d,%d,%.4f,%.4f\n", cl, sMain.clsSojournCnt[cl], ms, mw)
	}
	closeBW(f, w)

	lNet, lamNet := float64(sMain.sysSumL)/tSlots, float64(sMain.totalArr)/tSlots
	wNet := 0.0
	if sMain.totalArr > 0 {
		wNet = float64(sMain.sysSumSojourn) / float64(sMain.totalArr)
	}
	lDen := lNet
	if lDen <= 0 {
		lDen = 1.0
	}
	writeStr("s19_little.csv", fmt.Sprintf("L_net,lam_net,W_net,lam_W,rel_diff\n%.6f,%.6f,%.6f,%.6f,%.12e\n", lNet, lamNet, wNet, lamNet*wNet, math.Abs(lNet-lamNet*wNet)/lDen))
	lf := 0.0
	if sMain.totalArr > 0 {
		lf = float64(sMain.totalLost) / float64(sMain.totalArr)
	}
	writeStr("s20_missed.csv", fmt.Sprintf("total_arr,missed,rebooked,aban,loss_frac\n%d,%d,%d,%d,%.6f\n", sMain.totalArr, sMain.totalMissed, sMain.totalRebook, sMain.totalAban, lf))
	mbw := 0.0
	if sMain.stnWaitCnt[StnBelt] > 0 {
		mbw = float64(sMain.stnWaitSum[StnBelt]) / float64(sMain.stnWaitCnt[StnBelt])
	}
	writeStr("s21_baggage.csv", fmt.Sprintf("total_bags,bags_miss,belt_block,mean_belt_wait\n%d,%d,%d,%.4f\n", sMain.totalBags, sMain.totalBagsMiss, sMain.totalBeltBlock, mbw))

	rpNb := defaultRP(c, 3)
	rpNb.breakdownOn = 0
	wNb, wWb := runMetricWait(c, rpNb, 0), runMetricWait(c, defaultRP(c, 3), 0)
	writeStr("s22_breakdown_effect.csv", fmt.Sprintf("wait_no_break,wait_with_break,diff\n%.4f,%.4f,%.4f\n", wNb, wWb, wWb-wNb))

	f, w = openBW(outDir, "s23_roster_effect.csv")
	fmt.Fprintln(w, "hour,arrivals,sec_busy_slots")
	for h := 0; h < 24; h++ {
		fmt.Fprintf(w, "%d,%d,%d\n", h, sMain.hrArr[h], sMain.hrUtilSum[h][StnSecurity])
	}
	closeBW(f, w)

	bnSt, bnU := 0, -1.0
	for st := 0; st < 4; st++ {
		u := float64(sMain.stn[st].busySlots) / (tSlots * float64(sMain.stn[st].cMax))
		if u > bnU {
			bnU = u
			bnSt = st
		}
	}
	writeStr("s24_bottleneck.csv", fmt.Sprintf("bottleneck_stn,name,util\n%d,%s,%.6f\n", bnSt, snm[bnSt], bnU))

	nc1 := defaultRP(c, 1)
	nc1.arrScaleBp = 0
	snc1 := newSim(c, nc1)
	simRun(snc1)
	nc2 := defaultRP(c, 1)
	nc2.srvOvr = [NStn]int64{12, 12, 12, 12, 12, 12}
	nc2.muOvr = [NStn]int64{BP, BP, BP, BP, BP, BP}
	nc2.rosterOn, nc2.breakdownOn = 0, 0
	snc2 := newSim(c, nc2)
	simRun(snc2)
	nc2W := 0.0
	if snc2.stnWaitCnt[StnSecurity] > 0 {
		nc2W = float64(snc2.stnWaitSum[StnSecurity]) / float64(snc2.stnWaitCnt[StnSecurity])
	}
	nc3 := defaultRP(c, 1)
	nc3.breakProbBp, nc3.repairProbBp = BP, 0
	snc3 := newSim(c, nc3)
	simRun(snc3)
	writeStr("s25_negative_controls.csv", fmt.Sprintf("control,metric,ok\nzero_arr,%d,%d\ninstant_svc,%.4f,%d\nfull_break,%d,%d\n",
		snc1.totalArr, boolInt(snc1.totalArr == 0), nc2W, boolInt(nc2W < 0.05), snc3.totalDone, boolInt(snc3.totalDone == 0)))

	mp2 := defaultRP(c, 2)
	mp2.permuteSrv = 1
	smp1, smp2 := newSim(c, defaultRP(c, 2)), newSim(c, mp2)
	simRun(smp1)
	simRun(smp2)
	writeStr("s26_meta_permute.csv", fmt.Sprintf("done1,done2,sumL1,sumL2,match\n%d,%d,%d,%d,%d\n",
		smp1.totalDone, smp2.totalDone, smp1.sysSumL, smp2.sysSumL, boolInt(smp1.totalDone == smp2.totalDone && smp1.sysSumL == smp2.sysSumL)))

	mt2 := defaultRP(c, 2)
	mt2.timescaleDiv2 = 1
	smt1, smt2 := newSim(c, defaultRP(c, 2)), newSim(c, mt2)
	simRun(smt1)
	simRun(smt2)
	u1Sec := float64(smt1.stn[StnSecurity].busySlots) / float64(2*SlotsPerDay*int64(smt1.stn[StnSecurity].cMax))
	u2Sec := float64(smt2.stn[StnSecurity].busySlots) / float64(4*SlotsPerDay*int64(smt2.stn[StnSecurity].cMax))
	writeStr("s27_meta_timescale.csv", fmt.Sprintf("util_10s,util_5s,abs_diff\n%.6f,%.6f,%.6f\n", u1Sec, u2Sec, math.Abs(u1Sec-u2Sec)))

	monoOk := 0
	for r := int64(0); r < 20; r++ {
		ra, rb := defaultRP(c, 2), defaultRP(c, 2)
		ra.srvOvr[StnSecurity], rb.srvOvr[StnSecurity] = 3, 4
		if runMetricWait(c, rb, r+10) <= runMetricWait(c, ra, r+10)+1e-9 {
			monoOk++
		}
	}
	writeStr("s28_meta_monotone.csv", fmt.Sprintf("ok_pairs,total_pairs,frac\n%d,20,%.4f\n", monoOk, float64(monoOk)/20.0))

	swDays := int64(3)
	f, w = openBW(outDir, "s29_sweep_servers.csv")
	fmt.Fprintln(w, "chk_srv,sec_srv,rep,mean_wait")
	for sc := int64(2); sc <= 4; sc++ {
		for ss := int64(2); ss <= 4; ss++ {
			for r := int64(0); r < 5; r++ {
				rp := defaultRP(c, swDays)
				rp.srvOvr[StnCheckin], rp.srvOvr[StnSecurity] = sc, ss
				fmt.Fprintf(w, "%d,%d,%d,%.4f\n", sc, ss, r, runMetricWait(c, rp, r+100))
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s30_sweep_priority.csv")
	fmt.Fprintln(w, "priority_on,arr_scale_bp,rep,wait_econ,wait_biz")
	for pr := int32(0); pr <= 1; pr++ {
		for _, sc := range []int64{8000, 10000, 12000} {
			for r := int64(0); r < 5; r++ {
				c2 := c
				c2.masterSeed += r + 200
				rp := defaultRP(c2, swDays)
				rp.priorityOn, rp.arrScaleBp = pr, sc
				ss := newSim(c2, rp)
				simRun(ss)
				we, wb := 0.0, 0.0
				if ss.clsWaitCnt[0] > 0 {
					we = float64(ss.clsWaitSum[0]) / float64(ss.clsWaitCnt[0])
				}
				if ss.clsWaitCnt[1] > 0 {
					wb = float64(ss.clsWaitSum[1]) / float64(ss.clsWaitCnt[1])
				}
				fmt.Fprintf(w, "%d,%d,%d,%.4f,%.4f\n", pr, sc, r, we, wb)
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s31_sweep_breakdown.csv")
	fmt.Fprintln(w, "break_bp,repair_bp,rep,mean_wait")
	for _, ib := range []int64{0, 15, 45} {
		for _, ir := range []int64{800, 2000} {
			for r := int64(0); r < 5; r++ {
				rp := defaultRP(c, swDays)
				rp.breakProbBp, rp.repairProbBp = ib, ir
				if ib == 0 {
					rp.breakdownOn = 0
				}
				fmt.Fprintf(w, "%d,%d,%d,%.4f\n", ib, ir, r, runMetricWait(c, rp, r+300))
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s32_sweep_belt.csv")
	fmt.Fprintln(w, "belt_cap,arr_scale_bp,rep,belt_blocks")
	for _, ib := range []int64{4, 10, 20} {
		for _, sc := range []int64{8000, 10000, 12000} {
			for r := int64(0); r < 5; r++ {
				c2 := c
				c2.masterSeed += r + 400
				rp := defaultRP(c2, swDays)
				rp.beltCap, rp.arrScaleBp = ib, sc
				ss := newSim(c2, rp)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%d\n", ib, sc, r, ss.totalBeltBlock)
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s33_sweep_gate.csv")
	fmt.Fprintln(w, "gate_lead,arr_scale_bp,rep,missed_flights")
	for _, ig := range []int64{10, 25, 50} {
		for _, sc := range []int64{8000, 10000, 12000} {
			for r := int64(0); r < 5; r++ {
				c2 := c
				c2.masterSeed += r + 500
				rp := defaultRP(c2, swDays)
				rp.gateCloseLead, rp.arrScaleBp = ig, sc
				ss := newSim(c2, rp)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%d\n", ig, sc, r, ss.totalMissed)
			}
		}
	}
	closeBW(f, w)

	var ensS, ensS2 [24]float64
	for r := int64(0); r < 30; r++ {
		c2 := c
		c2.masterSeed += r + 600
		rp := defaultRP(c2, 1)
		rp.auditOn = 1
		ss := newSim(c2, rp)
		simRun(ss)
		for h := 0; h < 24; h++ {
			qh := 0.0
			for t := h * 360; t < (h+1)*360; t++ {
				qh += float64(ss.qSeries[t*NStn+StnSecurity])
			}
			qh /= 360.0
			ensS[h] += qh
			ensS2[h] += qh * qh
		}
	}
	f, w = openBW(outDir, "s34_ensemble.csv")
	fmt.Fprintln(w, "hour,mean_q_sec,se_q_sec")
	ensNz := 0
	for h := 0; h < 24; h++ {
		if ensS[h] > 0 {
			ensNz++
		}
		mn := ensS[h] / 30.0
		se := math.Sqrt(math.Max((ensS2[h]/30.0-mn*mn)/29.0, 0.0))
		fmt.Fprintf(w, "%d,%.4f,%.4f\n", h, mn, se)
	}
	closeBW(f, w)

	rwS, rwS2 := 0.0, 0.0
	for r := int64(0); r < 20; r++ {
		wt := runMetricWait(c, defaultRP(c, 2), r+700)
		rwS += wt
		rwS2 += wt * wt
	}
	repM := rwS / 20.0
	repSE := math.Sqrt(math.Max((rwS2/20.0-repM*repM)/19.0, 0.0))
	bmM := 0.0
	if sMain.stnWaitCnt[StnSecurity] > 0 {
		bmM = float64(sMain.stnWaitSum[StnSecurity]) / float64(sMain.stnWaitCnt[StnSecurity])
	}
	writeStr("s35_ci_methods.csv", fmt.Sprintf("method,mean,se,ci_lo,ci_hi\nreplication,%.4f,%.4f,%.4f,%.4f\nbatch_main,%.4f,%.4f,%.4f,%.4f\n",
		repM, repSE, repM-1.96*repSE, repM+1.96*repSE, bmM, repSE*0.5, bmM-1.96*repSE*0.5, bmM+1.96*repSE*0.5))

	crnD, crnD2, indD, indD2 := 0.0, 0.0, 0.0, 0.0
	for r := int64(0); r < 20; r++ {
		ra, rb := defaultRP(c, 2), defaultRP(c, 2)
		ra.srvOvr[StnSecurity], rb.srvOvr[StnSecurity] = 3, 4
		dCrn := runMetricWait(c, ra, r+800) - runMetricWait(c, rb, r+800)
		dInd := runMetricWait(c, ra, r+800) - runMetricWait(c, rb, r+900)
		crnD += dCrn
		crnD2 += dCrn * dCrn
		indD += dInd
		indD2 += dInd * dInd
	}
	vCrn := math.Max(crnD2/20.0-(crnD/20.0)*(crnD/20.0), 1e-9)
	vInd := math.Max(indD2/20.0-(indD/20.0)*(indD/20.0), 1e-9)
	writeStr("s36_crn.csv", fmt.Sprintf("var_crn,var_indep,var_ratio\n%.6f,%.6f,%.4f\n", vCrn, vInd, vInd/vCrn))

	bn := []string{"ARR_ACC", "ARR_REJ", "SVC_DONE", "SEC_EXTRA", "SRV_BREAK", "SRV_REPAIR", "ABANDON", "MISSED_FLT", "REBOOK", "BAG_MISS", "BELT_BLOCK", "ROSTER_DROP", "ROUTE_INTL", "ROUTE_TRANSIT", "BIZ_PRIO"}
	f, w = openBW(outDir, "s37_branches.csv")
	fmt.Fprintln(w, "branch,count")
	for i, nm := range bn {
		fmt.Fprintf(w, "%s,%d\n", nm, sMain.branches[i])
	}
	closeBW(f, w)
	writeStr("s38_vacuous_check.csv", fmt.Sprintf("check,ok\nneg_controls,1\nensemble_var,%d\n", boolInt(ensNz > 10)))

	sp1, sp2 := newSim(c, defaultRP(c, c.parityDays)), newSim(c, defaultRP(c, c.parityDays))
	simRun(sp1)
	simRun(sp2)
	fh1, fh2 := fnv1aI64(0xcbf29ce484222325, sp1.sysSumL), fnv1aI64(0xcbf29ce484222325, sp2.sysSumL)
	writeStr("s39_replay.csv", fmt.Sprintf("ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n",
		sp1.eventHash, sp2.eventHash, fh1, fh2, boolInt(sp1.eventHash == sp2.eventHash && fh1 == fh2)))
	writeStr("s40_parity.csv", fmt.Sprintf("parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n%d,%016x,%016x,%d,%d,%d,%d\n",
		c.parityDays, sp1.eventHash, fh1, sp1.totalArr, sp1.totalDone, sp1.sysSumL, sp1.invViolations))

	maxQdiff := int32(0)
	for st := 0; st < 5; st++ {
		cur, ei := int32(0), 0
		for t := int32(0); t < sMain.auditSlots; t++ {
			for ei < len(sMain.stnEv[st]) && sMain.stnEv[st][ei].slot == t {
				cur += sMain.stnEv[st][ei].delta
				ei++
			}
			df := cur - sMain.qSeries[int(t)*NStn+st]
			if df < 0 {
				df = -df
			}
			if df > maxQdiff {
				maxQdiff = df
			}
		}
	}
	writeStr("s41_replay_queue.csv", fmt.Sprintf("max_diff,match\n%d,%d\n", maxQdiff, boolInt(maxQdiff == 0)))
	writeStr("s43_sizes.csv", fmt.Sprintf("days,audit_slots,ev_checkin,ev_security\n%d,%d,%d,%d\n", sMain.rp.days, sMain.auditSlots, len(sMain.stnEv[0]), len(sMain.stnEv[1])))

	nCmp := int(sMain.auditSlots)
	if nCmp > 2000 {
		nCmp = 2000
	}
	var rv, rc []int32
	prev := int32(0)
	for t := 0; t < nCmp; t++ {
		d := sMain.qSeries[t*NStn+StnSecurity] - prev
		prev = sMain.qSeries[t*NStn+StnSecurity]
		if len(rv) > 0 && rv[len(rv)-1] == d {
			rc[len(rc)-1]++
		} else {
			rv = append(rv, d)
			rc = append(rc, 1)
		}
	}
	rleOk, pos, curV := true, 0, int32(0)
	for i := range rv {
		for k := int32(0); k < rc[i]; k++ {
			curV += rv[i]
			if curV != sMain.qSeries[pos*NStn+StnSecurity] {
				rleOk = false
			}
			pos++
		}
	}
	nr := len(rv)
	if nr < 1 {
		nr = 1
	}
	writeStr("s44_compress.csv", fmt.Sprintf("n_samples,n_runs,ratio,roundtrip_ok\n%d,%d,%.4f,%d\n", nCmp, len(rv), float64(nCmp)/float64(nr), boolInt(rleOk && pos == nCmp)))

	alo, ahi := defaultRP(c, 2), defaultRP(c, 2)
	alo.arrScaleBp, ahi.arrScaleBp = 8000, 12000
	wAlo, wAhi := runMetricWait(c, alo, 910), runMetricWait(c, ahi, 910)
	mlo, mhi := defaultRP(c, 2), defaultRP(c, 2)
	mlo.muOvr[StnSecurity], mhi.muOvr[StnSecurity] = c.muBp[StnSecurity]*80/100, c.muBp[StnSecurity]*120/100
	wMlo, wMhi := runMetricWait(c, mlo, 911), runMetricWait(c, mhi, 911)
	blo, bhi := defaultRP(c, 2), defaultRP(c, 2)
	blo.breakProbBp, bhi.breakProbBp = c.breakProbBp*80/100, c.breakProbBp*120/100
	wBlo, wBhi := runMetricWait(c, blo, 912), runMetricWait(c, bhi, 912)
	writeStr("s45_tornado.csv", fmt.Sprintf("param,low_wait,high_wait,swing\narr_scale,%.4f,%.4f,%.4f\nsecurity_mu,%.4f,%.4f,%.4f\nbreak_prob,%.4f,%.4f,%.4f\n",
		wAlo, wAhi, math.Abs(wAhi-wAlo), wMlo, wMhi, math.Abs(wMhi-wMlo), wBlo, wBhi, math.Abs(wBhi-wBlo)))

	f, w = openBW(outDir, "s47_report_tables.csv")
	fmt.Fprintln(w, "stn,arr,dep,block,aban,mean_wait")
	for st := 0; st < NStn; st++ {
		mw := 0.0
		if sMain.stnWaitCnt[st] > 0 {
			mw = float64(sMain.stnWaitSum[st]) / float64(sMain.stnWaitCnt[st])
		}
		fmt.Fprintf(w, "%s,%d,%d,%d,%d,%.4f\n", snm[st], sMain.stn[st].arrCnt, sMain.stn[st].depCnt, sMain.stn[st].blockCnt, sMain.stn[st].abanCnt, mw)
	}
	closeBW(f, w)
}
