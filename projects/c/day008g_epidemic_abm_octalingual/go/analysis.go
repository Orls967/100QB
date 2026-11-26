package main

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
)

func openBW(dir, name string) (*os.File, *bufio.Writer) {
	f, err := os.Create(dir + "/" + name)
	if err != nil {
		panic(err)
	}
	return f, bufio.NewWriter(f)
}

func closeBW(f *os.File, w *bufio.Writer) {
	_ = w.Flush()
	_ = f.Close()
}

func boolInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

func netDegrees(s *Sim) (float64, float64, float64) {
	shh, ssw, scom := 0.0, 0.0, 0.0
	for i := range s.pop {
		p := &s.pop[i]
		h := int(p.hhID)
		shh += float64(s.hhOff[h+1] - s.hhOff[h] - 1)
		if p.ageGrp == 0 {
			ssw += 6.0
		} else if p.ageGrp == 1 {
			ssw += 5.0
		}
		if p.superFlag != 0 {
			scom += 12.0
		} else {
			scom += 3.0
		}
	}
	nf := float64(s.n)
	return shh / nf, ssw / nf, scom / nf
}

func writeS01ToS09(c Config, outDir string) {
	hrng := uint64(0xcbf29ce484222325)
	pairs := [3][2]uint32{{0, 0}, {1, 2}, {2, 1}}
	expHv := [3]uint64{0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb}
	rngOk := true
	for i := 0; i < 3; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(pairs[i][0], pairs[i][1])
		sd := splitmix64(&x)
		st := newXoshiro(sd)
		var buf [32]byte
		for j := 0; j < 4; j++ {
			binary.LittleEndian.PutUint64(buf[j*8:(j+1)*8], st.nextU64())
		}
		hv := fnv1aBytes(buf[:])
		if hv != expHv[i] {
			rngOk = false
		}
		hrng = fnv1aI64(hrng, int64(hv))
	}
	st0 := newSim(c, defaultRunParams(200, 10, 450, 0), false)
	simRun(st0)
	invOk := st0.invViolations == 0

	f, w := openBW(outDir, "s01_selftest.txt")
	fmt.Fprintf(w, "ok,h_rng,inv_ok,tot_trans\n%d,%016x,%d,%d\n", boolInt(rngOk && invOk), hrng, boolInt(invOk), st0.totalTrans)
	closeBW(f, w)

	f, w = openBW(outDir, "s02_config.csv")
	fmt.Fprintf(w, "key,value\nmaster_seed,%d\nparity_N,%d\nparity_days,%d\nmain_days,%d\nbeta_bp_init,%d\n",
		c.masterSeed, c.parityN, c.parityDays, c.mainDays, c.betaBpInit)
	closeBW(f, w)

	sp := newSim(c, defaultRunParams(1000, 1, 450, 0), false)
	f, w = openBW(outDir, "s03_population.csv")
	fmt.Fprintf(w, "age_grp,count\n0,%d\n1,%d\n2,%d\n", sp.ageTot[0], sp.ageTot[1], sp.ageTot[2])
	closeBW(f, w)

	f, w = openBW(outDir, "s03_households.csv")
	fmt.Fprintln(w, "hh_id,size")
	for h := 0; h < sp.nHH && h < 50; h++ {
		fmt.Fprintf(w, "%d,%d\n", h, sp.hhOff[h+1]-sp.hhOff[h])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s03_schools_work.csv")
	fmt.Fprintln(w, "sw_id,size")
	for sw := 0; sw < sp.nSW && sw < 50; sw++ {
		fmt.Fprintf(w, "%d,%d\n", sw, sp.swOff[sw+1]-sp.swOff[sw])
	}
	closeBW(f, w)

	chh, csw, ccom := netDegrees(sp)
	f, w = openBW(outDir, "s04_network.csv")
	fmt.Fprintf(w, "n_hh,n_sw,deg_hh,deg_sw,deg_comm\n%d,%d,%.4f,%.4f,%.4f\n", sp.nHH, sp.nSW, chh, csw, ccom)
	closeBW(f, w)

	f, w = openBW(outDir, "s05_contact_script.csv")
	fmt.Fprintln(w, "case_id,vax,mask,quar,eff_bp")
	for k := 0; k < 12; k++ {
		vax, mask, quar := k%3, (k/3)%2, (k/6)%2
		p := int64(450)
		if quar != 0 {
			p = 0
		}
		if mask != 0 {
			p = (p * (BP - c.maskEffBp)) / BP
		}
		if vax == 1 {
			p = (p * (BP - c.vaxEff1Bp)) / BP
		} else if vax == 2 {
			p = (p * (BP - c.vaxEff2Bp)) / BP
		}
		fmt.Fprintf(w, "%d,%d,%d,%d,%d\n", k, vax, mask, quar, p)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s06_progression_script.csv")
	fmt.Fprintln(w, "case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp")
	for k := 0; k < 12; k++ {
		ag, rej := k%3, (k/3)%2
		ph := 1800
		if ag == 0 {
			ph = 120
		} else if ag == 1 {
			ph = 450
		}
		pd := 2000
		if ag == 0 {
			pd = 100
		} else if ag == 1 {
			pd = 500
		}
		if rej != 0 {
			pd *= 2
		}
		fmt.Fprintf(w, "%d,%d,%d,%d,%d\n", k, ag, rej, ph, pd)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s07_interventions.csv")
	fmt.Fprintf(w, "name,start_day,param_bp\nschool_close,%d,10000\nlockdown,%d,5000\nmask,%d,%d\ntest,%d,%d\nvax,%d,%d\nvariant,%d,%d\n",
		c.schoolCloseDay, c.lockdownDay, c.maskDay, c.maskEffBp, c.testDay, c.testProbBp, c.vaxDay, c.vaxDailyBp, c.variantDay, c.variantMultBp)
	closeBW(f, w)

	f, w = openBW(outDir, "s08_calibration.csv")
	fmt.Fprintln(w, "beta_bp,r0_analytic,attack_rate,calibrated")
	for _, bv := range []int64{300, 380, 450, 520, 600} {
		sc := newSim(c, defaultRunParams(1000, 120, bv, 0), false)
		simRun(sc)
		ch, cs2, cc2 := netDegrees(sc)
		r0 := analyticalR0(c, bv, ch, cs2, cc2)
		ar := float64(sc.totalTrans) / float64(sc.n)
		cal := boolInt(r0 >= 1.5 && r0 <= 3.5 && ar >= 0.20 && ar <= 0.80)
		fmt.Fprintf(w, "%d,%.4f,%.4f,%d\n", bv, r0, ar, cal)
	}
	closeBW(f, w)

	s9 := newSim(c, defaultRunParams(500, 30, 450, 0), false)
	simRun(s9)
	f, w = openBW(outDir, "s09_baseline_small.csv")
	fmt.Fprintf(w, "N,days,total_trans,deaths,inv_ok\n500,30,%d,%d,%d\n", s9.totalTrans, s9.dailyCnt[29][StD], boolInt(s9.invViolations == 0))
	closeBW(f, w)
}

func writeDailyCSV(s *Sim, outDir, fname string) {
	f, w := openBW(outDir, fname)
	fmt.Fprintln(w, "day,S,E,Ia,Is,H,ICU,R,D")
	for d := 0; d < int(s.rp.days); d++ {
		fmt.Fprintf(w, "%d", d)
		for k := 0; k < 8; k++ {
			fmt.Fprintf(w, ",%d", s.dailyCnt[d][k])
		}
		fmt.Fprintln(w)
	}
	closeBW(f, w)
}

func writeTransCSV(s *Sim, outDir, fname string) {
	f, w := openBW(outDir, fname)
	fmt.Fprintln(w, "day,src_id,dst_id,layer,variant,gen_int")
	for _, r := range s.trans {
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d\n", r.day, r.src, r.dst, r.lay, r.vari, r.genInt)
	}
	closeBW(f, w)
}

func writeS10ToS37(sNo, sMain *Sim, c Config, outDir string, betaCal int64) {
	dCnt := int(sMain.rp.days)
	writeDailyCSV(sNo, outDir, "daily_counts_nointerv.csv")
	writeTransCSV(sNo, outDir, "transmissions_log_nointerv.csv")
	writeDailyCSV(sMain, outDir, "daily_counts.csv")
	writeTransCSV(sMain, outDir, "transmissions_log.csv")

	f, w := openBW(outDir, "incidence_by_age.csv")
	fmt.Fprintln(w, "day,age_0_17,age_18_59,age_60p")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d,%d\n", d, sMain.dailyIncAge[d][0], sMain.dailyIncAge[d][1], sMain.dailyIncAge[d][2])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "hospital_daily.csv")
	fmt.Fprintln(w, "day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d,%d,%d\n", d, sMain.dailyHosp[d][0], sMain.dailyHosp[d][1], sMain.dailyHosp[d][2], sMain.dailyHosp[d][3])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "deaths_by_age.csv")
	fmt.Fprintln(w, "day,d_0_17,d_18_59,d_60p")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d,%d\n", d, sMain.dailyDeathsAge[d][0], sMain.dailyDeathsAge[d][1], sMain.dailyDeathsAge[d][2])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "tests.csv")
	fmt.Fprintln(w, "day,tests,positives")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d\n", d, sMain.dailyTests[d][0], sMain.dailyTests[d][1])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "quarantine.csv")
	fmt.Fprintln(w, "day,quarantined")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d\n", d, sMain.dailyQuar[d])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "vaccination.csv")
	fmt.Fprintln(w, "day,dose1,dose2")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d\n", d, sMain.dailyVax[d][0], sMain.dailyVax[d][1])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "variant_share.csv")
	fmt.Fprintln(w, "day,v0_inf,v1_inf")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d\n", d, sMain.dailyVar[d][0], sMain.dailyVar[d][1])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "generation_interval.csv")
	fmt.Fprintln(w, "idx,gen_int")
	for i := 0; i < len(sMain.trans) && i < 2000; i++ {
		fmt.Fprintf(w, "%d,%d\n", i, sMain.trans[i].genInt)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "agents_final.csv")
	fmt.Fprintln(w, "id,age_grp,status,vax_doses,offspring,ever_inf")
	for i := 0; i < len(sMain.pop) && i < 2000; i++ {
		p := &sMain.pop[i]
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d\n", i, p.ageGrp, p.status, p.vaxDoses, p.offspring, p.everInf)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "progression_log.csv")
	fmt.Fprintln(w, "day,pid,from_st,to_st")
	for _, r := range sMain.progs {
		fmt.Fprintf(w, "%d,%d,%d,%d\n", r.day, r.pid, r.fromSt, r.toSt)
	}
	closeBW(f, w)

	for a := 0; a < 3; a++ {
		f, w = openBW(outDir, fmt.Sprintf("daily_age_%d.csv", a))
		fmt.Fprintln(w, "day,inc,deaths")
		for d := 0; d < dCnt; d++ {
			fmt.Fprintf(w, "%d,%d,%d\n", d, sMain.dailyIncAge[d][a], sMain.dailyDeathsAge[d][a])
		}
		closeBW(f, w)
	}
	for l, nm := range []string{"layer_hh.csv", "layer_sw.csv", "layer_comm.csv"} {
		f, w = openBW(outDir, nm)
		fmt.Fprintf(w, "layer,infections\n%d,%d\n", l, sMain.layerInf[l])
		closeBW(f, w)
	}
	f, w = openBW(outDir, "hospital_cap_nointerv.csv")
	fmt.Fprintln(w, "day,hosp_occ,icu_occ")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d,%d\n", d, sNo.dailyHosp[d][0], sNo.dailyHosp[d][1])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "deaths_nointerv.csv")
	fmt.Fprintln(w, "day,deaths")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d\n", d, sNo.dailyCnt[d][StD])
	}
	closeBW(f, w)

	var ohBins [6]int64
	for i := range sNo.pop {
		if sNo.pop[i].everInf != 0 {
			b := int(sNo.pop[i].offspring)
			if b > 5 {
				b = 5
			}
			ohBins[b]++
		}
	}
	f, w = openBW(outDir, "offspring_hist.csv")
	fmt.Fprintln(w, "bin,count")
	for b := 0; b < 6; b++ {
		fmt.Fprintf(w, "%d,%d\n", b, ohBins[b])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "rt_nointerv.csv")
	fmt.Fprintln(w, "day,inc")
	for d := 0; d < dCnt; d++ {
		fmt.Fprintf(w, "%d,%d\n", d, sNo.dailyIncAge[d][0]+sNo.dailyIncAge[d][1]+sNo.dailyIncAge[d][2])
	}
	closeBW(f, w)

	totViol := sNo.invViolations + sMain.invViolations
	f, w = openBW(outDir, "s12_invariants.csv")
	fmt.Fprintf(w, "violations,ok\n%d,%d\n", totViol, boolInt(totViol == 0))
	closeBW(f, w)

	f, w = openBW(outDir, "s13_rt.csv")
	fmt.Fprintln(w, "day,rt")
	for d := 0; d < dCnt; d++ {
		num, den := int64(0), int64(0)
		for k := 0; k < 5; k++ {
			if d >= k {
				num += sMain.dailyIncAge[d-k][0] + sMain.dailyIncAge[d-k][1] + sMain.dailyIncAge[d-k][2]
			}
		}
		for k := 4; k < 9; k++ {
			if d >= k {
				den += sMain.dailyIncAge[d-k][0] + sMain.dailyIncAge[d-k][1] + sMain.dailyIncAge[d-k][2]
			}
		}
		rt := 0.0
		if den > 0 {
			rt = float64(num) / float64(den)
		}
		fmt.Fprintf(w, "%d,%.4f\n", d, rt)
	}
	closeBW(f, w)

	chh, csw, ccom := netDegrees(sNo)
	r0Ana := analyticalR0(c, betaCal, chh, csw, ccom)
	iEarly := float64(sNo.dailyIncAge[3][0]+sNo.dailyIncAge[3][1]+sNo.dailyIncAge[3][2]+sNo.dailyIncAge[4][0]+sNo.dailyIncAge[4][1]+sNo.dailyIncAge[4][2]+1) / 2.0
	iMid := 0.0
	for d := 10; d < 16 && d < dCnt; d++ {
		iMid += float64(sNo.dailyIncAge[d][0] + sNo.dailyIncAge[d][1] + sNo.dailyIncAge[d][2])
	}
	iMid /= 6.0
	rGrow := 0.0
	if iMid > iEarly {
		rGrow = math.Log(iMid/iEarly) / 9.0
	}
	dbl := 0.0
	if rGrow > 0 {
		dbl = math.Log(2.0) / rGrow
	}
	f, w = openBW(outDir, "s14_growth.csv")
	fmt.Fprintf(w, "growth_r,doubling_days\n%.6f,%.4f\n", rGrow, dbl)
	closeBW(f, w)

	r0Emp := 1.0 + rGrow*4.5
	ratio := 0.0
	if r0Ana > 0 {
		ratio = r0Emp / r0Ana
	}
	f, w = openBW(outDir, "s15_r0.csv")
	fmt.Fprintf(w, "r0_analytic,r0_empirical,ratio\n%.4f,%.4f,%.4f\n", r0Ana, r0Emp, ratio)
	closeBW(f, w)

	zAna := solveFinalSize(r0Ana)
	herd := 0.0
	if r0Ana > 1.0 {
		herd = 1.0 - 1.0/r0Ana
	}
	f, w = openBW(outDir, "s16_final_size.csv")
	fmt.Fprintf(w, "z_analytic,herd_thresh,attack_nointerv,attack_interv\n%.4f,%.4f,%.4f,%.4f\n",
		zAna, herd, float64(sNo.totalTrans)/float64(sNo.n), float64(sMain.totalTrans)/float64(sMain.n))
	closeBW(f, w)

	f, w = openBW(outDir, "s17_attack.csv")
	fmt.Fprintf(w, "category,count,frac\nHH,%d,%.4f\nSW,%d,%.4f\nCOMM,%d,%.4f\nAGE0,%d,%.4f\nAGE1,%d,%.4f\nAGE2,%d,%.4f\n",
		sMain.layerInf[0], float64(sMain.layerInf[0])/float64(sMain.n),
		sMain.layerInf[1], float64(sMain.layerInf[1])/float64(sMain.n),
		sMain.layerInf[2], float64(sMain.layerInf[2])/float64(sMain.n),
		sMain.ageInf[0], float64(sMain.ageInf[0])/float64(max64(sMain.ageTot[0], 1)),
		sMain.ageInf[1], float64(sMain.ageInf[1])/float64(max64(sMain.ageTot[1], 1)),
		sMain.ageInf[2], float64(sMain.ageInf[2])/float64(max64(sMain.ageTot[2], 1)))
	closeBW(f, w)

	m1, m2, nInf := 0.0, 0.0, int64(0)
	for i := range sNo.pop {
		if sNo.pop[i].everInf != 0 {
			o := float64(sNo.pop[i].offspring)
			m1 += o
			m2 += o * o
			nInf++
		}
	}
	if nInf > 0 {
		m1 /= float64(nInf)
		m2 /= float64(nInf)
	}
	varO := m2 - m1*m1
	kDisp := 10.0
	if varO > m1 && m1 > 0 {
		kDisp = (m1 * m1) / (varO - m1)
	}
	f, w = openBW(outDir, "s18_offspring.csv")
	fmt.Fprintf(w, "n_inf,mean_offspring,var_offspring,k_disp\n%d,%.4f,%.4f,%.4f\n", nInf, m1, varO, kDisp)
	closeBW(f, w)

	pkHNo, pkINo, pkHM, pkIM := int64(0), int64(0), int64(0), int64(0)
	for d := 0; d < dCnt; d++ {
		if sNo.dailyHosp[d][0] > pkHNo {
			pkHNo = sNo.dailyHosp[d][0]
		}
		if sNo.dailyHosp[d][1] > pkINo {
			pkINo = sNo.dailyHosp[d][1]
		}
		if sMain.dailyHosp[d][0] > pkHM {
			pkHM = sMain.dailyHosp[d][0]
		}
		if sMain.dailyHosp[d][1] > pkIM {
			pkIM = sMain.dailyHosp[d][1]
		}
	}
	f, w = openBW(outDir, "s19_hospital.csv")
	fmt.Fprintf(w, "hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n%d,%d,%d,%d,%d,%d,%d,%d\n",
		sMain.hospCap, sMain.icuCap, pkHNo, pkINo, sNo.branches[9]+sNo.branches[11], pkHM, pkIM, sMain.branches[9]+sMain.branches[11])
	closeBW(f, w)

	f, w = openBW(outDir, "s20_ifr.csv")
	fmt.Fprintln(w, "age_grp,infections,deaths,ifr")
	for a := 0; a < 3; a++ {
		dt := int64(0)
		for d := 0; d < dCnt; d++ {
			dt += sMain.dailyDeathsAge[d][a]
		}
		ifr := 0.0
		if sMain.ageInf[a] > 0 {
			ifr = float64(dt) / float64(sMain.ageInf[a])
		}
		fmt.Fprintf(w, "%d,%d,%d,%.6f\n", a, sMain.ageInf[a], dt, ifr)
	}
	closeBW(f, w)

	tTot, tPos := int64(0), int64(0)
	for d := 0; d < dCnt; d++ {
		tTot += sMain.dailyTests[d][0]
		tPos += sMain.dailyTests[d][1]
	}
	f, w = openBW(outDir, "s21_testing.csv")
	fmt.Fprintf(w, "tests,positives,quar_traces\n%d,%d,%d\n", tTot, tPos, sMain.branches[17])
	closeBW(f, w)

	var vn, vi [3]int64
	for i := range sMain.pop {
		vd := int(sMain.pop[i].vaxDoses)
		if vd > 2 {
			vd = 2
		}
		vn[vd]++
		if sMain.pop[i].everInf != 0 {
			vi[vd]++
		}
	}
	f, w = openBW(outDir, "s22_vaccine.csv")
	fmt.Fprintln(w, "doses,n_people,n_inf,attack_rate")
	for d := 0; d < 3; d++ {
		ar := 0.0
		if vn[d] > 0 {
			ar = float64(vi[d]) / float64(vn[d])
		}
		fmt.Fprintf(w, "%d,%d,%d,%.6f\n", d, vn[d], vi[d], ar)
	}
	closeBW(f, w)

	v0, v1 := int64(0), int64(0)
	for d := 0; d < dCnt; d++ {
		v0 += sMain.dailyVar[d][0]
		v1 += sMain.dailyVar[d][1]
	}
	v1Sh := 0.0
	if v0+v1 > 0 {
		v1Sh = float64(v1) / float64(v0+v1)
	}
	f, w = openBW(outDir, "s23_variant.csv")
	fmt.Fprintf(w, "v0_total,v1_total,v1_share\n%d,%d,%.6f\n", v0, v1, v1Sh)
	closeBW(f, w)

	f, w = openBW(outDir, "s24_waning.csv")
	fmt.Fprintf(w, "waning_events,reinfections\n%d,%d\n", sMain.branches[14], sMain.totalReinf)
	closeBW(f, w)

	r1 := defaultRunParams(500, 30, betaCal, 0)
	r1.contactScaleBp = 0
	nc1 := newSim(c, r1, false)
	simRun(nc1)
	r2 := defaultRunParams(500, 30, 0, 0)
	r2.contactScaleBp = -1
	nc2 := newSim(c, r2, false)
	simRun(nc2)
	r3 := defaultRunParams(500, 30, betaCal, 1)
	r3.vaxPerfect = 1
	nc3 := newSim(c, r3, false)
	simRun(nc3)
	f, w = openBW(outDir, "s25_negative_controls.csv")
	fmt.Fprintf(w, "control,infections,ok\nzero_contact,%d,%d\nzero_seed,%d,%d\nperfect_vax,%d,%d\n",
		nc1.totalTrans, boolInt(nc1.totalTrans == 0),
		nc2.totalTrans, boolInt(nc2.totalTrans == 0),
		nc3.totalTrans, boolInt(nc3.totalTrans == 0))
	closeBW(f, w)

	nsw, dsw := int64(1000), int64(90)
	f, w = openBW(outDir, "s26_sweep_comp_test.csv")
	fmt.Fprintln(w, "comp_bp,test_bp,rep,attack_rate")
	for _, cv := range []int64{2000, 6000, 9500} {
		for _, tv := range []int64{1000, 5000, 9000} {
			for rp := int64(0); rp < 3; rp++ {
				c2 := c
				c2.masterSeed += rp + 1
				r := defaultRunParams(nsw, dsw, betaCal, 1)
				r.complianceOvr = cv
				r.testOvr = tv
				ss := newSim(c2, r, false)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%.6f\n", cv, tv, rp, float64(ss.totalTrans)/float64(nsw))
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s27_sweep_vax.csv")
	fmt.Fprintln(w, "vax_day,vax_cov_bp,rep,attack_rate")
	for _, vd := range []int64{10, 45, 80} {
		for _, vc := range []int64{2000, 5500, 9000} {
			for rp := int64(0); rp < 3; rp++ {
				c2 := c
				c2.masterSeed += rp + 10
				r := defaultRunParams(nsw, dsw, betaCal, 1)
				r.vaxDayOvr = vd
				r.vaxCovOvr = vc
				ss := newSim(c2, r, false)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%.6f\n", vd, vc, rp, float64(ss.totalTrans)/float64(nsw))
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s28_sweep_school_lock.csv")
	fmt.Fprintln(w, "school_day,lock_day,rep,attack_rate")
	for _, sc := range []int64{15, 200} {
		for _, lk := range []int64{15, 45, 200} {
			for rp := int64(0); rp < 3; rp++ {
				c2 := c
				c2.masterSeed += rp + 20
				r := defaultRunParams(nsw, dsw, betaCal, 1)
				r.schoolOvr = sc
				r.lockOvr = lk
				ss := newSim(c2, r, false)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%.6f\n", sc, lk, rp, float64(ss.totalTrans)/float64(nsw))
			}
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s29_sweep_variant_mask.csv")
	fmt.Fprintln(w, "var_day,mask_day,rep,attack_rate")
	for _, vr := range []int64{15, 45, 200} {
		for _, mk := range []int64{10, 40, 200} {
			for rp := int64(0); rp < 3; rp++ {
				c2 := c
				c2.masterSeed += rp + 30
				r := defaultRunParams(nsw, dsw, betaCal, 1)
				r.varDayOvr = vr
				r.maskOvr = mk
				ss := newSim(c2, r, false)
				simRun(ss)
				fmt.Fprintf(w, "%d,%d,%d,%.6f\n", vr, mk, rp, float64(ss.totalTrans)/float64(nsw))
			}
		}
	}
	closeBW(f, w)

	var ensS, ensS2 [90]float64
	for rp := int64(0); rp < 20; rp++ {
		c2 := c
		c2.masterSeed += rp + 100
		ss := newSim(c2, defaultRunParams(nsw, dsw, betaCal, 1), false)
		simRun(ss)
		for d := 0; d < 90; d++ {
			inc := float64(ss.dailyIncAge[d][0] + ss.dailyIncAge[d][1] + ss.dailyIncAge[d][2])
			ensS[d] += inc
			ensS2[d] += inc * inc
		}
	}
	f, w = openBW(outDir, "s30_ensemble.csv")
	fmt.Fprintln(w, "day,mean_inc,se_inc")
	ensNz := 0
	for d := 0; d < 90; d++ {
		if ensS[d] > 0 {
			ensNz++
		}
		mn := ensS[d] / 20.0
		v := (ensS2[d]/20.0 - mn*mn) / 19.0
		if v < 0 {
			v = 0
		}
		fmt.Fprintf(w, "%d,%.4f,%.4f\n", d, mn, math.Sqrt(v))
	}
	closeBW(f, w)

	bn := []string{"INF_HH", "INF_SW", "INF_COMM", "FAIL_VAX", "FAIL_MASK", "FAIL_QUAR", "PROG_IA", "PROG_IS",
		"HOSP_ADM", "HOSP_REJ", "ICU_ADM", "ICU_REJ", "DEATH", "RECOVER", "WANING",
		"TEST_POS", "TEST_NEG", "QUAR_TRACE", "VAX_DOSE", "VAX_OUT", "VAR_IMPORT"}
	f, w = openBW(outDir, "s31_branches.csv")
	fmt.Fprintln(w, "branch,count")
	for i, nm := range bn {
		fmt.Fprintf(w, "%s,%d\n", nm, sMain.branches[i])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s32_vacuous_check.csv")
	fmt.Fprintf(w, "check,ok\nneg_controls,%d\nensemble_var,%d\n",
		boolInt(nc1.totalTrans == 0 && nc2.totalTrans == 0 && nc3.totalTrans == 0), boolInt(ensNz > 5))
	closeBW(f, w)

	sp1 := newSim(c, defaultRunParams(c.parityN, c.parityDays, betaCal, 1), false)
	simRun(sp1)
	sp2 := newSim(c, defaultRunParams(c.parityN, c.parityDays, betaCal, 1), false)
	simRun(sp2)
	h1, h2 := uint64(0xcbf29ce484222325), uint64(0xcbf29ce484222325)
	for i := 0; i < int(sp1.n); i++ {
		h1 = fnv1aI64(h1, int64(sp1.pop[i].status))
		h2 = fnv1aI64(h2, int64(sp2.pop[i].status))
		h1 = fnv1aI64(h1, int64(sp1.pop[i].offspring))
		h2 = fnv1aI64(h2, int64(sp2.pop[i].offspring))
	}
	f, w = openBW(outDir, "s33_replay.csv")
	fmt.Fprintf(w, "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n",
		sp1.eventHash, sp2.eventHash, h1, h2, boolInt(sp1.eventHash == sp2.eventHash && h1 == h2))
	closeBW(f, w)

	f, w = openBW(outDir, "s34_parity.csv")
	fmt.Fprintf(w, "parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n%d,%d,%016x,%016x,%d,%d,%d\n",
		c.parityN, c.parityDays, sp1.eventHash, h1, sp1.totalTrans, sp1.dailyCnt[c.parityDays-1][StD], sp1.invViolations)
	closeBW(f, w)

	f, w = openBW(outDir, "s36_sizes.csv")
	fmt.Fprintf(w, "N,days,n_trans,n_progs\n%d,%d,%d,%d\n", sMain.n, sMain.rp.days, len(sMain.trans), len(sMain.progs))
	closeBW(f, w)

	cntR := [8]int64{sMain.n, 0, 0, 0, 0, 0, 0, 0}
	pi := 0
	maxDiff := int64(0)
	for d := 0; d < dCnt; d++ {
		for pi < len(sMain.progs) && sMain.progs[pi].day == int32(d) {
			cntR[sMain.progs[pi].fromSt]--
			cntR[sMain.progs[pi].toSt]++
			pi++
		}
		for k := 0; k < 8; k++ {
			diff := cntR[k] - sMain.dailyCnt[d][k]
			if diff < 0 {
				diff = -diff
			}
			if diff > maxDiff {
				maxDiff = diff
			}
		}
	}
	f, w = openBW(outDir, "s37_replay_counts.csv")
	fmt.Fprintf(w, "max_diff,match\n%d,%d\n", maxDiff, boolInt(maxDiff == 0))
	closeBW(f, w)
}
