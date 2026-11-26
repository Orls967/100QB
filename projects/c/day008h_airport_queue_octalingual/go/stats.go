package main

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"os"
)

func openBW(dir, nm string) (*os.File, *bufio.Writer) {
	f, err := os.Create(dir + "/" + nm)
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

func runSingleSim(rng *Xoshiro256, c, K int, lamBp, muBp, slots int64) (float64, float64, float64) {
	q := 0
	var arr, blk, sumL, sumB int64
	for t := int64(0); t < slots; t++ {
		a := int(rng.bernoulli(lamBp))
		m := q
		if c < m {
			m = c
		}
		d := 0
		for s := 0; s < m; s++ {
			if rng.bernoulli(muBp) != 0 {
				d++
			}
		}
		if q < K {
			q += a
			if a != 0 {
				arr++
			}
		} else if a != 0 {
			blk++
		}
		q -= d
		if t >= slots/5 {
			mb := q
			if c < mb {
				mb = c
			}
			sumL += int64(q)
			sumB += int64(mb)
		}
	}
	ns := float64(slots - slots/5)
	den := float64(arr + blk)
	if den <= 0 {
		den = 1.0
	}
	return float64(sumL) / ns, float64(sumB) / (ns * float64(c)), float64(blk) / den
}

func runTandemSim(rng *Xoshiro256, K1, K2 int, lamBp, mu1Bp, mu2Bp, slots int64) (float64, float64, float64, float64) {
	q1, q2 := 0, 0
	var sL1, sL2, sU1, sU2 int64
	for t := int64(0); t < slots; t++ {
		a := int(rng.bernoulli(lamBp))
		d1, d2 := 0, 0
		if q1 > 0 {
			d1 = int(rng.bernoulli(mu1Bp))
		}
		if q2 > 0 {
			d2 = int(rng.bernoulli(mu2Bp))
		}
		if q1 < K1 {
			q1 += a
		}
		q1 -= d1
		q2 -= d2
		if d1 != 0 && q2 < K2 {
			q2++
		}
		if t >= slots/5 {
			sL1 += int64(q1)
			sL2 += int64(q2)
			if q1 > 0 {
				sU1++
			}
			if q2 > 0 {
				sU2++
			}
		}
	}
	ns := float64(slots - slots/5)
	return float64(sL1) / ns, float64(sL2) / ns, float64(sU1) / ns, float64(sU2) / ns
}

func writeS01ToS13(c Config, outDir string) {
	hrng := uint64(0xcbf29ce484222325)
	pairs := [3][2]uint32{{0, 0}, {1, 2}, {2, 1}}
	expHv := [3]uint64{0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb}
	rngOk := true
	for i := 0; i < 3; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(pairs[i][0], pairs[i][1])
		st := newXoshiro(splitmix64(&x))
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
	s0 := newSim(c, defaultRP(c, 1))
	simRun(s0)
	invOk := s0.invViolations == 0

	f, w := openBW(outDir, "s01_selftest.txt")
	fmt.Fprintf(w, "ok,h_rng,inv_ok,tot_arr,tot_done\n%d,%016x,%d,%d,%d\n", boolInt(rngOk && invOk), hrng, boolInt(invOk), s0.totalArr, s0.totalDone)
	closeBW(f, w)

	f, w = openBW(outDir, "s02_config.csv")
	fmt.Fprintf(w, "key,value\nmaster_seed,%d\nslots_per_day,%d\nparity_days,%d\narr_scale_bp,%d\nbelt_cap,%d\n", c.masterSeed, c.slotsPerDay, c.parityDays, c.arrScaleBp, c.beltCap)
	closeBW(f, w)

	f, w = openBW(outDir, "s03_flights.csv")
	fmt.Fprintln(w, "flight_id,dep_slot,gate_close_slot,seats")
	for h := 0; h < 24; h++ {
		fmt.Fprintf(w, "%d,%d,%d,180\n", h, h*360+300, h*360+300-int(c.gateCloseLead))
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s03_arrival_profile.csv")
	fmt.Fprintln(w, "hour,arr_bp")
	for h := 0; h < 24; h++ {
		fmt.Fprintf(w, "%d,%d\n", h, s0.arrProf[h])
	}
	closeBW(f, w)

	snm := []string{"CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"}
	f, w = openBW(outDir, "s04_stations.csv")
	fmt.Fprintln(w, "stn_id,name,c_max,mu_bp,cap,has_prio")
	for i := 0; i < NStn; i++ {
		fmt.Fprintf(w, "%d,%s,%d,%d,%d,%d\n", i, snm[i], s0.stn[i].cMax, c.muBp[i], s0.stn[i].cap, s0.stn[i].hasPrio)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s04_roster.csv")
	fmt.Fprintln(w, "hour,checkin,security,imigrasi,gate,belt,customs")
	for h := 0; h < 24; h++ {
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d,%d\n", h, s0.roster[h][0], s0.roster[h][1], s0.roster[h][2], s0.roster[h][3], s0.roster[h][4], s0.roster[h][5])
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s05_queue_script.csv")
	fmt.Fprintln(w, "case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost")
	for k := 0; k < 12; k++ {
		q0, cap, brk, gc := (k%4)+1, 10, boolInt(k%3 == 2), boolInt(k == 11)
		if k >= 6 {
			cap = 3
		}
		acc := q0
		if cap < acc {
			acc = cap
		}
		srv := acc
		if gc != 0 {
			srv = 0
		} else if brk != 0 {
			srv = acc - 1
		}
		lost := q0 - acc
		if gc != 0 {
			lost = acc
		}
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", k, k%3, q0, cap, k%2, brk, gc, acc, srv, lost)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s06_bag_script.csv")
	fmt.Fprintln(w, "case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed")
	for k := 0; k < 12; k++ {
		nb, bl, cap := k%3, (k/3)*4, 10
		acc := nb
		if bl+nb > cap {
			acc = cap - bl
			if acc < 0 {
				acc = 0
			}
		}
		blk := nb - acc
		ms := boolInt(k%4 == 3 && acc > 0)
		fmt.Fprintf(w, "%d,%d,%d,%d,%d,%d,%d,%d\n", k, nb, bl, cap, c.bagMissBp, acc, blk, ms)
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s07_dtmc_single.csv")
	fmt.Fprintln(w, "c,K,lam,mu,p_block,L,L_q,U,W,W_q")
	lam, mu := 0.08, 0.12
	for _, cc := range []int{1, 2, 3} {
		for _, kk := range []int{5, 10, 60} {
			dr := solveDtmcSingle(cc, kk, lam, mu)
			fmt.Fprintf(w, "%d,%d,%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", cc, kk, lam, mu, dr.pBlock, dr.L, dr.Lq, dr.U, dr.W, dr.Wq)
		}
	}
	closeBW(f, w)

	t1, t2, pbt := solveDtmcTandem(6, 6, 0.08, 0.14, 0.12)
	f, w = openBW(outDir, "s08_dtmc_tandem.csv")
	fmt.Fprintf(w, "K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n6,6,0.0800,0.1400,0.1200,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", t1.pBlock, t1.L, t1.U, t2.pBlock, t2.L, t2.U, pbt)
	closeBW(f, w)

	cf, d60 := solveGeoClosedForm(lam, mu), solveDtmcSingle(1, 60, lam, mu)
	f, w = openBW(outDir, "s09_geo_closed_form.csv")
	fmt.Fprintf(w, "lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12e\n", lam, mu, cf.L, cf.W, cf.U, d60.L, math.Abs(cf.L-d60.L)/cf.L)
	closeBW(f, w)

	f, w = openBW(outDir, "s10_single_validation.csv")
	fmt.Fprintln(w, "c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se")
	valSlots := int64(20000)
	for ic, cc := range []int{1, 2, 3} {
		kk := 5
		if cc == 1 {
			kk = 10
		}
		dr := solveDtmcSingle(cc, kk, lam, mu)
		var sl, sl2, su, su2, sp, sp2 float64
		for r := uint32(0); r < 50; r++ {
			x := uint64(c.masterSeed) ^ fnv1aU32Pair(500+uint32(ic), r)
			rng := newXoshiro(splitmix64(&x))
			vl, vu, vp := runSingleSim(&rng, cc, kk, 800, 1200, valSlots)
			sl += vl
			sl2 += vl * vl
			su += vu
			su2 += vu * vu
			sp += vp
			sp2 += vp * vp
		}
		for _, item := range []struct {
			nm         string
			dt, s1, s2 float64
		}{{"L", dr.L, sl, sl2}, {"U", dr.U, su, su2}, {"p_block", dr.pBlock, sp, sp2}} {
			mn := item.s1 / 50.0
			se := math.Sqrt(math.Max((item.s2/50.0-mn*mn)/49.0, 1e-18))
			z := (mn - item.dt) / se
			fmt.Fprintf(w, "%d,%d,%s,%.6f,%.6f,%.6f,%.3f,%d\n", cc, kk, item.nm, item.dt, mn, se, z, boolInt(math.Abs(z) <= 3.0))
		}
	}
	closeBW(f, w)

	f, w = openBW(outDir, "s11_tandem_validation.csv")
	fmt.Fprintln(w, "metric,dtmc,sim_mean,sim_se,z_score,within_3se")
	var sl1, sl12, sl2, sl22, su1, su12, su2, su22 float64
	for r := uint32(0); r < 50; r++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(600, r)
		rng := newXoshiro(splitmix64(&x))
		l1, l2, u1, u2 := runTandemSim(&rng, 6, 6, 800, 1400, 1200, valSlots)
		sl1 += l1
		sl12 += l1 * l1
		sl2 += l2
		sl22 += l2 * l2
		su1 += u1
		su12 += u1 * u1
		su2 += u2
		su22 += u2 * u2
	}
	for _, item := range []struct {
		nm         string
		dt, s1, s2 float64
	}{{"L1", t1.L, sl1, sl12}, {"L2", t2.L, sl2, sl22}, {"U1", t1.U, su1, su12}, {"U2", t2.U, su2, su22}} {
		mn := item.s1 / 50.0
		se := math.Sqrt(math.Max((item.s2/50.0-mn*mn)/49.0, 1e-18))
		z := (mn - item.dt) / se
		fmt.Fprintf(w, "%s,%.6f,%.6f,%.6f,%.3f,%d\n", item.nm, item.dt, mn, se, z, boolInt(math.Abs(z) <= 3.0))
	}
	closeBW(f, w)

	rpw := defaultRP(c, 2)
	rpw.auditOn = 1
	sw := newSim(c, rpw)
	simRun(sw)
	nBlk := int(sw.auditSlots) / 5
	blkY := make([]float64, nBlk)
	for b := 0; b < nBlk; b++ {
		sm := 0.0
		for k := 0; k < 5; k++ {
			sm += float64(sw.qSeries[(b*5+k)*NStn+StnSecurity])
		}
		blkY[b] = sm / 5.0
	}
	bestD, bestMser := 0, 1e30
	for d := 0; d < nBlk/2; d += 10 {
		m := float64(nBlk - d)
		mn, s2 := 0.0, 0.0
		for i := d; i < nBlk; i++ {
			mn += blkY[i]
		}
		mn /= m
		for i := d; i < nBlk; i++ {
			df := blkY[i] - mn
			s2 += df * df
		}
		mser := s2 / (m * m)
		if mser < bestMser {
			bestMser = mser
			bestD = d * 5
		}
	}
	f, w = openBW(outDir, "s12_warmup.csv")
	fmt.Fprintf(w, "warmup_slots,mser_min\n%d,%.6f\n", bestD, bestMser)
	closeBW(f, w)

	f, w = openBW(outDir, "s13_calibration.csv")
	fmt.Fprintln(w, "arr_scale_bp,max_util,loss_frac,calibrated")
	for _, scale := range []int64{15000, 18000, 22000, 26000, 30000} {
		rpc := defaultRP(c, 2)
		rpc.arrScaleBp = scale
		sc := newSim(c, rpc)
		simRun(sc)
		maxU := 0.0
		for st := 0; st < 4; st++ {
			u := float64(sc.stn[st].busySlots) / float64(2*SlotsPerDay*int64(sc.stn[st].cMax))
			if u > maxU {
				maxU = u
			}
		}
		lf := 0.0
		if sc.totalArr > 0 {
			lf = float64(sc.totalLost) / float64(sc.totalArr)
		}
		fmt.Fprintf(w, "%d,%.4f,%.4f,%d\n", scale, maxU, lf, boolInt(maxU >= 0.55 && maxU <= 0.92 && lf < 0.15))
	}
	closeBW(f, w)
}
