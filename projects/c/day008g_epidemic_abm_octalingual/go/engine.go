package main

import "math"

func analyticalR0(c Config, betaBp int64, cHH, cSW, cComm float64) float64 {
	beta := float64(betaBp) / float64(BP)
	dInf := 0.5 * float64(c.infDaysMin+c.infDaysMax)
	pIA := float64(c.asympProbBp) / float64(BP)
	rIA := float64(c.asympRelBp) / float64(BP)
	return beta * (cHH + cSW + cComm) * dInf * ((1.0 - pIA) + pIA*rIA)
}

func solveFinalSize(r0 float64) float64 {
	if r0 <= 1.0 {
		return 0.0
	}
	z := 0.5
	for i := 0; i < 80; i++ {
		z = 1.0 - math.Exp(-r0*z)
	}
	return z
}

func checkInvariants(s *Sim, day int) bool {
	var cnt [8]int64
	for i := range s.pop {
		cnt[s.pop[i].status]++
	}
	sum := int64(0)
	for k := 0; k < 8; k++ {
		sum += cnt[k]
		s.dailyCnt[day][k] = cnt[k]
	}
	ok := sum == s.n
	if day > 0 && s.dailyCnt[day][StD] < s.dailyCnt[day-1][StD] {
		ok = false
	}
	if s.hospOcc > s.hospCap || s.icuOcc > s.icuCap {
		ok = false
	}
	if s.logEnabled && int64(len(s.trans)) != s.totalTrans {
		ok = false
	}
	if !ok {
		s.invViolations++
	}
	return ok
}

func simRun(s *Sim) {
	n := int(s.n)
	nSeed := s.cfg.initInfected
	if s.rp.betaBp == 0 && s.rp.contactScaleBp < 0 {
		nSeed = 0
	}
	if s.rp.vaxPerfect != 0 {
		for i := 0; i < n; i++ {
			s.pop[i].vaxDoses = 2
			s.pop[i].vaxDay = 0
		}
	}
	if nSeed > 0 {
		seeded := int64(0)
		for k := 0; k < n && seeded < nSeed; k++ {
			pid := int(s.rng[2].randInt(0, s.n-1))
			if s.pop[pid].status == StS {
				s.pop[pid].everInf = 1
				s.pop[pid].infDay = 0
				simSetState(s, pid, StIS, 0, 4)
				seeded++
			}
		}
	}
	for d := 0; d < int(s.rp.days); d++ {
		di := int32(d)
		simIntervStep(s, di)
		simContactStep(s, di)
		simDiseaseStep(s, di)
		s.dailyHosp[d][0] = s.hospOcc
		s.dailyHosp[d][1] = s.icuOcc
		s.dailyHosp[d][2] = s.branches[9]
		s.dailyHosp[d][3] = s.branches[11]
		q := int64(0)
		for i := 0; i < n; i++ {
			if s.pop[i].quarUntil >= di {
				q++
			}
		}
		s.dailyQuar[d] = q
		checkInvariants(s, d)
	}
}
