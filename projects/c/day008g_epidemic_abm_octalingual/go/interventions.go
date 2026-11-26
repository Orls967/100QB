package main

func simIntervStep(s *Sim, day int32) {
	n := int(s.n)
	vd := s.cfg.variantDay
	if s.rp.varDayOvr >= 0 {
		vd = s.rp.varDayOvr
	}
	if int64(day) == vd {
		imported := 0
		for k := 0; k < n && imported < 5; k++ {
			pid := int(s.rng[3].randInt(0, s.n-1))
			if s.pop[pid].status == StS {
				s.pop[pid].variant = 1
				s.pop[pid].everInf = 1
				s.pop[pid].infDay = day
				simSetState(s, pid, StIS, day, 4)
				s.branches[20]++
				imported++
			}
		}
	}
	if s.rp.intervOn == 0 {
		return
	}

	if int64(day) >= s.cfg.testDay {
		tp := s.cfg.testProbBp
		if s.rp.testOvr >= 0 {
			tp = s.rp.testOvr
		}
		cp := s.rp.complianceOvr
		for i := 0; i < n; i++ {
			if s.pop[i].status == StIS && s.pop[i].quarUntil < day {
				effComp := int64(s.pop[i].compliance)
				if cp >= 0 {
					effComp = cp
				}
				if s.rng[3].randInt(1, BP) <= (tp*effComp)/BP {
					s.dailyTests[day][0]++
					if s.rng[3].randInt(1, BP) <= 8500 {
						s.dailyTests[day][1]++
						s.branches[15]++
						s.pop[i].quarUntil = day + 7
						h := int(s.pop[i].hhID)
						for idx := s.hhOff[h]; idx < s.hhOff[h+1]; idx++ {
							j := s.hhMem[idx]
							if j != i && s.pop[j].quarUntil < day && s.rng[3].randInt(1, BP) <= s.cfg.traceProbBp {
								s.pop[j].quarUntil = day + 5
								s.branches[17]++
							}
						}
					} else {
						s.branches[16]++
					}
				}
			}
		}
	}

	vday := s.cfg.vaxDay
	if s.rp.vaxDayOvr >= 0 {
		vday = s.rp.vaxDayOvr
	}
	if int64(day) >= vday {
		covCap := (s.n * 8000) / BP
		if s.rp.vaxCovOvr >= 0 {
			covCap = (s.n * s.rp.vaxCovOvr) / BP
		}
		dosesToday := max64((s.n*s.cfg.vaxDailyBp)/BP, 1)
		totalVaxed := int64(0)
		for i := 0; i < n; i++ {
			if s.pop[i].vaxDoses > 0 {
				totalVaxed++
			}
		}
		if totalVaxed >= covCap {
			s.branches[19]++
			return
		}
		given := int64(0)
		for grp := int32(2); grp >= 0 && given < dosesToday; grp-- {
			for i := 0; i < n && given < dosesToday; i++ {
				if s.pop[i].ageGrp != grp || s.pop[i].status == StD {
					continue
				}
				if s.pop[i].vaxDoses == 0 && totalVaxed < covCap {
					s.pop[i].vaxDoses = 1
					s.pop[i].vaxDay = day
					given++
					totalVaxed++
					s.dailyVax[day][0]++
					s.branches[18]++
				} else if s.pop[i].vaxDoses == 1 && day-s.pop[i].vaxDay >= 21 {
					s.pop[i].vaxDoses = 2
					s.pop[i].vaxDay = day
					given++
					s.dailyVax[day][1]++
					s.branches[18]++
				}
			}
		}
		if given == 0 {
			s.branches[19]++
		}
	}
}
