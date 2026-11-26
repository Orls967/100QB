package main

func vaxEffBp(s *Sim, p *Person, day int32) int64 {
	if p.vaxDoses <= 0 || p.vaxDay < 0 {
		return 0
	}
	if s.rp.vaxPerfect != 0 {
		return BP
	}
	base := s.cfg.vaxEff2Bp
	if p.vaxDoses == 1 {
		base = s.cfg.vaxEff1Bp
	}
	elapsed := int64(day - p.vaxDay)
	if elapsed <= 30 {
		return base
	}
	return max64(base-(elapsed-30)*35, 2000)
}

func tryInfect(s *Sim, src, tgt, lay int, day int32, betaEff int64) {
	if s.pop[tgt].status != StS || s.seenDay[tgt] == day {
		return
	}
	if s.pop[tgt].quarUntil >= day {
		s.branches[5]++
		return
	}
	pBp := betaEff
	if s.rp.intervOn != 0 {
		mDay := s.cfg.maskDay
		if s.rp.maskOvr >= 0 {
			mDay = s.rp.maskOvr
		}
		if int64(day) >= mDay {
			red := (s.cfg.maskEffBp * int64(s.pop[src].compliance)) / BP
			nextBp := (pBp * (BP - red)) / BP
			if nextBp < pBp {
				s.branches[4]++
			}
			pBp = nextBp
		}
	}
	ve := vaxEffBp(s, &s.pop[tgt], day)
	if ve > 0 {
		pBp = (pBp * (BP - ve)) / BP
		if pBp <= 0 {
			s.branches[3]++
			return
		}
	}
	if s.rng[1].randInt(1, BP) <= pBp {
		s.seenDay[tgt] = day
		s.candTgt = append(s.candTgt, tgt)
		s.candSrc = append(s.candSrc, src)
		s.candLay = append(s.candLay, lay)
		s.candVar = append(s.candVar, s.pop[src].variant)
	} else if ve > 0 {
		s.branches[3]++
	}
}

func simContactStep(s *Sim, day int32) {
	s.candTgt = s.candTgt[:0]
	s.candSrc = s.candSrc[:0]
	s.candLay = s.candLay[:0]
	s.candVar = s.candVar[:0]
	if s.rp.contactScaleBp == 0 {
		return
	}
	n := int(s.n)
	scDay := s.cfg.schoolCloseDay
	if s.rp.schoolOvr >= 0 {
		scDay = s.rp.schoolOvr
	}
	lkDay := s.cfg.lockdownDay
	if s.rp.lockOvr >= 0 {
		lkDay = s.rp.lockOvr
	}
	schoolClosed := s.rp.intervOn != 0 && int64(day) >= scDay
	locked := s.rp.intervOn != 0 && int64(day) >= lkDay

	for i := 0; i < n; i++ {
		st := s.pop[i].status
		if st != StIA && st != StIS {
			continue
		}
		baseBp := s.rp.betaBp
		if st == StIA {
			baseBp = (baseBp * s.cfg.asympRelBp) / BP
		}
		if s.pop[i].variant == 1 {
			baseBp = (baseBp * s.cfg.variantMultBp) / BP
		}
		if baseBp > BP {
			baseBp = BP
		}

		h := int(s.pop[i].hhID)
		for idx := s.hhOff[h]; idx < s.hhOff[h+1]; idx++ {
			j := s.hhMem[idx]
			if j != i {
				tryInfect(s, i, j, 0, day, baseBp)
			}
		}
		if s.pop[i].quarUntil >= day {
			s.branches[5]++
			continue
		}

		if s.pop[i].swID >= 0 && !(s.pop[i].ageGrp == 0 && schoolClosed) {
			sw := int(s.pop[i].swID)
			lo, hi := s.swOff[sw], s.swOff[sw+1]
			sz := hi - lo
			if sz > 1 {
				nSwC := 5
				if s.pop[i].ageGrp == 0 {
					nSwC = 6
				}
				if locked {
					nSwC = (nSwC + 1) / 2
				}
				for c := 0; c < nSwC; c++ {
					j := s.swMem[lo+int(s.rng[1].randInt(0, int64(sz)-1))]
					if j != i {
						tryInfect(s, i, j, 1, day, baseBp)
					}
				}
			}
		}
		nComm := 3
		if s.pop[i].superFlag != 0 {
			nComm = 12
		}
		if locked {
			nComm = (nComm + 1) / 3
		}
		for c := 0; c < nComm; c++ {
			j := int(s.rng[1].randInt(0, s.n-1))
			if j != i {
				tryInfect(s, i, j, 2, day, baseBp)
			}
		}
	}
}
