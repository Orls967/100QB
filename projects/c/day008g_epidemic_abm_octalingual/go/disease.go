package main

func simSetState(s *Sim, pid int, newSt, day, dur int32) {
	oldSt := s.pop[pid].status
	s.pop[pid].status = newSt
	s.pop[pid].stateDay = day
	s.pop[pid].durDays = dur
	if s.logEnabled {
		s.progs = append(s.progs, ProgRec{day: day, pid: pid, fromSt: oldSt, toSt: newSt})
	}
}

func simDiseaseStep(s *Sim, day int32) {
	n := int(s.n)
	for i := 0; i < n; i++ {
		elapsed := day - s.pop[i].stateDay
		st := s.pop[i].status
		if st == StE && elapsed >= s.pop[i].durDays {
			dur := int32(s.rng[2].randInt(s.cfg.infDaysMin, s.cfg.infDaysMax))
			if s.rng[2].randInt(1, BP) <= s.cfg.asympProbBp {
				s.branches[6]++
				simSetState(s, i, StIA, day, dur)
			} else {
				s.branches[7]++
				simSetState(s, i, StIS, day, dur)
			}
		} else if (st == StIA || st == StIS) && elapsed >= s.pop[i].durDays {
			if st == StIS {
				pHosp := int64(1800)
				if s.pop[i].ageGrp == 0 {
					pHosp = 120
				} else if s.pop[i].ageGrp == 1 {
					pHosp = 450
				}
				if s.pop[i].vaxDoses > 0 {
					pHosp = (pHosp * (BP - s.cfg.vaxHospBp)) / BP
				}
				if s.rng[4].randInt(1, BP) <= pHosp {
					if s.hospOcc < s.hospCap {
						s.hospOcc++
						s.pop[i].hospRej = 0
						s.branches[8]++
					} else {
						s.pop[i].hospRej = 1
						s.branches[9]++
					}
					dur := int32(s.rng[4].randInt(4, 8))
					simSetState(s, i, StH, day, dur)
					continue
				}
			}
			s.branches[13]++
			simSetState(s, i, StR, day, int32(s.cfg.waningMinDays))
		} else if st == StH && elapsed >= s.pop[i].durDays {
			if s.pop[i].hospRej == 0 && s.hospOcc > 0 {
				s.hospOcc--
			}
			pICU := int64(1800)
			if s.pop[i].ageGrp == 2 {
				pICU = 3500
			}
			if s.rng[4].randInt(1, BP) <= pICU {
				if s.icuOcc < s.icuCap {
					s.icuOcc++
					s.pop[i].hospRej = 0
					s.branches[10]++
				} else {
					s.pop[i].hospRej = 1
					s.branches[11]++
				}
				dur := int32(s.rng[4].randInt(4, 7))
				simSetState(s, i, StICU, day, dur)
			} else {
				pDie := int64(2000)
				if s.pop[i].ageGrp == 0 {
					pDie = 100
				} else if s.pop[i].ageGrp == 1 {
					pDie = 500
				}
				if s.pop[i].hospRej != 0 {
					pDie *= 2
				}
				if s.rng[4].randInt(1, BP) <= pDie {
					s.branches[12]++
					s.dailyDeathsAge[day][s.pop[i].ageGrp]++
					simSetState(s, i, StD, day, 0)
				} else {
					s.branches[13]++
					simSetState(s, i, StR, day, int32(s.cfg.waningMinDays))
				}
			}
		} else if st == StICU && elapsed >= s.pop[i].durDays {
			if s.pop[i].hospRej == 0 && s.icuOcc > 0 {
				s.icuOcc--
			}
			pDie := int64(4800)
			if s.pop[i].ageGrp == 0 {
				pDie = 1500
			} else if s.pop[i].ageGrp == 1 {
				pDie = 2800
			}
			if s.pop[i].hospRej != 0 {
				pDie *= 2
			}
			if pDie > BP {
				pDie = BP
			}
			if s.rng[4].randInt(1, BP) <= pDie {
				s.branches[12]++
				s.dailyDeathsAge[day][s.pop[i].ageGrp]++
				simSetState(s, i, StD, day, 0)
			} else {
				s.branches[13]++
				simSetState(s, i, StR, day, int32(s.cfg.waningMinDays))
			}
		} else if st == StR && elapsed >= s.pop[i].durDays {
			if s.rng[2].randInt(1, BP) <= s.cfg.waningBp {
				s.branches[14]++
				simSetState(s, i, StS, day, 0)
			}
		}
	}

	for c := 0; c < len(s.candTgt); c++ {
		tgt, src, lay, vari := s.candTgt[c], s.candSrc[c], s.candLay[c], s.candVar[c]
		if s.pop[tgt].status != StS {
			continue
		}
		if s.pop[tgt].everInf != 0 {
			s.totalReinf++
		}
		s.pop[tgt].everInf = 1
		s.pop[tgt].variant = vari
		s.pop[tgt].infDay = day
		s.pop[src].offspring++
		s.layerInf[lay]++
		s.branches[lay]++
		ag := s.pop[tgt].ageGrp
		s.dailyIncAge[day][ag]++
		s.ageInf[ag]++
		s.dailyVar[day][vari]++
		s.totalTrans++
		gi := day - s.pop[src].infDay
		if gi < 1 {
			gi = 1
		}
		if s.logEnabled {
			s.trans = append(s.trans, TransRec{day: day, src: src, dst: tgt, lay: lay, vari: vari, genInt: gi})
		}
		s.eventHash = fnv1aI64(s.eventHash, int64(day)*1000000+int64(tgt))
		simSetState(s, tgt, StE, day, int32(s.cfg.latentDays))
	}
}
