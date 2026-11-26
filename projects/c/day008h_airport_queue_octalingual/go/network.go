package main

const (
	BrArrAcc       = 0
	BrArrRej       = 1
	BrSvcDone      = 2
	BrSecExtra     = 3
	BrSrvBreak     = 4
	BrSrvRepair    = 5
	BrAbandon      = 6
	BrMissedFlt    = 7
	BrRebook       = 8
	BrBagMiss      = 9
	BrBeltBlock    = 10
	BrRosterDrop   = 11
	BrRouteIntl    = 12
	BrRouteTransit = 13
	BrBizPrio      = 14
)

type StnEvent struct {
	slot, stn, delta int32
}

type Sim struct {
	cfg             Config
	rp              RunParams
	rng             [8]Xoshiro256
	pool            [MaxPool]Agent
	freeStk         [MaxPool]int32
	freeTop         int
	stn             [NStn]Station
	arrProf         [24]int32
	roster          [24][NStn]int32
	branches        [15]int64
	totalArr        int64
	totalDone       int64
	totalLost       int64
	totalMissed     int64
	totalRebook     int64
	totalAban       int64
	totalBags       int64
	totalBagsMiss   int64
	totalBeltBlock  int64
	totalBreakdowns int64
	sysSumL         int64
	sysSumSojourn   int64
	invViolations   int64
	eventHash       uint64
	clsWaitSum      [3]int64
	clsWaitCnt      [3]int64
	clsSojournSum   [3]int64
	clsSojournCnt   [3]int64
	stnWaitSum      [NStn]int64
	stnWaitCnt      [NStn]int64
	hrArr           [24]int64
	hrUtilSum       [24][NStn]int64
	waitSamples     []int32
	sojournSamples  []int32
	qSeries         []int32
	srvSeries       []int32
	stnEv           [NStn][]StnEvent
	auditSlots      int32
}

func newSim(c Config, rp RunParams) *Sim {
	s := &Sim{cfg: c, rp: rp, eventHash: 0xcbf29ce484222325}
	for i := uint32(0); i < 8; i++ {
		x := uint64(c.masterSeed) ^ fnv1aU32Pair(200, i)
		s.rng[i] = newXoshiro(splitmix64(&x))
	}
	for i := MaxPool - 1; i >= 0; i-- {
		s.freeStk[s.freeTop] = int32(i)
		s.freeTop++
	}
	caps := [NStn]int32{1000, 1000, 1000, 1000, int32(rp.beltCap), 1000}
	prios := [NStn]int32{1, 1, 0, 0, 0, 0}
	for i := 0; i < NStn; i++ {
		cm := c.srvCnt[i]
		if rp.srvOvr[i] > 0 {
			cm = rp.srvOvr[i]
		}
		if cm > int64(MaxSrv) {
			cm = int64(MaxSrv)
		}
		initStation(&s.stn[i], int32(cm), caps[i], prios[i])
	}
	prof := [24]int32{180, 140, 120, 150, 260, 420, 680, 860, 920, 820, 700, 660, 680, 740, 820, 880, 940, 900, 760, 620, 480, 360, 280, 210}
	s.arrProf = prof
	for h := 0; h < 24; h++ {
		for st := 0; st < NStn; st++ {
			cm := s.stn[st].cMax
			if h >= 6 && h <= 20 {
				s.roster[h][st] = cm
			} else if cm > 1 {
				s.roster[h][st] = cm - 1
			} else {
				s.roster[h][st] = 1
			}
		}
	}
	if rp.auditOn != 0 {
		d := rp.days
		if AuditDays < d {
			d = AuditDays
		}
		s.auditSlots = int32(d * SlotsPerDay)
		s.qSeries = make([]int32, int(s.auditSlots)*NStn)
		s.srvSeries = make([]int32, int(s.auditSlots)*3)
	}
	return s
}

func (s *Sim) effBp(base int64) int64 {
	if s.rp.timescaleDiv2 != 0 {
		return base / 2
	}
	return base
}

func (s *Sim) enqueueAgent(stnID int, aidx int32, slot int64, toHead bool) {
	st := &s.stn[stnID]
	if st.len+st.busyCnt >= st.cap {
		s.branches[BrArrRej]++
		st.blockCnt++
		if stnID == StnBelt {
			s.totalBeltBlock++
		}
		if s.pool[aidx].isBag == 0 {
			s.totalLost++
			s.sysSumSojourn += slot - s.pool[aidx].arrSlot
		}
		poolFree(s, aidx)
		return
	}
	s.pool[aidx].curStn = int32(stnID)
	s.pool[aidx].qEnterSlot = slot
	st.arrCnt++
	if toHead {
		st.pushHead(aidx)
	} else {
		st.pushTail(aidx)
	}
	if s.rp.auditOn != 0 && int32(slot) < s.auditSlots {
		s.stnEv[stnID] = append(s.stnEv[stnID], StnEvent{slot: int32(slot), stn: int32(stnID), delta: 1})
	}
}

func simStepSlot(s *Sim, slot int64) {
	hr := int((slot % SlotsPerDay) / SlotsPerHour)
	var compAidx [NStn * MaxSrv]int32
	var compStn [NStn * MaxSrv]int
	nComp := 0

	arrP := (int64(s.arrProf[hr]) * s.rp.arrScaleBp) / BP
	if s.rp.timescaleDiv2 != 0 {
		arrP /= 2
	}
	if s.rng[0].bernoulli(arrP) != 0 {
		aidx := poolAlloc(s)
		if aidx >= 0 {
			ag := &s.pool[aidx]
			uCls := s.rng[0].randInt(1, BP)
			if uCls <= s.cfg.bizProbBp {
				ag.cls = ClsBiz
			} else if uCls <= s.cfg.bizProbBp+s.cfg.transitProbBp {
				ag.cls = ClsTransit
			} else {
				ag.cls = ClsEcon
			}
			if s.rng[0].bernoulli(s.cfg.intlProbBp) != 0 {
				ag.dest = DestIntl
			} else {
				ag.dest = DestDom
			}
			ag.nBags = int32(s.rng[0].randInt(0, 2))
			ag.arrSlot = slot
			ag.flightHr = int32(hr)
			ag.depSlot = slot + s.rng[0].randInt(110, 220)
			s.totalArr++
			s.hrArr[hr]++
			s.branches[BrArrAcc]++
			s.eventHash = fnv1aI64(s.eventHash, slot*100+int64(ag.cls)*10+int64(ag.dest))
			startStn := StnCheckin
			if ag.cls == ClsTransit {
				startStn = StnSecurity
				s.branches[BrRouteTransit]++
			}
			if ag.nBags > 0 && ag.cls != ClsTransit {
				for b := int32(0); b < ag.nBags; b++ {
					bidx := poolAlloc(s)
					if bidx >= 0 {
						s.pool[bidx].isBag = 1
						s.pool[bidx].cls = ag.cls
						s.pool[bidx].dest = ag.dest
						s.pool[bidx].arrSlot = slot
						s.pool[bidx].depSlot = ag.depSlot
						s.totalBags++
						s.enqueueAgent(StnBelt, bidx, slot, false)
					}
				}
			}
			s.enqueueAgent(startStn, aidx, slot, false)
		}
	}

	for st := 0; st < NStn; st++ {
		sn := &s.stn[st]
		baseMu := s.cfg.muBp[st]
		if s.rp.muOvr[st] > 0 {
			baseMu = s.rp.muOvr[st]
		}
		for idxK := int32(0); idxK < sn.cMax; idxK++ {
			k := idxK
			if s.rp.permuteSrv != 0 {
				k = (idxK + 1) % sn.cMax
			}
			sv := &sn.srv[k]
			if sv.status == SrvBusy && sv.agentIdx >= 0 {
				mu := baseMu
				if sv.extraSec != 0 {
					mu = s.cfg.securityExtraMuBp
				}
				if s.rng[1].bernoulli(s.effBp(mu)) != 0 {
					compAidx[nComp] = sv.agentIdx
					compStn[nComp] = st
					nComp++
					sv.agentIdx = -1
					sv.status = SrvIdle
					sv.extraSec = 0
					sn.busyCnt--
					sn.depCnt++
					s.branches[BrSvcDone]++
					if s.rp.auditOn != 0 && int32(slot) < s.auditSlots {
						s.stnEv[st] = append(s.stnEv[st], StnEvent{slot: int32(slot), stn: int32(st), delta: -1})
					}
				}
			}
		}
	}

	for i := 0; i < nComp; i++ {
		aidx := compAidx[i]
		st := compStn[i]
		ag := &s.pool[aidx]
		if ag.isBag != 0 {
			if s.rng[4].bernoulli(s.cfg.bagMissBp) != 0 {
				s.totalBagsMiss++
				s.branches[BrBagMiss]++
			}
			poolFree(s, aidx)
			continue
		}
		nextSt := -1
		if st == StnCheckin {
			nextSt = StnSecurity
		} else if st == StnSecurity {
			if ag.dest == DestIntl {
				nextSt = StnImigrasi
				s.branches[BrRouteIntl]++
			} else {
				nextSt = StnGate
			}
		} else if st == StnImigrasi {
			nextSt = StnGate
		} else if st == StnGate {
			lead := s.rp.gateCloseLead
			if s.rp.timescaleDiv2 != 0 {
				lead *= 2
			}
			if slot >= ag.depSlot-lead {
				s.totalMissed++
				s.branches[BrMissedFlt]++
				if s.rng[2].bernoulli(s.cfg.rebookProbBp) != 0 {
					s.totalRebook++
					s.branches[BrRebook]++
				}
				s.totalLost++
				s.sysSumSojourn += slot - ag.arrSlot
				poolFree(s, aidx)
				continue
			}
			if s.rp.customsOn != 0 && ag.dest == DestIntl {
				nextSt = StnCustoms
			}
		}
		if nextSt >= 0 {
			s.enqueueAgent(nextSt, aidx, slot, false)
		} else {
			s.totalDone++
			soj := slot - ag.arrSlot
			s.sysSumSojourn += soj
			s.clsSojournSum[ag.cls] += soj
			s.clsSojournCnt[ag.cls]++
			if len(s.sojournSamples) < 50000 {
				s.sojournSamples = append(s.sojournSamples, int32(soj))
			}
			poolFree(s, aidx)
		}
	}

	for st := 0; st < NStn; st++ {
		sn := &s.stn[st]
		for idxK := int32(0); idxK < sn.cMax; idxK++ {
			k := idxK
			if s.rp.permuteSrv != 0 {
				k = (idxK + 1) % sn.cMax
			}
			sv := &sn.srv[k]
			if sv.enabled != 0 && sv.status == SrvIdle && sn.len > 0 {
				aidx := sn.popNext(s.pool[:], s.rp.priorityOn, &s.branches[BrBizPrio])
				if aidx >= 0 {
					w := slot - s.pool[aidx].qEnterSlot
					sn.sumW += w
					sn.waitCnt++
					s.stnWaitSum[st] += w
					s.stnWaitCnt[st]++
					if s.pool[aidx].isBag == 0 {
						s.clsWaitSum[s.pool[aidx].cls] += w
						s.clsWaitCnt[s.pool[aidx].cls]++
						if len(s.waitSamples) < 50000 {
							s.waitSamples = append(s.waitSamples, int32(w))
						}
					}
					sv.status = SrvBusy
					sv.agentIdx = aidx
					sv.extraSec = 0
					if st == StnSecurity && s.rng[5].bernoulli(s.cfg.securityExtraProbBp) != 0 {
						sv.extraSec = 1
						s.branches[BrSecExtra]++
					}
					sn.busyCnt++
				}
			}
		}
	}

	pat := s.cfg.patienceSlots
	if s.rp.timescaleDiv2 != 0 {
		pat *= 2
	}
	for st := 0; st < NStn; st++ {
		if st == StnBelt {
			continue
		}
		sn := &s.stn[st]
		for pos := int32(0); pos < sn.len; {
			aidx := sn.q[(sn.head+pos)%int32(MaxQ)]
			if slot-s.pool[aidx].qEnterSlot >= pat {
				sn.removeAt(pos)
				sn.abanCnt++
				sn.depCnt++
				s.totalAban++
				s.totalLost++
				s.branches[BrAbandon]++
				s.sysSumSojourn += slot - s.pool[aidx].arrSlot
				if s.rp.auditOn != 0 && int32(slot) < s.auditSlots {
					s.stnEv[st] = append(s.stnEv[st], StnEvent{slot: int32(slot), stn: int32(st), delta: -1})
				}
				poolFree(s, aidx)
			} else {
				pos++
			}
		}
	}

	if s.rp.breakdownOn != 0 && s.rp.breakProbBp > 0 {
		bp := s.effBp(s.rp.breakProbBp)
		rpB := s.effBp(s.rp.repairProbBp)
		for st := 0; st < NStn; st++ {
			sn := &s.stn[st]
			for idxK := int32(0); idxK < sn.cMax; idxK++ {
				k := idxK
				if s.rp.permuteSrv != 0 {
					k = (idxK + 1) % sn.cMax
				}
				sv := &sn.srv[k]
				if sv.status == SrvDown {
					if s.rng[3].bernoulli(rpB) != 0 {
						sv.status = SrvIdle
						s.branches[BrSrvRepair]++
					}
				} else if sv.enabled != 0 && s.rng[3].bernoulli(bp) != 0 {
					if sv.status == SrvBusy && sv.agentIdx >= 0 {
						sn.pushHead(sv.agentIdx)
						sv.agentIdx = -1
						sn.busyCnt--
					}
					sv.status = SrvDown
					s.totalBreakdowns++
					s.branches[BrSrvBreak]++
				}
			}
		}
	}

	if s.rp.rosterOn != 0 && (slot%SlotsPerHour) == 0 {
		for st := 0; st < NStn; st++ {
			sn := &s.stn[st]
			target := s.roster[hr][st]
			if target < sn.cActive {
				s.branches[BrRosterDrop]++
			}
			sn.cActive = target
			for idxK := int32(0); idxK < sn.cMax; idxK++ {
				k := idxK
				if s.rp.permuteSrv != 0 {
					k = (idxK + 1) % sn.cMax
				}
				if idxK < target {
					sn.srv[k].enabled = 1
				} else {
					sn.srv[k].enabled = 0
				}
			}
		}
	}

	var ti, tb, td int32
	for st := 0; st < NStn; st++ {
		sn := &s.stn[st]
		inSt := sn.len + sn.busyCnt
		sn.sumL += int64(inSt)
		sn.busySlots += int64(sn.busyCnt)
		if st != StnBelt {
			s.sysSumL += int64(inSt)
		}
		s.hrUtilSum[hr][st] += int64(sn.busyCnt)
		if sn.busyCnt > sn.cMax || inSt > sn.cap {
			s.invViolations++
		}
		if s.rp.auditOn != 0 && int32(slot) < s.auditSlots {
			s.qSeries[int(slot)*NStn+st] = inSt
		}
		for k := int32(0); k < sn.cMax; k++ {
			if sn.srv[k].status == SrvIdle && sn.srv[k].enabled != 0 {
				ti++
			} else if sn.srv[k].status == SrvBusy {
				tb++
			} else if sn.srv[k].status == SrvDown {
				td++
			}
		}
	}
	if s.rp.auditOn != 0 && int32(slot) < s.auditSlots {
		s.srvSeries[int(slot)*3] = ti
		s.srvSeries[int(slot)*3+1] = tb
		s.srvSeries[int(slot)*3+2] = td
	}
}
