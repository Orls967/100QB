package main

import "fmt"

func doShip(s *Sim, src, dst, sku int, qty, day, tick, tAbs int64, bin int) {
	if qty <= 0 || s.acc[src].onHand[sku] <= 0 {
		return
	}
	if qty > s.acc[src].onHand[sku] {
		qty = s.acc[src].onHand[sku]
	}
	whCap := int64(20000) * s.capScale
	cCap := int64(12000) * s.capScale
	if dst >= 10 && dst < 14 {
		tot := int64(0)
		for k := 0; k < 4; k++ {
			tot += s.acc[dst].onHand[k]
		}
		if tot+qty > whCap {
			s.branches[4]++
			return
		}
	}
	cid := (src + dst + sku) % numCarriers
	if s.carrierDayVol[cid]+qty > cCap {
		cid = (cid + 1) % numCarriers
		if s.carrierDayVol[cid]+qty > cCap {
			s.branches[5]++
			return
		}
	}
	rate := carrierRate(s.carrierDayVol[cid])
	tier := 2
	if rate == 200 {
		tier = 0
	} else if rate == 160 {
		tier = 1
	}
	s.carrierTierCnt[tier]++
	s.carrierDayVol[cid] += qty
	s.carrierVolTot += qty
	if bin >= 0 && bin < 100 {
		s.binOrders[sku][bin] += qty
	}
	goods, freight := qty*skuWholesale[sku], qty*rate
	s.acc[src].onHand[sku] -= qty
	s.acc[dst].inTransit[sku] += qty
	transferCash(s, dst, src, goods, -1, day, tick)
	transferCash(s, dst, 17+cid, freight, -1, day, tick)
	if day < 6 {
		s.dailyCostEchelon[day][2] += freight
	}
	delay := s.rng[17+cid].randInt(0, 25)
	if delay > 18 {
		s.branches[15]++
	}
	lt := carrierBaseLt[cid] + delay
	if len(s.leadtimePairs) < 256 {
		s.leadtimePairs = append(s.leadtimePairs, [2]int64{lt, qty})
	}
	if s.logEnabled {
		s.ships = append(s.ships, ShipRec{day, tick, int32(src), int32(dst), int32(cid), int32(sku), int32(qty), 1, goods, freight})
	}
	heapPush(&s.heap, SimEvent{tAbs + lt, s.nextSeq, 2, int32(src), int32(dst), int32(sku), int32(qty), int32(cid)})
	s.nextSeq++
}
func hashEv(s *Sim, ev *SimEvent) {
	s.eventHash = fnv1aI64(s.eventHash, ev.time)
	s.eventHash = fnv1aI64(s.eventHash, int64(ev.etype))
	s.eventHash = fnv1aI64(s.eventHash, int64(ev.src))
	s.eventHash = fnv1aI64(s.eventHash, int64(ev.dst))
	s.eventHash = fnv1aI64(s.eventHash, int64(ev.sku))
	s.eventHash = fnv1aI64(s.eventHash, int64(ev.qty))
}
func simRunDays(s *Sim, days, evPerDay, safetyOvr, disruptOvr int64, transshipOn int) {
	T, evProc, cs := s.cfg.ticksPerDay, int64(0), s.capScale
	chkStep := max(evPerDay/8, 1)
	shkStep := max(evPerDay/250, 1)
	for d := int64(0); d < days; d++ {
		dayBase := d * T
		s.carrierCapTot += numCarriers * 12000 * cs
		s.carrierDayVol = [3]int64{}
		for f := 14; f < 17; f++ {
			if s.acc[f].plantState == 1 {
				s.acc[f].disruptRem--
				if s.acc[f].disruptRem <= 0 {
					s.acc[f].plantState = 2
				}
			} else if s.acc[f].plantState == 2 {
				s.acc[f].plantState = 0
			}
		}
		if d > 0 {
			for wh := 10; wh < 13; wh++ {
				for sku := 0; sku < 4; sku++ {
					var sts []int
					var reqs, seqs []int64
					var prio []int32
					sLvl := int64(450)
					if safetyOvr > 0 {
						sLvl = safetyOvr
					}
					sLvl *= cs
					for st := 0; st < numStores; st++ {
						if storeWh[st] != wh {
							continue
						}
						inv := s.acc[st].onHand[sku] + s.acc[st].inTransit[sku]
						if inv < sLvl {
							sts = append(sts, st)
							reqs = append(reqs, sLvl-inv)
							seqs = append(seqs, int64(s.nextSeq))
							s.nextSeq++
							p := int32(1)
							if st%3 == 0 {
								p = 2
							}
							prio = append(prio, p)
						}
					}
					if len(sts) > 0 {
						fills := runAllocAuction(s, wh, sku, sts, reqs, seqs, prio)
						b0 := int((d * 100) / days)
						for i, st := range sts {
							if s.fAlloc != nil {
								fmt.Fprintf(s.fAlloc, "%d,%d,%d,%d,%d,%d\n", d, wh, sku, st, reqs[i], fills[i])
							}
							if fills[i] > 0 {
								doShip(s, wh, st, sku, fills[i], d, 0, dayBase, b0)
							}
						}
					}
				}
			}
		}
		for e := int64(0); e < evPerDay; e++ {
			tAbs := dayBase + (e*T)/evPerDay
			st := int32(s.rng[numAccounts].randInt(0, 9))
			sku := int32(s.rng[numAccounts].randInt(0, 3))
			var qty int32
			if s.rng[numAccounts].randInt(0, 9) == 0 {
				qty = int32(s.rng[numAccounts].randInt(12, 24))
			} else {
				qty = int32(s.rng[numAccounts].randInt(1, 5))
			}
			heapPush(&s.heap, SimEvent{tAbs, s.nextSeq, 0, accCust, st, sku, qty, 0})
			s.nextSeq++
			if e%20 == 0 {
				f := int32((e / 20) % 3)
				heapPush(&s.heap, SimEvent{tAbs + 1, s.nextSeq, 3, 14 + f, 10 + f, f, 130, 0})
				s.nextSeq++
			}
			if e%40 == 0 {
				kq := int32(s.rng[13].randInt(15, 35))
				heapPush(&s.heap, SimEvent{tAbs + 2, s.nextSeq, 4, 13, 13, 3, kq, 0})
				s.nextSeq++
			}
			if e%60 == 0 && transshipOn != 0 {
				tq := int32(s.rng[13].randInt(30, 50))
				heapPush(&s.heap, SimEvent{tAbs + 3, s.nextSeq, 5, 10, 13, 0, tq, 0})
				s.nextSeq++
			}
			if e%chkStep == 0 {
				dp := s.cfg.disruptProbPct
				if disruptOvr >= 0 {
					dp = disruptOvr
				}
				if s.rng[numAccounts+1].randInt(1, 100) <= dp {
					f := int32(s.rng[numAccounts+1].randInt(0, 2))
					heapPush(&s.heap, SimEvent{tAbs + 4, s.nextSeq, 7, 14 + f, 14 + f, f, 1, 0})
					s.nextSeq++
				}
				heapPush(&s.heap, SimEvent{tAbs + 6, s.nextSeq, 8, 0, 0, 0, 0, 0})
				s.nextSeq++
			}
			if e%shkStep == 0 {
				if s.rng[numAccounts+1].randInt(1, 100) <= s.cfg.shockProbPct {
					stS := int32(s.rng[numAccounts+1].randInt(0, 9))
					cls := int32(s.rng[numAccounts+1].randInt(0, 2))
					heapPush(&s.heap, SimEvent{tAbs + 5, s.nextSeq, 6, accCust, stS, 0, (cls + 1) * 10, cls})
					s.nextSeq++
				}
			}
		}
		dayEnd := (d + 1) * T
		for len(s.heap) > 0 && s.heap[0].time < dayEnd {
			ev := heapPop(&s.heap)
			tick := ev.time - dayBase
			bin := int(((d*T + tick) * 100) / (days * T))
			if bin < 0 {
				bin = 0
			} else if bin > 99 {
				bin = 99
			}
			hashEv(s, &ev)
			evProc++
			if evProc <= 5000 && s.fOrders != nil {
				fmt.Fprintf(s.fOrders, "%d,%d,%d,%d,%d,%d,%d,%d\n", d, tick, ev.seq, ev.etype, ev.src, ev.dst, ev.sku, ev.qty)
			}
			switch ev.etype {
			case 0:
				st, sku, qty := int(ev.dst), int(ev.sku), int64(ev.qty)
				s.demEventsTot++
				s.demUnitsTot += qty
				s.binDemand[sku][bin] += qty
				s.binStoreDem[st][bin] += qty
				if d < 6 {
					s.dailyDemand[d][sku] += qty
					s.dailyFillStore[st][0] += qty
				}
				s.acc[st].ewmaMilli[sku] = int32((s.cfg.ewmaAlphaMilli*qty*1000 + (1000-s.cfg.ewmaAlphaMilli)*int64(s.acc[st].ewmaMilli[sku])) / 1000)
				if s.acc[st].status == 3 {
					s.branches[2]++
				} else if s.acc[st].onHand[sku] >= qty {
					s.acc[st].onHand[sku] -= qty
					s.totalSold[sku] += qty
					s.demUnitsFilled += qty
					rev := qty * skuRetail[sku]
					transferCash(s, accCust, st, rev, -1, d, tick)
					if s.logEnabled {
						s.sales = append(s.sales, SaleRec{d, tick, ev.dst, ev.sku, ev.qty, rev})
					}
					s.branches[0]++
					if d < 6 {
						s.dailyFillStore[st][1] += qty
					}
				} else {
					part := s.acc[st].onHand[sku]
					if part > 0 {
						s.acc[st].onHand[sku] = 0
						s.totalSold[sku] += part
						s.demUnitsFilled += part
						rev := part * skuRetail[sku]
						transferCash(s, accCust, st, rev, -1, d, tick)
						if s.logEnabled {
							s.sales = append(s.sales, SaleRec{d, tick, ev.dst, ev.sku, int32(part), rev})
						}
						if d < 6 {
							s.dailyFillStore[st][1] += part
						}
					}
					s.acc[st].backorder[sku] += qty - part
					s.demBoEvents++
					s.branches[1]++
				}
				capS := int64(450)
				if safetyOvr > 0 {
					capS = safetyOvr
				}
				capS *= cs
				sLvl := (capS * 40) / 100
				invPos := s.acc[st].onHand[sku] + s.acc[st].inTransit[sku] - s.acc[st].backorder[sku]
				if invPos < sLvl && s.acc[st].status == 0 {
					s.branches[3]++
					srcWh := 10 + sku
					if sku == 3 || st == 9 {
						srcWh = 13
					}
					doShip(s, srcWh, st, sku, capS-invPos, d, tick, ev.time, bin)
				}
			case 2:
				dst, sku, qty := int(ev.dst), int(ev.sku), int64(ev.qty)
				if s.acc[dst].inTransit[sku] >= qty {
					s.acc[dst].inTransit[sku] -= qty
					s.acc[dst].onHand[sku] += qty
				}
				if s.logEnabled {
					s.ships = append(s.ships, ShipRec{d, tick, ev.src, ev.dst, ev.aux, ev.sku, ev.qty, 2, 0, 0})
				}
			case 3:
				fac, sku, qty := int(ev.src), int(ev.sku), int64(ev.qty)
				if s.acc[fac].plantState == 1 {
					qty = 0
				} else if s.acc[fac].plantState == 2 {
					qty /= 2
				}
				if qty > 0 {
					cost := qty * skuProd[sku]
					s.acc[fac].onHand[sku] += qty
					s.totalProduced[sku] += qty
					transferCash(s, fac, accSupp, cost, -1, d, tick)
					if d < 6 {
						s.dailyCostEchelon[d][0] += cost
					}
					if s.logEnabled {
						s.prods = append(s.prods, ProdRec{d, tick, ev.src, ev.sku, int32(qty), cost})
					}
					doShip(s, fac, int(ev.dst), sku, qty, d, tick, ev.time, bin)
				}
			case 4:
				m := min(s.acc[13].onHand[0], min(s.acc[13].onHand[1], s.acc[13].onHand[2]))
				mx := max(s.acc[13].onHand[0], max(s.acc[13].onHand[1], s.acc[13].onHand[2]))
				imb := mx - m
				if len(s.kitImbSeries) < 64 {
					s.kitImbSeries = append(s.kitImbSeries, imb)
				}
				q := int64(ev.qty)
				if m >= q {
					s.branches[8]++
				} else {
					q = max(m, 0)
					s.branches[9]++
				}
				if q > 0 {
					for k := 0; k < 3; k++ {
						s.acc[13].onHand[k] -= q
					}
					s.acc[13].onHand[3] += q
					s.totalKits += q
					fee := q * s.cfg.kitFee
					transferCash(s, 13, accFee, fee, -1, d, tick)
					if d < 6 {
						s.dailyCostEchelon[d][1] += fee
					}
					if s.logEnabled {
						s.kits = append(s.kits, KitRec{d, tick, 13, int32(q), int32(imb), fee})
					}
				}
			case 5:
				minK := 0
				for k := 1; k < 3; k++ {
					if s.acc[13].onHand[k] < s.acc[13].onHand[minK] {
						minK = k
					}
				}
				bestW := 10
				for w := 11; w < 13; w++ {
					if s.acc[w].onHand[minK] > s.acc[bestW].onHand[minK] {
						bestW = w
					}
				}
				s.branches[10]++
				s.nTransships++
				doShip(s, bestW, 13, minK, int64(ev.qty), d, tick, ev.time, bin)
				if s.fTrans != nil {
					fmt.Fprintf(s.fTrans, "%d,%d,%d,13,%d,%d\n", d, tick, bestW, minK, ev.qty)
				}
			case 7:
				s.acc[ev.src].plantState, s.acc[ev.src].disruptRem = 1, 1
				s.branches[11]++
				if s.fDisrupt != nil {
					fmt.Fprintf(s.fDisrupt, "%d,%d,%d,1,1\n", d, tick, ev.src)
				}
			case 6:
				cls := int(ev.aux)
				mag := int64(cls+1) * 15 * cs
				s.branches[12]++
				s.shockCnt[cls]++
				s.shockResp[cls] += mag
				if len(s.shockEvs) < 128 {
					s.shockEvs = append(s.shockEvs, [3]int{bin, cls, int(ev.dst)})
				}
				if s.fShocks != nil {
					fmt.Fprintf(s.fShocks, "%d,%d,%d,%d,%d\n", d, tick, cls, ev.dst, mag)
				}
				heapPush(&s.heap, SimEvent{ev.time + 1, s.nextSeq, 0, accCust, ev.dst, 0, int32(min(mag, 5000)), 0})
				s.nextSeq++
			case 8:
				simRiskCheck(s, d, tick)
			}
			if evProc%s.cfg.kChk == 0 {
				checkInvariants(s)
			}
		}
		simEodClearing(s, d)
		checkInvariants(s)
	}
}
