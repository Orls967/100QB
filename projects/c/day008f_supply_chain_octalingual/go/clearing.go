package main

import "fmt"

func simRiskCheck(s *Sim, day, tick int64) {
	for i := 0; i < numStores+numWarehouses; i++ {
		boVal := int64(0)
		for k := 0; k < numSKUs; k++ {
			boVal += s.acc[i].backorder[k] * skuWholesale[k]
		}
		reqMargin := (boVal * s.cfg.creditMarginBps) / 10000
		if s.acc[i].cash < reqMargin && s.acc[i].status == 0 {
			s.acc[i].status = 1
			s.branches[13]++
			if s.fCredit != nil {
				fmt.Fprintf(s.fCredit, "%d,%d,%d,%d,1\n", day, i, s.acc[i].cash, reqMargin)
			}
		} else if s.acc[i].status == 1 && s.acc[i].cash >= reqMargin {
			s.acc[i].status = 0
		}
		if s.acc[i].cash < 0 && s.acc[i].status != 2 {
			s.acc[i].status = 2
			s.branches[14]++
			transferCash(s, i, accFee, 50000, 2, day, tick)
			if s.fRisk != nil {
				fmt.Fprintf(s.fRisk, "%d,%d,%d,LIQUIDATE,50000\n", day, tick, i)
			}
		}
	}
}
func simEodClearing(s *Sim, day int64) {
	for i := 0; i < numStores+numWarehouses; i++ {
		ohTot, itTot, boTot := int64(0), int64(0), int64(0)
		for k := 0; k < numSKUs; k++ {
			ohTot += s.acc[i].onHand[k]
			itTot += s.acc[i].inTransit[k]
			boTot += s.acc[i].backorder[k]
			if day < 6 {
				s.dayInvSnap[day][k] += s.acc[i].onHand[k]
			}
		}
		hold, pen := ohTot*s.cfg.holdingCostPerDay, boTot*s.cfg.backorderPenalty
		if hold > 0 {
			transferCash(s, i, accFee, hold, 0, day, s.cfg.ticksPerDay)
			if day < 6 {
				s.dailyCostEchelon[day][3] += hold
			}
		}
		if pen > 0 {
			transferCash(s, i, accFee, pen, 1, day, s.cfg.ticksPerDay)
			if day < 6 {
				s.dailyCostEchelon[day][3] += pen
			}
		}
		if s.fAccSer != nil {
			fmt.Fprintf(s.fAccSer, "%d,%d,%d,%d,%d,%d\n", day, i, s.acc[i].cash, ohTot, itTot, s.acc[i].status)
		}
	}
	simRiskCheck(s, day, s.cfg.ticksPerDay)
}
