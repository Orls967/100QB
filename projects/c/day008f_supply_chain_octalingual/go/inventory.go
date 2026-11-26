package main

const (
	numSKUs       = 4
	numFactories  = 3
	numWarehouses = 4
	numStores     = 10
	numCarriers   = 3
	numAccounts   = 23
	accCust       = 20
	accSupp       = 21
	accFee        = 22
)

var (
	skuRetail     = [4]int64{10000, 10200, 9800, 31500}
	skuWholesale  = [4]int64{8000, 8200, 7800, 26000}
	skuProd       = [4]int64{6000, 6200, 5800, 0}
	storeWh       = [10]int{10, 10, 10, 11, 11, 11, 12, 12, 12, 12}
	carrierBaseLt = [3]int64{120, 200, 160}
)

func carrierRate(dayVol int64) int64 {
	if dayVol < 500 {
		return 200
	}
	if dayVol < 2000 {
		return 160
	}
	return 120
}
func transferCash(s *Sim, fromID, toID int, amt int64, kind int32, day, tick int64) {
	if amt <= 0 || fromID == toID {
		return
	}
	s.acc[fromID].cash -= amt
	s.acc[toID].cash += amt
	if kind >= 0 {
		s.fees = append(s.fees, FeeRec{day, tick, int32(fromID), int32(toID), kind, amt})
	}
}
func checkInvariants(s *Sim) int {
	totCash := int64(0)
	for i := 0; i < numAccounts; i++ {
		totCash += s.acc[i].cash
	}
	if totCash != s.initCash {
		s.invViolations++
		return 0
	}
	for k := 0; k < numSKUs; k++ {
		oh, it := int64(0), int64(0)
		for i := 0; i < numAccounts; i++ {
			if s.acc[i].onHand[k] < 0 || s.acc[i].inTransit[k] < 0 || s.acc[i].backorder[k] < 0 {
				s.invViolations++
				return 0
			}
			oh += s.acc[i].onHand[k]
			it += s.acc[i].inTransit[k]
		}
		lhs := oh + it + s.totalSold[k]
		rhs := s.initStock[k]
		if k < 3 {
			lhs += s.totalKits
			rhs += s.totalProduced[k]
		} else {
			rhs += s.totalKits
		}
		if lhs != rhs {
			s.invViolations++
			return 0
		}
	}
	return 1
}
