package main

import "sort"

func runAllocAuction(s *Sim, wh, sku int, stores []int, reqs, seqs []int64, prio []int32) []int64 {
	n := len(stores)
	fills := make([]int64, n)
	totReq := int64(0)
	for _, r := range reqs {
		totReq += r
	}
	if totReq <= 0 {
		return fills
	}
	avail := s.acc[wh].onHand[sku]
	s.allocCases++
	s.allocTotalReq += totReq
	if avail <= 0 {
		return fills
	}
	if avail >= totReq {
		copy(fills, reqs)
		s.allocTotalFill += totReq
		s.allocNonzero++
		s.branches[6]++
		return fills
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(i, j int) bool {
		a, b := idx[i], idx[j]
		if prio[a] != prio[b] {
			return prio[a] > prio[b]
		}
		if seqs[a] != seqs[b] {
			return seqs[a] < seqs[b]
		}
		return stores[a] < stores[b]
	})
	rem := avail
	for rem > 0 {
		maxShort := int64(-1)
		for _, i := range idx {
			if sh := reqs[i] - fills[i]; sh > maxShort {
				maxShort = sh
			}
		}
		if maxShort <= 0 {
			break
		}
		prog := false
		for _, i := range idx {
			if rem <= 0 {
				break
			}
			if reqs[i]-fills[i] == maxShort {
				fills[i]++
				rem--
				prog = true
			}
		}
		if !prog {
			break
		}
	}
	filled := avail - rem
	s.allocTotalFill += filled
	if filled > 0 {
		s.allocNonzero++
	}
	s.branches[7]++
	return fills
}
