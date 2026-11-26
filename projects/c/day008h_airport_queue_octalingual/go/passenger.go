package main

type Agent struct {
	active, cls, dest, flightHr, nBags, curStn, extraSec, isBag int32
	arrSlot, qEnterSlot, depSlot                                int64
}

func poolAlloc(s *Sim) int32 {
	if s.freeTop <= 0 {
		return -1
	}
	s.freeTop--
	idx := s.freeStk[s.freeTop]
	s.pool[idx] = Agent{active: 1}
	return idx
}

func poolFree(s *Sim, idx int32) {
	if idx < 0 || int(idx) >= MaxPool || s.pool[idx].active == 0 {
		return
	}
	s.pool[idx].active = 0
	s.freeStk[s.freeTop] = idx
	s.freeTop++
}
