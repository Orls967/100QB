package main

func simRun(s *Sim) {
	mult := int64(1)
	if s.rp.timescaleDiv2 != 0 {
		mult = 2
	}
	totalSlots := s.rp.days * SlotsPerDay * mult
	for t := int64(0); t < totalSlots; t++ {
		simStepSlot(s, t)
	}
	remPax := int64(0)
	for i := 0; i < MaxPool; i++ {
		if s.pool[i].active != 0 && s.pool[i].isBag == 0 {
			remPax++
			s.sysSumSojourn += totalSlots - s.pool[i].arrSlot
		}
	}
	stnRem := int64(0)
	for st := 0; st < NStn; st++ {
		if st != StnBelt {
			stnRem += int64(s.stn[st].len + s.stn[st].busyCnt)
		}
	}
	if s.totalArr != s.totalDone+s.totalLost+remPax || stnRem != remPax {
		s.invViolations++
	}
	if s.sysSumL != s.sysSumSojourn {
		s.invViolations++
	}
}
