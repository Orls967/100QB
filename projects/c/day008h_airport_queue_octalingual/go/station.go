package main

type Server struct {
	status, enabled, agentIdx, extraSec int32
}

type Station struct {
	q                                             [MaxQ]int32
	head, tail, len, cap, hasPrio                 int32
	cMax, cActive, busyCnt                        int32
	srv                                           [MaxSrv]Server
	arrCnt, depCnt, blockCnt, abanCnt, sumL, sumW int64
	waitCnt, busySlots                            int64
}

func initStation(st *Station, cMax, cap, hasPrio int32) {
	*st = Station{cMax: cMax, cActive: cMax, cap: cap, hasPrio: hasPrio}
	for i := 0; i < int(cMax); i++ {
		st.srv[i].enabled = 1
		st.srv[i].agentIdx = -1
	}
}

func (st *Station) pushTail(idx int32) {
	if int(st.len) >= MaxQ {
		return
	}
	st.q[st.tail] = idx
	st.tail = (st.tail + 1) % int32(MaxQ)
	st.len++
}

func (st *Station) pushHead(idx int32) {
	if int(st.len) >= MaxQ {
		return
	}
	st.head = (st.head - 1 + int32(MaxQ)) % int32(MaxQ)
	st.q[st.head] = idx
	st.len++
}

func (st *Station) removeAt(pos int32) int32 {
	idx := st.q[(st.head+pos)%int32(MaxQ)]
	for i := pos; i < st.len-1; i++ {
		st.q[(st.head+i)%int32(MaxQ)] = st.q[(st.head+i+1)%int32(MaxQ)]
	}
	st.tail = (st.tail - 1 + int32(MaxQ)) % int32(MaxQ)
	st.len--
	return idx
}

func (st *Station) popNext(pool []Agent, priorityOn int32, prioBr *int64) int32 {
	if st.len == 0 {
		return -1
	}
	if priorityOn != 0 && st.hasPrio != 0 {
		for i := int32(0); i < st.len; i++ {
			aidx := st.q[(st.head+i)%int32(MaxQ)]
			if pool[aidx].cls == ClsBiz {
				if i > 0 {
					*prioBr++
					return st.removeAt(i)
				}
				break
			}
		}
	}
	res := st.q[st.head]
	st.head = (st.head + 1) % int32(MaxQ)
	st.len--
	return res
}
