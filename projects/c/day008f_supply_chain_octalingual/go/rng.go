package main

import "math/bits"

type Xoshiro256 struct{ s [4]uint64 }

func splitmix64(x *uint64) uint64 {
	*x += 0x9E3779B97F4A7C15
	z := *x
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}
func fnv1aBytes(buf []byte) uint64 {
	h := uint64(0xCBF29CE484222325)
	for _, b := range buf {
		h = (h ^ uint64(b)) * 0x100000001B3
	}
	return h
}
func fnv1aU32Pair(a, b uint32) uint64 {
	var buf [8]byte
	for i := 0; i < 4; i++ {
		buf[i] = byte((a >> (8 * i)) & 0xFF)
		buf[4+i] = byte((b >> (8 * i)) & 0xFF)
	}
	return fnv1aBytes(buf[:])
}
func fnv1aI64(h uint64, v int64) uint64 {
	u := uint64(v)
	for i := 0; i < 8; i++ {
		h = (h ^ ((u >> (8 * i)) & 0xFF)) * 0x100000001B3
	}
	return h
}
func x256Seed(seed uint64) Xoshiro256 {
	var st Xoshiro256
	x := seed
	for i := 0; i < 4; i++ {
		st.s[i] = splitmix64(&x)
	}
	return st
}
func (st *Xoshiro256) next() uint64 {
	res := bits.RotateLeft64(st.s[0]+st.s[3], 23) + st.s[0]
	t := st.s[1] << 17
	st.s[2] ^= st.s[0]
	st.s[3] ^= st.s[1]
	st.s[1] ^= st.s[2]
	st.s[0] ^= st.s[3]
	st.s[2] ^= t
	st.s[3] = bits.RotateLeft64(st.s[3], 45)
	return res
}
func (st *Xoshiro256) randInt(lo, hi int64) int64 {
	if hi <= lo {
		return lo
	}
	span := uint64(hi - lo + 1)
	return lo + int64((st.next()>>11)%span)
}

type SimEvent struct {
	time                           int64
	seq                            uint64
	etype, src, dst, sku, qty, aux int32
}

func evLt(a, b SimEvent) bool {
	if a.time != b.time {
		return a.time < b.time
	}
	return a.seq < b.seq
}
func heapPush(h *[]SimEvent, ev SimEvent) {
	*h = append(*h, ev)
	i := len(*h) - 1
	for i > 0 {
		p := (i - 1) >> 1
		if !evLt(ev, (*h)[p]) {
			break
		}
		(*h)[i] = (*h)[p]
		i = p
	}
	(*h)[i] = ev
}
func heapPop(h *[]SimEvent) SimEvent {
	top := (*h)[0]
	last := (*h)[len(*h)-1]
	*h = (*h)[:len(*h)-1]
	if len(*h) > 0 {
		i, n := 0, len(*h)
		for {
			l, r, best := 2*i+1, 2*i+2, i
			if l < n && evLt((*h)[l], last) {
				best = l
			}
			cmp := last
			if best != i {
				cmp = (*h)[best]
			}
			if r < n && evLt((*h)[r], cmp) {
				best = r
			}
			if best == i {
				break
			}
			(*h)[i] = (*h)[best]
			i = best
		}
		(*h)[i] = last
	}
	return top
}
