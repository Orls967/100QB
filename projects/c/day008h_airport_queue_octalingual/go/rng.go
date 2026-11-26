package main

import (
	"encoding/binary"
	"math/bits"
)

const BP int64 = 10000

func splitmix64(x *uint64) uint64 {
	*x += 0x9e3779b97f4a7c15
	z := *x
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

func fnv1aBytes(buf []byte) uint64 {
	h := uint64(0xcbf29ce484222325)
	for _, b := range buf {
		h = (h ^ uint64(b)) * 0x100000001b3
	}
	return h
}

func fnv1aU32Pair(a, b uint32) uint64 {
	var buf [8]byte
	binary.LittleEndian.PutUint32(buf[0:4], a)
	binary.LittleEndian.PutUint32(buf[4:8], b)
	return fnv1aBytes(buf[:])
}

func fnv1aI64(h uint64, v int64) uint64 {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], uint64(v))
	for _, b := range buf {
		h = (h ^ uint64(b)) * 0x100000001b3
	}
	return h
}

type Xoshiro256 struct {
	s [4]uint64
}

func newXoshiro(sd uint64) Xoshiro256 {
	x := sd
	var st Xoshiro256
	for i := 0; i < 4; i++ {
		st.s[i] = splitmix64(&x)
	}
	return st
}

func (st *Xoshiro256) nextU64() uint64 {
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
	return lo + int64((st.nextU64()>>11)%uint64(hi-lo+1))
}

func (st *Xoshiro256) bernoulli(pBp int64) int32 {
	if pBp <= 0 {
		return 0
	}
	if pBp >= BP {
		return 1
	}
	if (st.nextU64()>>11)%uint64(BP) < uint64(pBp) {
		return 1
	}
	return 0
}
