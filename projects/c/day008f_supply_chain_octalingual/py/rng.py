import heapq, struct

M64 = 0xFFFFFFFFFFFFFFFF

def splitmix64(x):
    x = (x + 0x9E3779B97F4A7C15) & M64
    z = x
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return (z ^ (z >> 31)) & M64, x

def fnv1a_bytes(bs):
    h = 0xCBF29CE484222325
    for b in bs:
        h = ((h ^ b) * 0x100000001B3) & M64
    return h

def fnv1a_u32pair(a, b):
    return fnv1a_bytes(struct.pack("<II", a, b))

def fnv1a_i64(h, v):
    u = v & M64
    for _ in range(8):
        h = ((h ^ (u & 0xFF)) * 0x100000001B3) & M64
        u >>= 8
    return h

def rotl64(x, k):
    return ((x << k) | (x >> (64 - k))) & M64

class Xoshiro256:
    def __init__(self, seed):
        self.s = []
        x = seed & M64
        for _ in range(4):
            z, x = splitmix64(x)
            self.s.append(z)

    def next_u64(self):
        s = self.s
        res = (rotl64((s[0] + s[3]) & M64, 23) + s[0]) & M64
        t = (s[1] << 17) & M64
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]
        s[2] ^= t; s[3] = rotl64(s[3], 45)
        return res

    def rand_int(self, lo, hi):
        if hi <= lo: return lo
        return lo + ((self.next_u64() >> 11) % (hi - lo + 1))

class MinHeap:
    def __init__(self): self.h = []
    def empty(self): return not self.h
    def peek_time(self): return self.h[0][0]
    def push(self, ev): heapq.heappush(self.h, ev)
    def pop(self): return heapq.heappop(self.h)
