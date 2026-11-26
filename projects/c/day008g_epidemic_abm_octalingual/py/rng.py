import struct

MASK64 = 0xFFFFFFFFFFFFFFFF
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def splitmix64(x: int):
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    z = x
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return x, (z ^ (z >> 31)) & MASK64


def fnv1a_bytes(buf: bytes) -> int:
    h = FNV_OFFSET
    for b in buf:
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


def fnv1a_u32pair(a: int, b: int) -> int:
    return fnv1a_bytes(struct.pack("<II", a & 0xFFFFFFFF, b & 0xFFFFFFFF))


def fnv1a_i64(h: int, v: int) -> int:
    for b in struct.pack("<Q", v & MASK64):
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


class Xoshiro256:
    def __init__(self, sd: int):
        x = sd & MASK64
        self.s = [0] * 4
        for i in range(4):
            x, self.s[i] = splitmix64(x)

    def _rotl64(self, v: int, k: int) -> int:
        return ((v << k) | (v >> (64 - k))) & MASK64

    def next_u64(self) -> int:
        res = (self._rotl64((self.s[0] + self.s[3]) & MASK64, 23) + self.s[0]) & MASK64
        t = (self.s[1] << 17) & MASK64
        self.s[2] ^= self.s[0]
        self.s[3] ^= self.s[1]
        self.s[1] ^= self.s[2]
        self.s[0] ^= self.s[3]
        self.s[2] ^= t
        self.s[3] = self._rotl64(self.s[3], 45)
        return res

    def rand_int(self, lo: int, hi: int) -> int:
        if hi <= lo:
            return lo
        return lo + ((self.next_u64() >> 11) % (hi - lo + 1))
