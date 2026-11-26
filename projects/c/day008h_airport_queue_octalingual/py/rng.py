import struct

BP = 10000
MASK64 = 0xFFFFFFFFFFFFFFFF
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def splitmix64(x: int) -> tuple[int, int]:
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    z = x
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return x, (z ^ (z >> 31)) & MASK64


def fnv1a_bytes(data: bytes) -> int:
    h = FNV_OFFSET
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


def fnv1a_u32pair(a: int, b: int) -> int:
    return fnv1a_bytes(struct.pack("<II", a & 0xFFFFFFFF, b & 0xFFFFFFFF))


def fnv1a_i64(h: int, v: int) -> int:
    return fnv1a_bytes_with(h, struct.pack("<Q", v & MASK64))


def fnv1a_bytes_with(h: int, data: bytes) -> int:
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & MASK64
    return h


class Xoshiro256:
    def __init__(self, sd: int) -> None:
        x = sd & MASK64
        self.s = [0, 0, 0, 0]
        for i in range(4):
            x, v = splitmix64(x)
            self.s[i] = v

    @staticmethod
    def rotl64(v: int, k: int) -> int:
        return ((v << k) | (v >> (64 - k))) & MASK64

    def next_u64(self) -> int:
        s = self.s
        res = (self.rotl64((s[0] + s[3]) & MASK64, 23) + s[0]) & MASK64
        t = (s[1] << 17) & MASK64
        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]
        s[2] ^= t
        s[3] = self.rotl64(s[3], 45)
        return res

    def rand_int(self, lo: int, hi: int) -> int:
        if hi <= lo:
            return lo
        return lo + ((self.next_u64() >> 11) % (hi - lo + 1))

    def bernoulli(self, p_bp: int) -> int:
        if p_bp <= 0:
            return 0
        if p_bp >= BP:
            return 1
        return 1 if ((self.next_u64() >> 11) % BP) < p_bp else 0
