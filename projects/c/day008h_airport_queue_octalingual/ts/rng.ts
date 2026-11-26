export const BP = 10000;
export const MASK64 = 0xffffffffffffffffn;
export const FNV_OFFSET = 0xcbf29ce484222325n;
export const FNV_PRIME = 0x100000001b3n;

export function splitmix64(xIn: bigint): [bigint, bigint] {
  const x = (xIn + 0x9e3779b97f4a7c15n) & MASK64;
  let z = x;
  z = ((z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n) & MASK64;
  z = ((z ^ (z >> 27n)) * 0x94d049bb133111ebn) & MASK64;
  return [x, (z ^ (z >> 31n)) & MASK64];
}

export function fnv1aBytes(data: Uint8Array): bigint {
  let h = FNV_OFFSET;
  for (let i = 0; i < data.length; i++) {
    h = ((h ^ BigInt(data[i])) * FNV_PRIME) & MASK64;
  }
  return h;
}

export function fnv1aU32Pair(a: number, b: number): bigint {
  const buf = new Uint8Array(8);
  const dv = new DataView(buf.buffer);
  dv.setUint32(0, a >>> 0, true);
  dv.setUint32(4, b >>> 0, true);
  return fnv1aBytes(buf);
}

export function fnv1aU64(hIn: bigint, vIn: bigint): bigint {
  let h = hIn;
  let v = vIn & MASK64;
  for (let i = 0; i < 8; i++) {
    h = ((h ^ (v & 0xffn)) * FNV_PRIME) & MASK64;
    v >>= 8n;
  }
  return h;
}

export function fnv1aI64(h: bigint, v: number): bigint {
  return fnv1aU64(h, BigInt.asUintN(64, BigInt(v)));
}

export function hex64(v: bigint): string {
  return (v & MASK64).toString(16).padStart(16, "0");
}

export class Xoshiro256 {
  s: [bigint, bigint, bigint, bigint];

  constructor(sd: bigint) {
    let x = sd & MASK64;
    const s: [bigint, bigint, bigint, bigint] = [0n, 0n, 0n, 0n];
    for (let i = 0; i < 4; i++) {
      const [nx, v] = splitmix64(x);
      x = nx;
      s[i] = v;
    }
    this.s = s;
  }

  static rotl64(v: bigint, k: bigint): bigint {
    return ((v << k) | (v >> (64n - k))) & MASK64;
  }

  nextU64(): bigint {
    const s = this.s;
    const res = (Xoshiro256.rotl64((s[0] + s[3]) & MASK64, 23n) + s[0]) & MASK64;
    const t = (s[1] << 17n) & MASK64;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = Xoshiro256.rotl64(s[3], 45n);
    return res;
  }

  randInt(lo: number, hi: number): number {
    if (hi <= lo) return lo;
    const span = BigInt(hi - lo + 1);
    return lo + Number((this.nextU64() >> 11n) % span);
  }

  bernoulli(pBp: number): number {
    if (pBp <= 0) return 0;
    if (pBp >= BP) return 1;
    return Number((this.nextU64() >> 11n) % BigInt(BP)) < pBp ? 1 : 0;
  }
}
