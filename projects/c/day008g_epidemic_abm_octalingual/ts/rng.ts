export const MASK64 = 0xffffffffffffffffn;
export const FNV_OFFSET = 0xcbf29ce484222325n;
export const FNV_PRIME = 0x100000001b3n;

export function splitmix64(x: bigint): [bigint, bigint] {
  const nx = (x + 0x9e3779b97f4a7c15n) & MASK64;
  let z = nx;
  z = ((z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n) & MASK64;
  z = ((z ^ (z >> 27n)) * 0x94d049bb133111ebn) & MASK64;
  return [nx, (z ^ (z >> 31n)) & MASK64];
}

export function fnv1aBytes(buf: Uint8Array): bigint {
  let h = FNV_OFFSET;
  for (let i = 0; i < buf.length; i++) {
    h = ((h ^ BigInt(buf[i])) * FNV_PRIME) & MASK64;
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

export function fnv1aU64(h: bigint, v: bigint): bigint {
  let uv = v & MASK64;
  for (let i = 0; i < 8; i++) {
    h = ((h ^ (uv & 0xffn)) * FNV_PRIME) & MASK64;
    uv >>= 8n;
  }
  return h;
}

export function fnv1aI64(h: bigint, v: number): bigint {
  return fnv1aU64(h, BigInt(v));
}

export function hex64(v: bigint): string {
  return (v & MASK64).toString(16).padStart(16, "0");
}

export class Xoshiro256 {
  s: bigint[];

  constructor(sd: bigint) {
    let x = sd & MASK64;
    this.s = [0n, 0n, 0n, 0n];
    for (let i = 0; i < 4; i++) {
      const [nx, val] = splitmix64(x);
      x = nx;
      this.s[i] = val;
    }
  }

  private rotl64(v: bigint, k: bigint): bigint {
    return ((v << k) | (v >> (64n - k))) & MASK64;
  }

  nextU64(): bigint {
    const res = (this.rotl64((this.s[0] + this.s[3]) & MASK64, 23n) + this.s[0]) & MASK64;
    const t = (this.s[1] << 17n) & MASK64;
    this.s[2] ^= this.s[0];
    this.s[3] ^= this.s[1];
    this.s[1] ^= this.s[2];
    this.s[0] ^= this.s[3];
    this.s[2] ^= t;
    this.s[3] = this.rotl64(this.s[3], 45n);
    return res;
  }

  randInt(lo: number, hi: number): number {
    if (hi <= lo) return lo;
    const span = BigInt(hi - lo + 1);
    return lo + Number((this.nextU64() >> 11n) % span);
  }
}
