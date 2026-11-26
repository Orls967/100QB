export const M64 = 0xffffffffffffffffn;

export function splitmix64(x: bigint): [bigint, bigint] {
  x = (x + 0x9e3779b97f4a7c15n) & M64;
  let z = x;
  z = ((z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n) & M64;
  z = ((z ^ (z >> 27n)) * 0x94d049bb133111ebn) & M64;
  return [(z ^ (z >> 31n)) & M64, x];
}

export function fnv1aBytes(buf: Uint8Array): bigint {
  let h = 0xcbf29ce484222325n;
  for (let i = 0; i < buf.length; i++) h = ((h ^ BigInt(buf[i])) * 0x100000001b3n) & M64;
  return h;
}

export function fnv1aU32Pair(a: number, b: number): bigint {
  const buf = new Uint8Array(8);
  const dv = new DataView(buf.buffer);
  dv.setUint32(0, a >>> 0, true);
  dv.setUint32(4, b >>> 0, true);
  return fnv1aBytes(buf);
}

export function fnv1aI64(h: bigint, v: number | bigint): bigint {
  let u = BigInt.asUintN(64, typeof v === "bigint" ? v : BigInt(v));
  for (let i = 0; i < 8; i++) {
    h = ((h ^ (u & 0xffn)) * 0x100000001b3n) & M64;
    u >>= 8n;
  }
  return h;
}

const rotl64 = (x: bigint, k: bigint): bigint => ((x << k) | (x >> (64n - k))) & M64;

export class Xoshiro256 {
  s: [bigint, bigint, bigint, bigint] = [0n, 0n, 0n, 0n];
  constructor(seed: bigint) {
    let x = seed & M64;
    for (let i = 0; i < 4; i++) {
      const [z, nx] = splitmix64(x);
      this.s[i] = z; x = nx;
    }
  }
  nextU64(): bigint {
    const s = this.s;
    const res = (rotl64((s[0] + s[3]) & M64, 23n) + s[0]) & M64;
    const t = (s[1] << 17n) & M64;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl64(s[3], 45n);
    return res;
  }
  randInt(lo: number, hi: number): number {
    if (hi <= lo) return lo;
    const span = BigInt(hi - lo + 1);
    return lo + Number((this.nextU64() >> 11n) % span);
  }
}

export interface SimEvent { time: number; seq: number; etype: number; src: number; dst: number; sku: number; qty: number; aux: number; }
const evLt = (a: SimEvent, b: SimEvent): boolean => a.time !== b.time ? a.time < b.time : a.seq < b.seq;

export class MinHeap {
  h: SimEvent[] = [];
  empty(): boolean { return this.h.length === 0; }
  peekTime(): number { return this.h[0].time; }
  push(ev: SimEvent): void {
    this.h.push(ev);
    let i = this.h.length - 1;
    while (i > 0) {
      const p = (i - 1) >> 1;
      if (!evLt(ev, this.h[p])) break;
      this.h[i] = this.h[p]; i = p;
    }
    this.h[i] = ev;
  }
  pop(): SimEvent {
    const top = this.h[0], last = this.h.pop()!;
    if (this.h.length > 0) {
      let i = 0; const n = this.h.length;
      for (;;) {
        const l = 2 * i + 1, r = l + 1; let best = i;
        if (l < n && evLt(this.h[l], last)) best = l;
        if (r < n && evLt(this.h[r], best === i ? last : this.h[best])) best = r;
        if (best === i) break;
        this.h[i] = this.h[best]; i = best;
      }
      this.h[i] = last;
    }
    return top;
  }
}
