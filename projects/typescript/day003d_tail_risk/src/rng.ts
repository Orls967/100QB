/**
 * Seedable 64-bit xoshiro256++ PRNG implemented with 32-bit unsigned limbs
 * for bit-exact equivalence to 64-bit BigInt xoshiro256++ without BigInt heap churn,
 * paired with Box-Muller normals and exact Marsaglia-Tsang Gamma / Student-t sampling.
 */

const INV_2_32 = 1.0 / 4294967296.0;
const INV_2_53 = 1.0 / 9007199254740992.0; // 2^-53
const POW_2_21 = 2097152.0; // 2^21
const TWO_PI = 2.0 * Math.PI;
const U64_MASK = 0xffff_ffff_ffff_ffffn;

/** SplitMix64 mixer over BigInt for initializing state vectors cleanly from a 64-bit seed. */
export function splitmix64(x: bigint): bigint {
  let z = (x + 0x9e37_79b9_7f4a_7c15n) & U64_MASK;
  z = ((z ^ (z >> 30n)) * 0xbf58_476d_1ce4_e5b9n) & U64_MASK;
  z = ((z ^ (z >> 27n)) * 0x94d0_49bb_1331_11ebn) & U64_MASK;
  return (z ^ (z >> 31n)) & U64_MASK;
}

/** Deterministic per-chunk 64-bit seed derived from main seed, task salt, and chunk index. */
export function deriveChunkSeed(
  mainSeed: number | bigint,
  taskSalt: number,
  chunkIdx: number,
): bigint {
  const s = typeof mainSeed === "bigint" ? mainSeed & U64_MASK : BigInt(mainSeed >>> 0);
  const tag = (BigInt(taskSalt >>> 0) << 32n) | BigInt(chunkIdx >>> 0);
  return splitmix64(s ^ splitmix64(tag));
}

function rotl64LoHi(lo: number, hi: number, k: number): [number, number] {
  if (k === 0) {
    return [lo >>> 0, hi >>> 0];
  }
  if (k === 32) {
    return [hi >>> 0, lo >>> 0];
  }
  if (k < 32) {
    const newLo = ((lo << k) | (hi >>> (32 - k))) >>> 0;
    const newHi = ((hi << k) | (lo >>> (32 - k))) >>> 0;
    return [newLo, newHi];
  }
  const r = k - 32;
  const newLo = ((hi << r) | (lo >>> (32 - r))) >>> 0;
  const newHi = ((lo << r) | (hi >>> (32 - r))) >>> 0;
  return [newLo, newHi];
}

function add64LoHi(aLo: number, aHi: number, bLo: number, bHi: number): [number, number] {
  const loSum = aLo + bLo;
  const lo = loSum >>> 0;
  const carry = loSum >= 4294967296.0 ? 1 : 0;
  const hi = (aHi + bHi + carry) >>> 0;
  return [lo, hi];
}

/** Fast bit-exact 64-bit xoshiro256++ PRNG with Box-Muller, Gamma, and Student-t samplers. */
export class Xoshiro256PlusPlus {
  private s0Lo: number;
  private s0Hi: number;
  private s1Lo: number;
  private s1Hi: number;
  private s2Lo: number;
  private s2Hi: number;
  private s3Lo: number;
  private s3Hi: number;
  private hasSpare = false;
  private spare = 0.0;

  constructor(seed: bigint | number) {
    let z = typeof seed === "bigint" ? seed & U64_MASK : BigInt(seed >>> 0);
    const words: bigint[] = [];
    for (let i = 0; i < 4; i += 1) {
      z = splitmix64(z);
      words.push(z);
    }
    const w0 = words[0] ?? 1n;
    const w1 = words[1] ?? 1n;
    const w2 = words[2] ?? 1n;
    const w3 = words[3] ?? 1n;
    this.s0Lo = Number(w0 & 0xffff_ffffn) >>> 0;
    this.s0Hi = Number((w0 >> 32n) & 0xffff_ffffn) >>> 0;
    this.s1Lo = Number(w1 & 0xffff_ffffn) >>> 0;
    this.s1Hi = Number((w1 >> 32n) & 0xffff_ffffn) >>> 0;
    this.s2Lo = Number(w2 & 0xffff_ffffn) >>> 0;
    this.s2Hi = Number((w2 >> 32n) & 0xffff_ffffn) >>> 0;
    this.s3Lo = Number(w3 & 0xffff_ffffn) >>> 0;
    this.s3Hi = Number((w3 >> 32n) & 0xffff_ffffn) >>> 0;
    if (
      (this.s0Lo |
        this.s0Hi |
        this.s1Lo |
        this.s1Hi |
        this.s2Lo |
        this.s2Hi |
        this.s3Lo |
        this.s3Hi) ===
      0
    ) {
      this.s0Lo = 1;
    }
  }

  /** Advance state by one step and return raw 64-bit output as [lo32, hi32]. */
  nextUint64Parts(): [number, number] {
    const [sumLo, sumHi] = add64LoHi(this.s0Lo, this.s0Hi, this.s3Lo, this.s3Hi);
    const [rotLo, rotHi] = rotl64LoHi(sumLo, sumHi, 23);
    const [resLo, resHi] = add64LoHi(rotLo, rotHi, this.s0Lo, this.s0Hi);

    const tLo = (this.s1Lo << 17) >>> 0;
    const tHi = ((this.s1Hi << 17) | (this.s1Lo >>> 15)) >>> 0;

    this.s2Lo = (this.s2Lo ^ this.s0Lo) >>> 0;
    this.s2Hi = (this.s2Hi ^ this.s0Hi) >>> 0;
    this.s3Lo = (this.s3Lo ^ this.s1Lo) >>> 0;
    this.s3Hi = (this.s3Hi ^ this.s1Hi) >>> 0;
    this.s1Lo = (this.s1Lo ^ this.s2Lo) >>> 0;
    this.s1Hi = (this.s1Hi ^ this.s2Hi) >>> 0;
    this.s0Lo = (this.s0Lo ^ this.s3Lo) >>> 0;
    this.s0Hi = (this.s0Hi ^ this.s3Hi) >>> 0;

    this.s2Lo = (this.s2Lo ^ tLo) >>> 0;
    this.s2Hi = (this.s2Hi ^ tHi) >>> 0;

    const [r3Lo, r3Hi] = rotl64LoHi(this.s3Lo, this.s3Hi, 45);
    this.s3Lo = r3Lo;
    this.s3Hi = r3Hi;

    return [resLo, resHi];
  }

  /** Uniform double in [0, 1) using the upper 53 bits of the 64-bit output. */
  nextUniform(): number {
    const [lo, hi] = this.nextUint64Parts();
    return (hi * POW_2_21 + (lo >>> 11)) * INV_2_53;
  }

  /** Uniform double strictly in (0, 1) safe for log and quantile evaluation. */
  nextUniformOpen(): number {
    const [lo, hi] = this.nextUint64Parts();
    return (hi * POW_2_21 + (lo >>> 11) + 0.5) * INV_2_53;
  }

  /** Unsigned 32-bit integer from the high half of next output. */
  nextUint32(): number {
    const [, hi] = this.nextUint64Parts();
    return hi >>> 0;
  }

  /** Bernoulli draw with probability p. */
  nextBernoulli(p: number): number {
    return this.nextUniform() < p ? 1 : 0;
  }

  /** Standard normal N(0, 1) variate via classic Box-Muller transform. */
  nextNormal(): number {
    if (this.hasSpare) {
      this.hasSpare = false;
      return this.spare;
    }
    let u1 = this.nextUniform();
    if (u1 <= 0.0) {
      u1 = INV_2_32 * INV_2_32;
    }
    const u2 = this.nextUniform();
    const r = Math.sqrt(-2.0 * Math.log(u1));
    const theta = TWO_PI * u2;
    this.spare = r * Math.sin(theta);
    this.hasSpare = true;
    return r * Math.cos(theta);
  }

  /** Gamma(shape, scale=1) variate via Marsaglia-Tsang (2000). */
  nextGamma(shape: number): number {
    if (!(shape > 0.0)) {
      throw new Error(`Gamma shape must be positive, got ${shape}`);
    }
    if (shape < 1.0) {
      const u = this.nextUniformOpen();
      return this.nextGamma(shape + 1.0) * Math.pow(u, 1.0 / shape);
    }
    const d = shape - 1.0 / 3.0;
    const c = 1.0 / Math.sqrt(9.0 * d);
    while (true) {
      let x = 0.0;
      let v = 0.0;
      do {
        x = this.nextNormal();
        v = 1.0 + c * x;
      } while (v <= 0.0);
      v = v * v * v;
      const u = this.nextUniformOpen();
      const x2 = x * x;
      if (u < 1.0 - 0.0331 * x2 * x2) {
        return d * v;
      }
      if (Math.log(u) < 0.5 * x2 + d * (1.0 - v + Math.log(v))) {
        return d * v;
      }
    }
  }

  /** Chi-square(df) variate = 2 * Gamma(df / 2, 1). */
  nextChiSquare(df: number): number {
    return 2.0 * this.nextGamma(0.5 * df);
  }

  /** Raw Student-t(nu) variate with variance nu/(nu-2). */
  nextStudentT(nu: number): number {
    const z = this.nextNormal();
    const chi2 = this.nextChiSquare(nu);
    return z / Math.sqrt(chi2 / nu);
  }

  /** Standardized Student-t(nu) variate scaled to variance 1 (requires nu > 2). */
  nextStdT(nu: number): number {
    if (!(nu > 2.0)) {
      throw new Error(`Standardized t requires nu > 2, got ${nu}`);
    }
    return this.nextStudentT(nu) * Math.sqrt((nu - 2.0) / nu);
  }
}

/**
 * Pure BigInt reference implementation of 64-bit xoshiro256++ exported for unit testing
 * bit-for-bit equality against Xoshiro256PlusPlus.
 */
export function xoshiro256ppBigIntSequence(seed: bigint | number, count: number): bigint[] {
  let z = typeof seed === "bigint" ? seed & U64_MASK : BigInt(seed >>> 0);
  const s = [0n, 0n, 0n, 0n];
  for (let i = 0; i < 4; i += 1) {
    z = splitmix64(z);
    s[i] = z;
  }
  const rotl = (x: bigint, k: bigint): bigint => ((x << k) | (x >> (64n - k))) & U64_MASK;
  const out: bigint[] = new Array<bigint>(count);
  for (let i = 0; i < count; i += 1) {
    const s0 = s[0] ?? 0n;
    const s1 = s[1] ?? 0n;
    const s2 = s[2] ?? 0n;
    const s3 = s[3] ?? 0n;
    const res = (rotl((s0 + s3) & U64_MASK, 23n) + s0) & U64_MASK;
    out[i] = res;

    const t = (s1 << 17n) & U64_MASK;
    const n2 = s2 ^ s0;
    const n3 = s3 ^ s1;
    const n1 = s1 ^ n2;
    const n0 = s0 ^ n3;
    s[0] = n0;
    s[1] = n1;
    s[2] = n2 ^ t;
    s[3] = rotl(n3, 45n);
  }
  return out;
}
