/** Seedable Xoshiro128++ PRNG with Polar Marsaglia normal sampling. */

const MASK32 = 0xffff_ffffn;
const INV_2_32 = 1.0 / 4294967296.0;

function splitmix64(state: { value: bigint }): bigint {
  state.value = (state.value + 0x9e37_79b9_7f4a_7c15n) & 0xffff_ffff_ffff_ffffn;
  let z = state.value;
  z = ((z ^ (z >> 30n)) * 0xbf58_476d_1ce4_e5b9n) & 0xffff_ffff_ffff_ffffn;
  z = ((z ^ (z >> 27n)) * 0x94d0_49bb_1331_11ebn) & 0xffff_ffff_ffff_ffffn;
  return (z ^ (z >> 31n)) & 0xffff_ffff_ffff_ffffn;
}

export function deriveChunkSeed(seed: number, chunkIdx: number): bigint {
  const s = BigInt(seed >>> 0);
  const c = BigInt(chunkIdx + 1);
  return (s ^ (c * 0x9e37_79b9_7f4a_7c15n)) & 0xffff_ffff_ffff_ffffn;
}

function rotl32(x: number, k: number): number {
  return ((x << k) | (x >>> (32 - k))) >>> 0;
}

export class SeededRng {
  private s0: number;
  private s1: number;
  private s2: number;
  private s3: number;
  private hasSpare = false;
  private spare = 0.0;

  constructor(seed: bigint | number) {
    const state = { value: typeof seed === "bigint" ? seed : BigInt(seed >>> 0) };
    const w01 = splitmix64(state);
    const w23 = splitmix64(state);
    this.s0 = Number(w01 & MASK32) >>> 0;
    this.s1 = Number((w01 >> 32n) & MASK32) >>> 0;
    this.s2 = Number(w23 & MASK32) >>> 0;
    this.s3 = Number((w23 >> 32n) & MASK32) >>> 0;
    if ((this.s0 | this.s1 | this.s2 | this.s3) === 0) {
      this.s0 = 1;
    }
  }

  /** Returns a uniform double in [0, 1). */
  nextUniform(): number {
    const result = (rotl32((this.s0 + this.s3) >>> 0, 7) + this.s0) >>> 0;
    const t = (this.s1 << 9) >>> 0;

    this.s2 ^= this.s0;
    this.s3 ^= this.s1;
    this.s1 ^= this.s2;
    this.s0 ^= this.s3;
    this.s2 ^= t;
    this.s3 = rotl32(this.s3, 11);

    return result * INV_2_32;
  }

  /** Returns a standard normal N(0, 1) sample via Polar Marsaglia. */
  nextNormal(): number {
    if (this.hasSpare) {
      this.hasSpare = false;
      return this.spare;
    }
    let u = 0.0;
    let v = 0.0;
    let s = 0.0;
    do {
      u = 2.0 * this.nextUniform() - 1.0;
      v = 2.0 * this.nextUniform() - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s === 0.0);

    const mul = Math.sqrt((-2.0 * Math.log(s)) / s);
    this.spare = v * mul;
    this.hasSpare = true;
    return u * mul;
  }
}
