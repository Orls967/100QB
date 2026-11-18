import { describe, expect, it } from "vitest";
import { CHUNK_SIZE, runChunk } from "../src/core.js";
import {
  type OptionParams,
  bsCall,
  bsPut,
  deltaAsianBump,
  deltaAsianPathwise,
  geometricAsianCall,
  priceAsianControlVariate,
  priceAsianNaive,
} from "../src/options.js";
import { SeededRng } from "../src/rng.js";

function defaultParams(): OptionParams {
  return {
    s0: 100.0,
    k: 100.0,
    r: 0.05,
    sigma: 0.2,
    t: 1.0,
    n_steps: 252,
  };
}

describe("day001f_asian_options", () => {
  it("satisfies Black-Scholes put-call parity", () => {
    const p = { ...defaultParams(), k: 105.0 };
    const c = bsCall(p);
    const put = bsPut(p);
    const rhs = p.s0 - p.k * Math.exp(-p.r * p.t);
    expect(Math.abs(c - put - rhs)).toBeLessThan(1e-10);
  });

  it("matches geometric Asian closed form with large MC", () => {
    const p = defaultParams();
    const exact = geometricAsianCall(p);
    const nPaths = 40_000;
    const nChunks = Math.ceil(nPaths / CHUNK_SIZE);
    let sy = 0.0;
    let syy = 0.0;
    for (let c = 0; c < nChunks; c += 1) {
      const len = Math.min(nPaths - c * CHUNK_SIZE, CHUNK_SIZE);
      const acc = runChunk("cv", p, c, len, 42);
      sy += acc[1] ?? 0;
      syy += acc[3] ?? 0;
    }
    const mean = sy / nPaths;
    const se = Math.sqrt((syy - nPaths * mean * mean) / (nPaths - 1) / nPaths);
    expect(Math.abs(mean - exact)).toBeLessThan(3.0 * se);
  });

  it("gives much smaller SE with control variate than naive MC", async () => {
    const p = defaultParams();
    const naive = await priceAsianNaive(p, 20_000, 42);
    const cv = await priceAsianControlVariate(p, 20_000, 42);
    expect(cv.se).toBeLessThan(naive.se / 10.0);
  });

  it("is deterministic across same seed and 1 vs 4 workers", async () => {
    const p = { ...defaultParams(), n_steps: 50 };
    const a = await priceAsianControlVariate(p, 10_000, 7, 1);
    const b = await priceAsianControlVariate(p, 10_000, 7, 4);
    expect(a.price).toBe(b.price);
    expect(a.se).toBe(b.se);
  });

  it("matches pathwise delta and bump-and-revalue within 3 SE", async () => {
    const p = defaultParams();
    const pw = await deltaAsianPathwise(p, 50_000, 42);
    const bmp = await deltaAsianBump(p, 50_000, 0.01, 42);
    const seComb = Math.hypot(pw.se, bmp.se);
    expect(Math.abs(pw.delta - bmp.delta)).toBeLessThan(3.0 * seComb);
  });

  it("produces standard normal samples with mean ~ 0 and var ~ 1", () => {
    const rng = new SeededRng(42);
    const n = 40_000;
    let sum = 0.0;
    let sumSq = 0.0;
    for (let i = 0; i < n; i += 1) {
      const z = rng.nextNormal();
      sum += z;
      sumSq += z * z;
    }
    const mean = sum / n;
    const variance = sumSq / n - mean * mean;
    expect(Math.abs(mean)).toBeLessThan(0.02);
    expect(Math.abs(variance - 1.0)).toBeLessThan(0.03);
  });
});
