import { describe, expect, it } from "vitest";
import {
  gpdCdf,
  gpdPdf,
  gpdPotEs,
  gpdPotVar,
  gpdQuantile,
  locationScaleTCdf,
  locationScaleTQuantile,
  normalEs,
  stdTCdf,
  stdTEs,
  stdTEsNumerical,
  stdTPdf,
  stdTQuantile,
  studentTCdf,
  studentTCdfNumerical,
  studentTQuantile,
} from "../src/distributions.js";
import { fitGpdMle, hillEstimator, potFromSample } from "../src/evt.js";
import { Xoshiro256PlusPlus, xoshiro256ppBigIntSequence } from "../src/rng.js";
import {
  adaptiveSimpson,
  chi2Sf,
  digamma,
  erfc,
  logGamma,
  normalCdf,
  normalQuantile,
  normalSf,
} from "../src/special.js";

describe("PRNG xoshiro256++ and special functions", () => {
  it("matches 64-bit BigInt xoshiro256++ bit-for-bit over 5000 steps", () => {
    const seed = 20261121;
    const rng = new Xoshiro256PlusPlus(seed);
    const ref = xoshiro256ppBigIntSequence(seed, 5000);
    for (let i = 0; i < 5000; i += 1) {
      const [lo, hi] = rng.nextUint64Parts();
      const got = (BigInt(hi) << 32n) | BigInt(lo);
      expect(got).toBe(ref[i]);
    }
  });

  it("matches known special values for logGamma, digamma, normal and chi2Sf", () => {
    expect(Math.abs(Math.exp(logGamma(0.5)) - Math.sqrt(Math.PI))).toBeLessThan(1e-14);
    expect(Math.abs(Math.exp(logGamma(6)) - 120)).toBeLessThan(1e-11);
    const eulerGamma = 0.5772156649015329;
    expect(Math.abs(digamma(1.0) + eulerGamma)).toBeLessThan(1e-13);
    expect(Math.abs(erfc(0.0) - 1.0)).toBeLessThan(1e-15);
    expect(Math.abs(normalCdf(normalQuantile(0.999)) - 0.999)).toBeLessThan(1e-14);
    expect(Math.abs(normalSf(3.0) - (1.0 - normalCdf(3.0)))).toBeLessThan(1e-15);
    expect(Math.abs(chi2Sf(3.841458820694124, 1.0) - 0.05)).toBeLessThan(1e-12);
    expect(Math.abs(chi2Sf(5.991464547107979, 2.0) - 0.05)).toBeLessThan(1e-12);
  });
});

describe("Student-t, standardized ES, GPD, and EVT estimators", () => {
  it("matches adaptive quadrature for Student-t CDF (< 1e-10) and quantile roundtrip", () => {
    for (const nu of [3.5, 5, 10, 30]) {
      for (const x of [-4.5, -1.7, 0.4, 2.5, 5.0]) {
        const cdfAnalytic = studentTCdf(x, nu);
        const cdfNum = studentTCdfNumerical(x, nu, 1e-13);
        expect(Math.abs(cdfAnalytic - cdfNum)).toBeLessThan(1e-10);
      }
      for (const p of [0.001, 0.05, 0.5, 0.95, 0.99, 0.999]) {
        const q = studentTQuantile(p, nu);
        expect(Math.abs(studentTCdf(q, nu) - p)).toBeLessThan(1e-12);
        const qStd = stdTQuantile(p, nu);
        expect(Math.abs(stdTCdf(qStd, nu) - p)).toBeLessThan(1e-12);
      }
    }
    // Location-scale Student-t consistency
    const loc = 0.001;
    const sc = 0.02;
    const qLs = locationScaleTQuantile(0.975, loc, sc, 5, true);
    expect(Math.abs(locationScaleTCdf(qLs, loc, sc, 5, true) - 0.975)).toBeLessThan(1e-12);
  });

  it("verifies exact standardized Student-t ES formula against numerical integration (< 1e-9)", () => {
    const nus = [3.5, 5, 10, 30];
    const ps = [0.95, 0.975, 0.99, 0.999];
    for (const nu of nus) {
      for (const p of ps) {
        const esClosed = stdTEs(p, nu);
        const esQuad = stdTEsNumerical(p, nu, 1e-12);
        expect(Math.abs(esClosed - esQuad)).toBeLessThan(1e-9);
      }
    }
    // High-nu ES approaches normal ES
    expect(Math.abs(stdTEs(0.99, 1000) - normalEs(0.99))).toBeLessThan(5e-3);
    // Direct integration of x * stdTPdf(x, 5) / (1 - p) on [q95, infinity) via x = q95 / t
    const q95 = stdTQuantile(0.95, 5);
    const tailMoment =
      adaptiveSimpson(
        (t) => {
          if (t <= 0.0) return 0.0;
          const x = q95 / t;
          return (q95 / (t * t)) * x * stdTPdf(x, 5);
        },
        0.0,
        1.0,
        1e-12,
        16,
      ) /
      (1.0 - 0.95);
    expect(Math.abs(tailMoment - stdTEs(0.95, 5))).toBeLessThan(1e-10);
  });

  it("verifies GPD PDF/CDF/quantile/ES formulas numerically across shapes", () => {
    const u = 1.5;
    const pExceed = 0.1;
    for (const xi of [-0.1, 0.0, 0.2, 0.4]) {
      const beta = 0.6;
      const y = 0.9;
      const numCdf = adaptiveSimpson((t) => gpdPdf(t, xi, beta), 0.0, y, 1e-12);
      expect(Math.abs(numCdf - gpdCdf(y, xi, beta))).toBeLessThan(1e-10);
      const qVal = gpdQuantile(0.75, xi, beta);
      expect(Math.abs(gpdCdf(qVal, xi, beta) - 0.75)).toBeLessThan(1e-12);

      // Verify ES = E[X | X > VaR_p] = (VaR_p + beta - xi * u) / (1 - xi) by tail quantile quadrature
      const var99 = gpdPotVar(u, pExceed, xi, beta, 0.99);
      const esFormula = gpdPotEs(u, var99, xi, beta);
      const condTailSurv = (1.0 - 0.99) / pExceed; // P(Y > var99 - u) = 0.1
      const conditionalTailMean =
        u +
        adaptiveSimpson(
          (v) => {
            if (v <= 0.0) return 0.0;
            const v5 = v ** 5;
            const tailProb = condTailSurv * v5 * v;
            if (tailProb <= 1e-300) return 0.0;
            const yVal =
              Math.abs(xi) < 1e-8
                ? -beta * Math.log(tailProb)
                : (beta / xi) * (Math.pow(tailProb, -xi) - 1.0);
            return 6.0 * v5 * yVal;
          },
          0.0,
          1.0,
          1e-11,
          16,
        );
      expect(Math.abs(esFormula - conditionalTailMean)).toBeLessThan(1e-8);
    }
  });

  it("recovers known GPD parameters via MLE and recovers 1/nu=0.2 on t(5) upper tail", () => {
    const rng = new Xoshiro256PlusPlus(123456);
    const nExceed = 10_000;
    const xiTrue = 0.25;
    const betaTrue = 0.7;
    const yGpd = new Float64Array(nExceed);
    for (let i = 0; i < nExceed; i += 1) {
      yGpd[i] = gpdQuantile(rng.nextUniformOpen(), xiTrue, betaTrue);
    }
    const fit = fitGpdMle(yGpd);
    expect(fit.converged).toBe(true);
    expect(Math.abs(fit.xi - xiTrue)).toBeLessThan(3.0 * fit.seXi);
    expect(Math.abs(fit.beta - betaTrue)).toBeLessThan(3.0 * fit.seBeta);

    // POT & Hill on t(5) standardized tail (N=100,000 fast test)
    const sample = new Float64Array(100_000);
    for (let i = 0; i < sample.length; i += 1) {
      sample[i] = stdTQuantile(rng.nextUniformOpen(), 5.0);
    }
    const pot = potFromSample(sample, 0.95, [0.99], [0.99]);
    expect(pot.gpd.xi).toBeGreaterThan(0.1);
    expect(pot.gpd.xi).toBeLessThan(0.32);

    const desc = Array.from(sample).sort((a, b) => b - a);
    const hill = hillEstimator(desc, 1000, true);
    expect(hill.xiHat).toBeGreaterThan(0.15);
    expect(hill.xiHat).toBeLessThan(0.35);
  });
});
