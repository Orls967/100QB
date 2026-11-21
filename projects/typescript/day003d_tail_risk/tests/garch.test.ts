import { describe, expect, it } from "vitest";
import {
  exactGarchSqAcf,
  fitGarchGauss,
  fitGarchT,
  garchTLogLikAndScore,
  reconstructSigma2,
  scoreVsCentralFd,
  simulateGarchT,
  theoreticalGarchMoments,
  TRUE_GARCH_PARAMS,
  type GarchTParams,
} from "../src/garch.js";
import { nelderMead } from "../src/optim.js";
import { Xoshiro256PlusPlus } from "../src/rng.js";
import { mean, sampleAcf, sampleExcessKurtosis, se, variance } from "../src/stats.js";

describe("GARCH(1,1)-t simulation, analytical score, reconstruction & MLE", () => {
  it("reconstructs sigma_t^2 identically to simulation up to 1e-12", () => {
    const rng = new Xoshiro256PlusPlus(20261121);
    const sim = simulateGarchT(2500, TRUE_GARCH_PARAMS, rng);
    const rec = reconstructSigma2(sim.returns, TRUE_GARCH_PARAMS);
    let maxDiff = 0.0;
    for (let t = 0; t < 2500; t += 1) {
      const d = Math.abs((rec.sigma2[t] ?? 0.0) - (sim.sigma2[t] ?? 0.0));
      if (d > maxDiff) maxDiff = d;
    }
    expect(maxDiff).toBeLessThan(1e-12);
    expect(Math.abs(rec.sigma2Next - sim.sigma2NextTrue)).toBeLessThan(1e-15);
  });

  it("verifies analytical GARCH-t score against central finite differences (< 1e-5)", () => {
    const rng = new Xoshiro256PlusPlus(777);
    const sim = simulateGarchT(1500, TRUE_GARCH_PARAMS, rng);
    const testPoints: GarchTParams[] = [
      TRUE_GARCH_PARAMS,
      { mu: -0.0002, omega: 3.5e-6, alpha: 0.11, beta: 0.85, nu: 6.2 },
    ];
    for (const pt of testPoints) {
      const check = scoreVsCentralFd(sim.returns, pt);
      expect(check.maxRelErr).toBeLessThan(1e-5);
    }
  });

  it("BFGS matches Nelder-Mead on GARCH-t log-likelihood and stays off parameter bounds", () => {
    const rng = new Xoshiro256PlusPlus(9999);
    const sim = simulateGarchT(2500, TRUE_GARCH_PARAMS, rng);
    const fit = fitGarchT(sim.returns, undefined, true);
    expect(fit.converged).toBe(true);
    expect(fit.hitBoundary).toBe(false);
    expect(fit.gradNormInf).toBeLessThan(1e-5);

    const nm = nelderMead(
      (theta) => {
        const p: GarchTParams = {
          mu: theta[0] ?? 0,
          omega: theta[1] ?? 1e-6,
          alpha: theta[2] ?? 0.05,
          beta: theta[3] ?? 0.9,
          nu: theta[4] ?? 5,
        };
        if (
          p.omega <= 0 ||
          p.alpha < 0 ||
          p.beta < 0 ||
          p.alpha + p.beta >= 0.9999 ||
          p.nu <= 2.01
        ) {
          return 1e12;
        }
        return -garchTLogLikAndScore(sim.returns, p).logLik / sim.returns.length;
      },
      [fit.params.mu, fit.params.omega, fit.params.alpha, fit.params.beta, fit.params.nu],
      400,
      1e-10,
    );
    expect(Math.abs(fit.logLik / 2500 + nm.value)).toBeLessThan(1e-7);

    // Also check Gaussian QMLE & warm-start equivalence
    const gFit = fitGarchGauss(sim.returns);
    expect(gFit.converged).toBe(true);
    expect(gFit.hitBoundary).toBe(false);

    const warmFit = fitGarchT(sim.returns, TRUE_GARCH_PARAMS, false);
    expect(Math.abs(warmFit.logLik - fit.logLik)).toBeLessThan(1e-5);
  });

  it("checks unconditional moments and proves non-existence vs existence of eps^4 ACF formula", () => {
    // Verify mathematical condition D < 0 for true model
    const theoTrue = theoreticalGarchMoments(TRUE_GARCH_PARAMS);
    expect(theoTrue.unconditionalVariance).toBeCloseTo(1e-4, 15);
    expect(theoTrue.innovationExcessKurtosisZ).toBeCloseTo(6.0, 12);
    expect(theoTrue.fourthMomentFinite).toBe(false); // D = -0.0116 < 0

    // Verify finite-4th-moment anchor: GARCH-t(12) with alpha=0.10, beta=0.80 => D = 0.1625 > 0
    const anchorParams: GarchTParams = { mu: 0.0003, omega: 1e-5, alpha: 0.1, beta: 0.8, nu: 12.0 };
    const theoAnchor = exactGarchSqAcf(
      anchorParams.alpha,
      anchorParams.beta,
      3.0 + 6.0 / (12.0 - 4.0),
    );
    expect(theoAnchor.D).toBeGreaterThan(0.15);

    const rng = new Xoshiro256PlusPlus(314159);
    const acf1s: number[] = [];
    const acf2s: number[] = [];
    const vars: number[] = [];
    const kurts: number[] = [];
    for (let rep = 0; rep < 30; rep += 1) {
      const sim = simulateGarchT(5000, anchorParams, rng);
      const eps = new Float64Array(5000);
      const eps2 = new Float64Array(5000);
      for (let t = 0; t < 5000; t += 1) {
        const e = (sim.returns[t] ?? 0.0) - anchorParams.mu;
        eps[t] = e;
        eps2[t] = e * e;
      }
      vars.push(variance(eps, 0));
      kurts.push(sampleExcessKurtosis(sim.z));
      acf1s.push(sampleAcf(eps2, 1));
      acf2s.push(sampleAcf(eps2, 2));
    }
    const trueVarAnchor = anchorParams.omega / (1 - anchorParams.alpha - anchorParams.beta);
    expect(Math.abs(mean(vars) - trueVarAnchor)).toBeLessThan(3.0 * se(vars));
    expect(Math.abs(mean(kurts) - 6.0 / (12.0 - 4.0))).toBeLessThan(3.0 * se(kurts));
    expect(Math.abs(mean(acf1s) - theoAnchor.rho1)).toBeLessThan(3.0 * se(acf1s));
    expect(Math.abs(mean(acf2s) - theoAnchor.rho2)).toBeLessThan(3.0 * se(acf2s));
  });
});
