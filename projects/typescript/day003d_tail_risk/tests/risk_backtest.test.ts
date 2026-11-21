import { describe, expect, it } from "vitest";
import {
  acerbiSzekelyZ2,
  christoffersenFromTransitions,
  christoffersenTest,
  kupiecExactAnalyticalSize,
  kupiecTest,
  simulateZ2Single,
} from "../src/backtest.js";
import { stdTEs, stdTQuantile } from "../src/distributions.js";
import { fitGarchT, simulateGarchT, TRUE_GARCH_PARAMS } from "../src/garch.js";
import { runParallelTask } from "../src/parallel.js";
import { estimateRiskMethodsOnSample } from "../src/risk.js";
import { Xoshiro256PlusPlus } from "../src/rng.js";
import { mean, quantile, se } from "../src/stats.js";

describe("Conditional VaR/ES, backtests, and deterministic parallelism", () => {
  it("computes loss VaR/ES consistently so true GARCH-t achieves exact 1% and 2.5% violation rate", () => {
    const rng = new Xoshiro256PlusPlus(20261121);
    const n = 20_000;
    const sim = simulateGarchT(n, TRUE_GARCH_PARAMS, rng);
    const q99 = stdTQuantile(0.99, TRUE_GARCH_PARAMS.nu);
    const q975 = stdTQuantile(0.975, TRUE_GARCH_PARAMS.nu);
    const es99 = stdTEs(0.99, TRUE_GARCH_PARAMS.nu);
    const es975 = stdTEs(0.975, TRUE_GARCH_PARAMS.nu);

    const losses = new Float64Array(n);
    const var99Arr = new Float64Array(n);
    const es99Arr = new Float64Array(n);
    const exc99 = new Uint8Array(n);
    let cnt975 = 0;

    for (let t = 0; t < n; t += 1) {
      const sig = Math.sqrt(sim.sigma2[t] ?? 1e-4);
      const loss = -(sim.returns[t] ?? 0.0);
      const v99 = -TRUE_GARCH_PARAMS.mu + sig * q99;
      const v975 = -TRUE_GARCH_PARAMS.mu + sig * q975;
      losses[t] = loss;
      var99Arr[t] = v99;
      es99Arr[t] = -TRUE_GARCH_PARAMS.mu + sig * es99;
      if (loss > v99) exc99[t] = 1;
      if (loss > v975) cnt975 += 1;
    }

    const rate99 = exc99.reduce((a, b) => a + b, 0) / n;
    const rate975 = cnt975 / n;
    const se99 = Math.sqrt((0.01 * 0.99) / n);
    const se975 = Math.sqrt((0.025 * 0.975) / n);
    expect(Math.abs(rate99 - 0.01)).toBeLessThan(3.0 * se99);
    expect(Math.abs(rate975 - 0.025)).toBeLessThan(3.0 * se975);

    const ku = kupiecTest(exc99, 0.99);
    expect(ku.pValueExact).toBeGreaterThan(0.01);
    const ch = christoffersenTest(exc99, 0.99);
    expect(ch.pCc).toBeGreaterThan(0.01);
    const z2 = acerbiSzekelyZ2(losses, var99Arr, es99Arr, 0.99);
    expect(Math.abs(z2)).toBeLessThan(0.25);
    void es975;
  });

  it("evaluates risk estimators on filtered GARCH-t residuals cleanly", () => {
    const rng = new Xoshiro256PlusPlus(4242);
    const sim = simulateGarchT(2500, TRUE_GARCH_PARAMS, rng);
    const fit = fitGarchT(sim.returns, undefined, false);
    const risk = estimateRiskMethodsOnSample(
      fit.params,
      fit.residualsHat,
      Math.sqrt(fit.sigma2NextHat),
      Math.sqrt(sim.sigma2NextTrue),
      TRUE_GARCH_PARAMS,
    );
    expect(risk.methods.gaussian.var999).toBeLessThan(risk.trueValues.var999);
    expect(risk.methods.fhs_evt.var999).toBeGreaterThan(risk.methods.gaussian.var999);
    expect(risk.thresholdSensitivity.length).toBe(3);
  });

  it("verifies Kupiec, Christoffersen power under Markov clustering, and Z2 null vs alternative", () => {
    // Kupiec analytical size vs 0.05 nominal
    const sz99 = kupiecExactAnalyticalSize(1000, 0.01, 0.05);
    expect(sz99.exactSize).toBeLessThanOrEqual(0.05);
    expect(sz99.exactSize).toBeGreaterThan(0.025);

    // Christoffersen rejects clustered transitions strongly
    const chCluster = christoffersenFromTransitions(970, 5, 5, 20, 25, 1000, 0.01);
    expect(chCluster.pInd).toBeLessThan(1e-10);
    expect(chCluster.pCc).toBeLessThan(1e-10);

    // Z2 under H0 has mean near 0, whereas Gaussian misspecification has negative Z2
    const rngNull = new Xoshiro256PlusPlus(555);
    const zNull: number[] = [];
    const zGauss: number[] = [];
    for (let i = 0; i < 400; i += 1) {
      zNull.push(simulateZ2Single(500, [0.975], rngNull, "true")[0] ?? 0.0);
      zGauss.push(simulateZ2Single(500, [0.975], rngNull, "gaussian")[0] ?? 0.0);
    }
    expect(Math.abs(mean(zNull))).toBeLessThan(3.0 * se(zNull));
    expect(mean(zGauss)).toBeLessThan(-0.12);
    const crit = quantile(zNull, 0.05);
    expect(crit).toBeLessThan(0.0);
  });

  it("guarantees bit-for-bit identical parallel execution across workers=1 vs workers=4", async () => {
    const seed = 20261121;
    const salt = 91;

    // 1) Pure numerical task (pbh_sample)
    const pbh1 = await runParallelTask({ kind: "pbh_sample", nu: 5 }, 5000, 1250, seed, salt, 1);
    const pbh4 = await runParallelTask({ kind: "pbh_sample", nu: 5 }, 5000, 1250, seed, salt, 4);
    expect(pbh1.workerCount).toBe(1);
    expect(pbh4.workerCount).toBe(4);
    for (let c = 0; c < pbh1.outputs.length; c += 1) {
      const o1 = pbh1.outputs[c];
      const o4 = pbh4.outputs[c];
      if (o1?.kind === "pbh_sample" && o4?.kind === "pbh_sample") {
        expect(Array.from(o1.samples)).toEqual(Array.from(o4.samples));
      }
    }

    // 2) Full GARCH simulation + BFGS fit task
    const g1 = await runParallelTask(
      { kind: "garch_fit_rep", nObs: 800, includeRisk: true },
      8,
      2,
      seed,
      salt + 1,
      1,
    );
    const g4 = await runParallelTask(
      { kind: "garch_fit_rep", nObs: 800, includeRisk: true },
      8,
      2,
      seed,
      salt + 1,
      4,
    );
    expect(g4.outputs).toEqual(g1.outputs);
  });
});
