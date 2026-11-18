import { existsSync, mkdirSync, writeFileSync } from "node:fs";
import path from "node:path";
import { performance } from "node:perf_hooks";
import {
  type OptionParams,
  bsCall,
  deltaAsianBump,
  deltaAsianPathwise,
  fitLogLogSlope,
  geometricAsianCall,
  priceAsianAntithetic,
  priceAsianControlVariate,
  priceAsianNaive,
  priceEuropeanMC,
} from "./options.js";
import {
  type ConvergencePoint,
  type EstimatorSummary,
  renderConvergenceSvg,
  renderEstimatorsSvg,
} from "./svg.js";

function findRepoRoot(): string {
  let dir = process.cwd();
  while (true) {
    if (existsSync(path.join(dir, "src", "qb100"))) {
      return dir;
    }
    const parent = path.dirname(dir);
    if (parent === dir) {
      throw new Error("Repository root not found");
    }
    dir = parent;
  }
}

async function main(): Promise<void> {
  const params: OptionParams = {
    s0: 100.0,
    k: 100.0,
    r: 0.05,
    sigma: 0.2,
    t: 1.0,
    n_steps: 252,
  };
  const nPaths = 200_000;
  const seed = 42;

  const tTotalStart = performance.now();

  let t0 = performance.now();
  const bsEuro = bsCall(params);
  const mcEuro = await priceEuropeanMC(params, nPaths, seed);
  const tSanity = performance.now() - t0;
  const diffEuro = mcEuro.price - bsEuro;
  const zEuro = Math.abs(diffEuro) / mcEuro.se;

  t0 = performance.now();
  const resNaive = await priceAsianNaive(params, nPaths, seed);
  const tNaive = performance.now() - t0;

  t0 = performance.now();
  const resAnti = await priceAsianAntithetic(params, nPaths, seed);
  const tAnti = performance.now() - t0;

  t0 = performance.now();
  const resCV = await priceAsianControlVariate(params, nPaths, seed);
  const tCV = performance.now() - t0;

  t0 = performance.now();
  const deltaPW = await deltaAsianPathwise(params, nPaths, seed);
  const tPW = performance.now() - t0;

  t0 = performance.now();
  const deltaBmp = await deltaAsianBump(params, nPaths, 0.01, seed);
  const tBmp = performance.now() - t0;

  const coreMs = performance.now() - tTotalStart;

  const deltaDiff = deltaPW.delta - deltaBmp.delta;
  const deltaZ = Math.abs(deltaDiff) / Math.hypot(deltaPW.se, deltaBmp.se);

  const varNaive = resNaive.se * resNaive.se;
  const estimators: Array<EstimatorSummary & { beta?: number; y_exact?: number }> = [
    {
      name: "Naive MC (200k)",
      price: resNaive.price,
      se: resNaive.se,
      ci_low: resNaive.ci_low,
      ci_high: resNaive.ci_high,
      vrf: 1.0,
    },
    {
      name: "Antithetic (100k pasang)",
      price: resAnti.price,
      se: resAnti.se,
      ci_low: resAnti.ci_low,
      ci_high: resAnti.ci_high,
      vrf: varNaive / (resAnti.se * resAnti.se),
    },
    {
      name: "Control Variate (Geom. Asia)",
      price: resCV.price,
      se: resCV.se,
      ci_low: resCV.ci_low,
      ci_high: resCV.ci_high,
      vrf: varNaive / (resCV.se * resCV.se),
      beta: resCV.beta,
      y_exact: resCV.y_exact,
    },
  ];

  t0 = performance.now();
  const nGrid = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000];
  const convPoints: ConvergencePoint[] = [];
  const seN: number[] = [];
  const seA: number[] = [];
  const seC: number[] = [];
  for (const n of nGrid) {
    const [rn, ra, rc] = await Promise.all([
      priceAsianNaive(params, n, seed),
      priceAsianAntithetic(params, n, seed),
      priceAsianControlVariate(params, n, seed),
    ]);
    seN.push(rn.se);
    seA.push(ra.se);
    seC.push(rc.se);
    convPoints.push({
      n,
      se_naive: rn.se,
      se_antithetic: ra.se,
      se_control_variate: rc.se,
    });
  }
  const tConv = performance.now() - t0;

  const slopes = {
    naive: fitLogLogSlope(nGrid, seN),
    antithetic: fitLogLogSlope(nGrid, seA),
    control_variate: fitLogLogSlope(nGrid, seC),
  };

  const results = {
    params,
    n_paths: nPaths,
    seed,
    sanity_check: {
      bs_call: bsEuro,
      mc_call: mcEuro,
      diff: diffEuro,
      z_score: zEuro,
    },
    estimators,
    delta: {
      pathwise: deltaPW,
      bump_1pct: deltaBmp,
      diff: deltaDiff,
      z_score: deltaZ,
    },
    convergence: convPoints,
    slopes,
    timings_ms: {
      sanity_check: tSanity,
      naive: tNaive,
      antithetic: tAnti,
      control_variate: tCV,
      delta_pathwise: tPW,
      delta_bump: tBmp,
      convergence: tConv,
      total_200k_core: coreMs,
    },
  };

  const root = findRepoRoot();
  const outDir = path.join(root, "outputs", "day001f");
  const figDir = path.join(root, "outputs", "figures");
  mkdirSync(outDir, { recursive: true });
  mkdirSync(figDir, { recursive: true });

  writeFileSync(path.join(outDir, "results.json"), JSON.stringify(results, null, 2));
  writeFileSync(
    path.join(figDir, "day001f_asian-convergence-loglog.svg"),
    renderConvergenceSvg(convPoints, slopes),
  );
  writeFileSync(
    path.join(figDir, "day001f_asian-estimators-ci.svg"),
    renderEstimatorsSvg(estimators),
  );

  console.log(
    `BS Eropa: ${bsEuro.toFixed(6)} | MC Eropa: ${mcEuro.price.toFixed(6)} ± ${mcEuro.se.toFixed(6)} (${zEuro.toFixed(2)} SE)`,
  );
  console.log(
    `Geom. Asian exact: ${geometricAsianCall(params).toFixed(6)} | beta: ${resCV.beta.toFixed(6)}`,
  );
  for (const e of estimators) {
    console.log(
      `${e.name.padEnd(28)}: ${e.price.toFixed(6)} ± ${e.se.toFixed(6)} (CI [${e.ci_low.toFixed(6)}, ${e.ci_high.toFixed(6)}]) | VRF = ${e.vrf.toFixed(2)}x`,
    );
  }
  console.log(
    `Delta PW: ${deltaPW.delta.toFixed(6)} ± ${deltaPW.se.toFixed(6)} | Bump: ${deltaBmp.delta.toFixed(6)} ± ${deltaBmp.se.toFixed(6)} | diff = ${deltaDiff >= 0 ? "+" : ""}${deltaDiff.toFixed(6)} (${deltaZ.toFixed(2)} SE)`,
  );
  console.log(
    `Slopes -> Naive: ${slopes.naive.toFixed(4)}, Antithetic: ${slopes.antithetic.toFixed(4)}, CV: ${slopes.control_variate.toFixed(4)}`,
  );
  console.log(`Total 200k core time: ${coreMs.toFixed(1)} ms`);
}

await main();
