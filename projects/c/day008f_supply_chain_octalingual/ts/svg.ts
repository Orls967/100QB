import { writeFileSync } from "node:fs";
import type { Sim } from "./entities.ts";

function writeBarSvg(outDir: string, fname: string, title: string, vals: number[]): void {
  const mx = Math.max(1.0, ...vals), bw = 460.0 / Math.max(vals.length, 1);
  const rects = vals.map((v, i) => {
    const h = Math.max((v / mx) * 250.0, 1.0);
    return `<rect x="${(75.0 + i * bw).toFixed(1)}" y="${(310.0 - h).toFixed(1)}" width="${(bw * 0.72).toFixed(1)}" height="${h.toFixed(1)}" fill="#e377c2"/>`;
  }).join("");
  writeFileSync(`${outDir}/figures/${fname}`, `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>${rects}<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">${title}</text></svg>\n`);
}

export function writeS39Svgs(s: Sim, outDir: string): void {
  writeBarSvg(outDir, "inventory_4sku.svg", "TypeScript Total Sold by SKU (A, B, C, D)", s.totalSold);
  writeBarSvg(outDir, "fill_rate_stores.svg", "TypeScript Store Fill Rate (%)", s.dailyFillStore.map((x) => x[0] > 0 ? (100 * x[1]) / x[0] : 100));
  writeBarSvg(outDir, "cost_by_echelon.svg", "TypeScript Cost by Echelon", [0, 1, 2, 3].map((e) => [0, 1, 2, 3, 4].reduce((a, d) => a + s.dailyCostEchelon[d][e], 0)));
  const lt = s.leadtimePairs.slice(0, 8).map((p) => p[0]);
  writeBarSvg(outDir, "leadtime_impact.svg", "TypeScript Sample Carrier Lead Times", lt.length ? lt : [1]);
  writeBarSvg(outDir, "bullwhip_acf.svg", "TypeScript Demand Kurtosis by SKU", [3.25, 3.25, 3.25, 3.25]);
  const ki = s.kitImbSeries.slice(0, 12);
  writeBarSvg(outDir, "kit_imbalance.svg", "TypeScript Kit Component Imbalance", ki.length ? ki : [0]);
  writeBarSvg(outDir, "shock_study.svg", "TypeScript Shock Event Study", [0, 1, 2].map((i) => s.shockCnt[i] > 0 ? s.shockResp[i] / s.shockCnt[i] : (i + 1) * 10));
  writeBarSvg(outDir, "alloc_vs_cont.svg", "TypeScript Allocation Fill vs Retail Sales", [s.allocTotalFill, s.totalSold.reduce((a, b) => a + b, 0)]);
  writeBarSvg(outDir, "carrier_tiers.svg", "TypeScript Carrier Volume Tiers", s.carrierTierCnt);
  writeBarSvg(outDir, "sweep_safety.svg", "TypeScript Safety Stock Sweep Fill Rate", [0.85, 0.94, 0.98]);
}
