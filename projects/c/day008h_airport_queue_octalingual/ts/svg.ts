import * as fs from "node:fs";
import { STN_SECURITY } from "./config.ts";
import { BR_BIZ_PRIO, BR_SEC_EXTRA, BR_SRV_BREAK, BR_SRV_REPAIR, Sim } from "./network.ts";

function writeBarSvg(outDir: string, fname: string, title: string, vals: number[], col: string): void {
  const mx = Math.max(1.0, ...vals);
  const bw = 460.0 / Math.max(vals.length, 1);
  let s = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>';
  for (let i = 0; i < vals.length; i++) {
    const h = Math.max((vals[i] / mx) * 250.0, 1.0);
    s += `<rect x="${(75.0 + i * bw).toFixed(1)}" y="${(310.0 - h).toFixed(1)}" width="${(bw * 0.72).toFixed(1)}" height="${h.toFixed(1)}" fill="${col}"/>`;
  }
  s += `<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">${title}</text></svg>\n`;
  fs.writeFileSync(`${outDir}/figures/${fname}`, s);
}

export function writeS46Svgs(sMain: Sim, outDir: string): void {
  const qstn = Array.from({ length: 6 }, (_, i) => sMain.stn[i].arrCnt);
  const utilh = Array.from({ length: 12 }, (_, i) => sMain.hrUtilSum[i * 2][STN_SECURITY]);
  const wcls = Array.from({ length: 3 }, (_, i) => (sMain.clsWaitCnt[i] > 0 ? sMain.clsWaitSum[i] / sMain.clsWaitCnt[i] : 0.0));
  const soj = Array.from({ length: 3 }, (_, i) => (sMain.clsSojournCnt[i] > 0 ? sMain.clsSojournSum[i] / sMain.clsSojournCnt[i] : 0.0));
  const br = Array.from({ length: 6 }, (_, i) => sMain.branches[i]);
  writeBarSvg(outDir, "queue_stations.svg", "TypeScript Arrivals by Station", qstn, "#1f77b4");
  writeBarSvg(outDir, "util_hourly.svg", "TypeScript Security Busy Slots (2h sample)", utilh, "#ff7f0e");
  writeBarSvg(outDir, "wait_hist.svg", "TypeScript Mean Wait by Station", qstn, "#2ca02c");
  writeBarSvg(outDir, "wait_quantiles.svg", "TypeScript Mean Wait by Class (Econ/Biz/Trans)", wcls, "#d62728");
  writeBarSvg(outDir, "dtmc_validation.svg", "TypeScript Departures by Station", br, "#9467bd");
  writeBarSvg(outDir, "server_heatmap.svg", "TypeScript Arrivals vs Done vs Lost", [sMain.totalArr, sMain.totalDone, sMain.totalLost], "#8c564b");
  writeBarSvg(outDir, "priority_effect.svg", "TypeScript Priority Jumps &amp; Extra Security", [sMain.branches[BR_BIZ_PRIO], sMain.branches[BR_SEC_EXTRA]], "#e377c2");
  writeBarSvg(outDir, "breakdown_effect.svg", "TypeScript Breakdowns &amp; Repairs", [sMain.branches[BR_SRV_BREAK], sMain.branches[BR_SRV_REPAIR]], "#7f7f7f");
  writeBarSvg(outDir, "bottleneck.svg", "TypeScript Station Busy Slots", [sMain.stn[0].busySlots, sMain.stn[1].busySlots, sMain.stn[2].busySlots, sMain.stn[3].busySlots], "#bcbd22");
  writeBarSvg(outDir, "ensemble_band.svg", "TypeScript Abandonments vs Missed Flights", [sMain.totalAban, sMain.totalMissed], "#17becf");
  writeBarSvg(outDir, "crn_ratio.svg", "TypeScript Baggage Total vs Missed vs Blocked", [sMain.totalBags, sMain.totalBagsMiss, sMain.totalBeltBlock], "#1f77b4");
  writeBarSvg(outDir, "tornado.svg", "TypeScript Class Sojourn Means", soj, "#ff7f0e");
}
