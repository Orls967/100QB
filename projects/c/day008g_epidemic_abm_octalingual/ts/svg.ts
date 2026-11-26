import { writeFileSync } from "node:fs";
import { ST_D, ST_IA, ST_IS, ST_R, Sim } from "./population.ts";

function writeBarSVG(outDir: string, fname: string, title: string, vals: number[], col: string): void {
  const mx = Math.max(1.0, ...vals);
  const bw = 460.0 / Math.max(vals.length, 1);
  let out = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>';
  for (let i = 0; i < vals.length; i++) {
    const h = Math.max((vals[i] / mx) * 250.0, 1.0);
    out += `<rect x="${(75.0 + i * bw).toFixed(1)}" y="${(310.0 - h).toFixed(1)}" width="${(bw * 0.72).toFixed(1)}" height="${h.toFixed(1)}" fill="${col}"/>`;
  }
  out += `<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">${title}</text></svg>\n`;
  writeFileSync(`${outDir}/figures/${fname}`, out);
}

export function writeS38SVGs(sNo: Sim, sMain: Sim, outDir: string): void {
  const dCnt = sMain.rp.days;
  const step = Math.max(Math.floor(dCnt / 12), 1);
  const epi: number[] = [], inc: number[] = [], rt: number[] = [], hosp: number[] = [], ens: number[] = [];
  for (let i = 0; i < 12; i++) {
    const d = Math.min(i * step, dCnt - 1);
    epi.push(sMain.daily_cnt[d][ST_IA] + sMain.daily_cnt[d][ST_IS]);
    inc.push(sMain.daily_inc_age[d][0] + sMain.daily_inc_age[d][1] + sMain.daily_inc_age[d][2]);
    rt.push(sNo.daily_cnt[d][ST_IA] + sNo.daily_cnt[d][ST_IS]);
    hosp.push(sMain.daily_hosp[d][0]);
    ens.push(sMain.daily_var[d][1]);
  }
  writeBarSVG(outDir, "epi_curves.svg", "TypeScript Active Infectious (Interv)", epi, "#1f77b4");
  writeBarSVG(outDir, "incidence_age.svg", "TypeScript Daily Incidence Sample", inc, "#ff7f0e");
  writeBarSVG(outDir, "rt_series.svg", "TypeScript Active Infectious (No Interv)", rt, "#d62728");
  writeBarSVG(outDir, "hospital_cap.svg", "TypeScript Hospital Occupancy Sample", hosp, "#9467bd");
  writeBarSVG(outDir, "layer_attack.svg", "TypeScript Infections by Layer (HH/SW/Comm)", sMain.layer_inf, "#2ca02c");
  writeBarSVG(outDir, "offspring_dist.svg", "TypeScript Infections by Age Group", sMain.age_inf, "#8c564b");
  writeBarSVG(outDir, "sweep_heatmap.svg", "TypeScript No-Interv vs Interv Total Infections", [sNo.total_trans, sMain.total_trans], "#e377c2");
  writeBarSVG(outDir, "vaccine_effect.svg", "TypeScript Vaccine Doses &amp; Blocked Transmissions", [sMain.branches[18], sMain.branches[3]], "#17becf");
  writeBarSVG(outDir, "variant_share.svg", "TypeScript Variant 1 Daily Incidence", ens, "#bcbd22");
  writeBarSVG(outDir, "ensemble_band.svg", "TypeScript Final Recovered vs Deaths", [sMain.daily_cnt[dCnt - 1][ST_R], sMain.daily_cnt[dCnt - 1][ST_D]], "#7f7f7f");
}
