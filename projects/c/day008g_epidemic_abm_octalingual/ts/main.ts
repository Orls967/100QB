import { mkdirSync, readdirSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { writeS01ToS09, writeS10ToS37 } from "./analysis.ts";
import { simRun } from "./engine.ts";
import { Sim, defaultRP, loadConfig } from "./population.ts";
import { fnv1aBytes, hex64 } from "./rng.ts";
import { writeS38SVGs } from "./svg.ts";

function manifestScan(dirPath: string, prefix: string, lines: string[]): void {
  for (const nm of readdirSync(dirPath).sort()) {
    if (nm.startsWith(".")) continue;
    const p = join(dirPath, nm);
    if (statSync(p).isDirectory()) {
      manifestScan(p, `${prefix}${nm}/`, lines);
    } else if (nm !== "files_manifest.csv") {
      const b = readFileSync(p);
      let lc = 0;
      for (let i = 0; i < b.length; i++) if (b[i] === 10) lc++;
      lines.push(`${prefix}${nm},${lc},${b.length},${hex64(fnv1aBytes(b))}`);
    }
  }
}

function nowSec(): number {
  return Number(process.hrtime.bigint()) / 1e9;
}

const outDir = "out/ts";
mkdirSync(`${outDir}/figures`, { recursive: true });
const cfg = loadConfig("config.kv");
if (process.argv[2] === "s01") {
  writeS01ToS09(cfg, outDir);
  process.exit(0);
}

let tot = 0;
for (let st = 1; st <= 9; st++) console.log(`[ts] tahap ${st}/39 (s${String(st).padStart(2, "0")})`);
let t0 = nowSec();
writeS01ToS09(cfg, outDir);
tot += nowSec() - t0;

const nMain = 20000, betaCal = 400;
console.log("[ts] tahap 10/39 (s10)");
t0 = nowSec();
const sNo = new Sim(cfg, defaultRP(nMain, cfg.main_days, betaCal, 0), true);
simRun(sNo);
tot += nowSec() - t0;

console.log("[ts] tahap 11/39 (s11)");
t0 = nowSec();
const sMain = new Sim(cfg, defaultRP(nMain, cfg.main_days, betaCal, 1), true);
simRun(sMain);
tot += nowSec() - t0;

for (let st = 12; st <= 34; st++) console.log(`[ts] tahap ${st}/39 (s${String(st).padStart(2, "0")})`);
t0 = nowSec();
writeS10ToS37(sNo, sMain, cfg, outDir, betaCal);
tot += nowSec() - t0;

console.log("[ts] tahap 35/39 (s35)");
const vals: number[] = [];
for (let i = 0; i < 5; i++) {
  const sb = new Sim(cfg, defaultRP(2000, 60, betaCal, 1), false);
  const tb = nowSec();
  simRun(sb);
  vals.push(120000.0 / Math.max(nowSec() - tb, 1e-9));
}
vals.sort((a, b) => a - b);
const med = vals[2];
const dev = vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b);
writeFileSync(`${outDir}/s35_bench.csv`, `person_days,pd_per_sec_med,pd_per_sec_mad\n120000,${med.toFixed(2)},${dev[2].toFixed(2)}\n`);

console.log("[ts] tahap 36/39 (s36)\n[ts] tahap 37/39 (s37)\n[ts] tahap 38/39 (s38)");
t0 = nowSec();
writeS38SVGs(sNo, sMain, outDir);
tot += nowSec() - t0;

console.log("[ts] tahap 39/39 (s39)");
writeFileSync(`${outDir}/s39_summary.json`, `{"lang":"ts","N":${nMain},"days":${cfg.main_days},"N_sw":1000,"toolchain":"node ${process.version}","total_sec":${tot.toFixed(3)}}\n`);
const mLines = ["file,lines,bytes,fnv1a"];
manifestScan(outDir, "", mLines);
writeFileSync(`${outDir}/files_manifest.csv`, mLines.join("\n") + "\n");
