import { mkdirSync, readdirSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { writeS01ToS09, writeS11ToS38 } from "./analysis.ts";
import { simRunDays } from "./engine.ts";
import { loadConfig, Sim } from "./entities.ts";
import { fnv1aBytes } from "./rng.ts";
import { writeS39Svgs } from "./svg.ts";

function manifestScan(dir: string, prefix: string, lines: string[]): void {
  for (const nm of readdirSync(dir).sort()) {
    if (nm.startsWith(".")) continue;
    const p = `${dir}/${nm}`;
    if (statSync(p).isDirectory()) manifestScan(p, `${prefix}${nm}/`, lines);
    else if (nm !== "files_manifest.csv") {
      const b = readFileSync(p);
      let lc = 0; for (let i = 0; i < b.length; i++) if (b[i] === 10) lc++;
      lines.push(`${prefix}${nm},${lc},${b.length},${fnv1aBytes(b).toString(16).padStart(16, "0")}`);
    }
  }
}

const outDir = "out/ts";
mkdirSync(`${outDir}/figures`, { recursive: true });
const cfg = loadConfig("config.kv");
const cmd = process.argv[2] ?? "all";
if (cmd === "s01") { writeS01ToS09(cfg, outDir); process.exit(0); }

let tot = 0.0;
for (let st = 1; st <= 9; st++) console.log(`[ts] tahap ${st}/40 (s${String(st).padStart(2, "0")})`);
let t0 = performance.now(); writeS01ToS09(cfg, outDir); tot += (performance.now() - t0) / 1000;
console.log("[ts] tahap 10/40 (s10)");
t0 = performance.now(); const sim = new Sim(cfg, true, outDir, 80000);
simRunDays(sim, cfg.num_days, 80000, 0, -1, 1); sim.closeLogs(); tot += (performance.now() - t0) / 1000;
for (let st = 11; st <= 38; st++) console.log(`[ts] tahap ${st}/40 (s${String(st).padStart(2, "0")})`);
t0 = performance.now(); writeS11ToS38(sim, cfg, outDir); tot += (performance.now() - t0) / 1000;
console.log("[ts] tahap 39/40 (s39)");
t0 = performance.now(); writeS39Svgs(sim, outDir); tot += (performance.now() - t0) / 1000;
console.log("[ts] tahap 40/40 (s40)");
writeFileSync(`${outDir}/s40_summary.json`, `{"lang":"ts","days":5,"ev_per_day":80000,"ev_total":${cfg.num_days * 80000},"toolchain":"node v24.21.0 (--experimental-strip-types)","total_sec":${tot.toFixed(3)}}\n`);
const mf = ["file,lines,bytes,fnv1a"]; manifestScan(outDir, "", mf);
writeFileSync(`${outDir}/files_manifest.csv`, mf.join("\n") + "\n");
