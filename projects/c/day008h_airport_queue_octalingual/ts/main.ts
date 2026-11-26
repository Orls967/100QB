import * as fs from "node:fs";
import * as path from "node:path";
import { writeS14ToS45 } from "./analysis.ts";
import { defaultRp, loadConfig } from "./config.ts";
import { simRun } from "./engine.ts";
import { Sim } from "./network.ts";
import { fnv1aBytes, hex64 } from "./rng.ts";
import { writeS01ToS13 } from "./stats.ts";
import { writeS46Svgs } from "./svg.ts";

function manifestScan(dir: string, prefix: string, lines: string[]): void {
  const entries = fs.readdirSync(dir).sort();
  for (const nm of entries) {
    if (nm.startsWith(".")) continue;
    const p = path.join(dir, nm);
    if (fs.statSync(p).isDirectory()) {
      manifestScan(p, `${prefix}${nm}/`, lines);
    } else if (nm !== "files_manifest.csv") {
      const b = fs.readFileSync(p);
      let nl = 0;
      for (let i = 0; i < b.length; i++) if (b[i] === 10) nl++;
      lines.push(`${prefix}${nm},${nl},${b.length},${hex64(fnv1aBytes(b))}`);
    }
  }
}

function main(): void {
  const outDir = "out/ts";
  fs.mkdirSync(`${outDir}/figures`, { recursive: true });
  const cfg = loadConfig("config.kv");
  if (process.argv[2] === "s01") {
    writeS01ToS13(cfg, outDir);
    return;
  }

  let tot = 0.0;
  for (let st = 1; st <= 13; st++) console.log(`[ts] tahap ${st}/48 (s${String(st).padStart(2, "0")})`);
  let t0 = performance.now();
  writeS01ToS13(cfg, outDir);
  tot += (performance.now() - t0) / 1000.0;

  const mainDays = 60;
  console.log("[ts] tahap 14/48 (s14)");
  t0 = performance.now();
  const rpm = defaultRp(cfg, mainDays);
  rpm.auditOn = 1;
  const sMain = new Sim(cfg, rpm);
  simRun(sMain);
  tot += (performance.now() - t0) / 1000.0;

  for (let st = 15; st <= 41; st++) console.log(`[ts] tahap ${st}/48 (s${String(st).padStart(2, "0")})`);
  t0 = performance.now();
  writeS14ToS45(cfg, outDir, sMain);
  tot += (performance.now() - t0) / 1000.0;

  console.log("[ts] tahap 42/48 (s42)");
  const vals: number[] = [];
  for (let i = 0; i < 5; i++) {
    const sb = new Sim(cfg, defaultRp(cfg, 3));
    const tb = performance.now();
    simRun(sb);
    const dt = Math.max((performance.now() - tb) / 1000.0, 1e-9);
    vals.push(sb.sysSumL / dt);
  }
  vals.sort((a, b) => a - b);
  const med = vals[2];
  const dev = vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b);
  fs.writeFileSync(`${outDir}/s42_bench.csv`, `slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,${med.toFixed(2)},${dev[2].toFixed(2)}\n`);

  console.log("[ts] tahap 43/48 (s43)\n[ts] tahap 44/48 (s44)\n[ts] tahap 45/48 (s45)\n[ts] tahap 46/48 (s46)");
  t0 = performance.now();
  writeS46Svgs(sMain, outDir);
  tot += (performance.now() - t0) / 1000.0;

  console.log("[ts] tahap 47/48 (s47)\n[ts] tahap 48/48 (s48)");
  fs.writeFileSync(`${outDir}/s48_summary.json`, `{"lang":"ts","days":${mainDays},"sw_days":3,"toolchain":"node ${process.version} (--experimental-strip-types)","total_sec":${tot.toFixed(3)}}\n`);
  const mfLines = ["file,lines,bytes,fnv1a"];
  manifestScan(outDir, "", mfLines);
  fs.writeFileSync(`${outDir}/files_manifest.csv`, mfLines.join("\n") + "\n");
}

main();
