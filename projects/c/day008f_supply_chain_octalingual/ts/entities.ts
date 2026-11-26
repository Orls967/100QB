import { openSync, readFileSync, writeSync, closeSync } from "node:fs";
import { fnv1aU32Pair, MinHeap, splitmix64, Xoshiro256 } from "./rng.ts";
import { ACC_CUST, ACC_FEE, ACC_SUPP, NUM_ACCOUNTS, NUM_STORES } from "./inventory.ts";

export type Config = Record<string, number>;
export interface Account { cash: number; onHand: number[]; inTransit: number[]; backorder: number[]; creditLimit: number; status: number; ewmaMilli: number[]; plantState: number; disruptRem: number; }

export function loadConfig(path: string): Config {
  const c: Config = {};
  for (const ln of readFileSync(path, "utf8").split("\n")) {
    const eq = ln.indexOf("=");
    if (eq > 0) c[ln.slice(0, eq).trim()] = Number(ln.slice(eq + 1).trim());
  }
  return c;
}

export class Sim {
  cfg: Config; logEnabled: boolean; nextSeq = 1; eventHash = 0xcbf29ce484222325n;
  heap = new MinHeap(); rng: Xoshiro256[] = []; acc: Account[] = [];
  initCash = 0; initStock = [0, 0, 0, 0]; totalProduced = [0, 0, 0, 0]; totalSold = [0, 0, 0, 0]; totalKits = 0;
  branches = new Array(16).fill(0); carrierDayVol = [0, 0, 0]; carrierTierCnt = [0, 0, 0];
  allocTotalReq = 0; allocTotalFill = 0; allocCases = 0; allocNonzero = 0; invViolations = 0;
  dailyDemand = Array.from({ length: 6 }, () => [0, 0, 0, 0]);
  dailyCostEchelon = Array.from({ length: 6 }, () => [0, 0, 0, 0]);
  dailyFillStore = Array.from({ length: 10 }, () => [0, 0]);
  shockCnt = [0, 0, 0]; shockResp = [0, 0, 0];
  dayInvSnap = Array.from({ length: 6 }, () => [0, 0, 0, 0]);
  kitImbSeries: number[] = []; leadtimePairs: [number, number][] = [];
  sales: number[][] = []; ships: number[][] = []; prods: number[][] = []; kits: number[][] = []; fees: number[][] = [];
  fds: Record<string, number> = {};

  constructor(c: Config, logEn: boolean, outDir: string) {
    this.cfg = c; this.logEnabled = logEn;
    for (let i = 0; i < NUM_ACCOUNTS + 2; i++) {
      const x = BigInt(c.master_seed) ^ fnv1aU32Pair(Math.floor(i / 10), i);
      const [sd] = splitmix64(x);
      this.rng.push(new Xoshiro256(sd));
    }
    for (let i = 0; i < NUM_ACCOUNTS; i++)
      this.acc.push({ cash: 0, onHand: [0,0,0,0], inTransit: [0,0,0,0], backorder: [0,0,0,0], creditLimit: 0, status: 0, ewmaMilli: [0,0,0,0], plantState: 0, disruptRem: 0 });
    for (let i = 0; i < NUM_STORES; i++) { this.acc[i].cash = 2000000000; this.acc[i].creditLimit = 500000000; this.acc[i].onHand = [400,400,400,120]; this.acc[i].ewmaMilli = [12000,12000,12000,12000]; }
    for (let i = 10; i < 14; i++) { this.acc[i].cash = 5000000000; this.acc[i].creditLimit = 1000000000; this.acc[i].onHand = [2500,2500,2500,800]; }
    for (let i = 14; i < 17; i++) { this.acc[i].cash = 5000000000; this.acc[i].creditLimit = 2000000000; this.acc[i].onHand[i - 14] = 4000; }
    for (let i = 17; i < 20; i++) { this.acc[i].cash = 500000000; this.acc[i].creditLimit = 200000000; }
    this.acc[ACC_CUST].cash = 50000000000; this.acc[ACC_SUPP].cash = 50000000000; this.acc[ACC_FEE].cash = 1000000000;
    for (const a of this.acc) { this.initCash += a.cash; for (let k = 0; k < 4; k++) this.initStock[k] += a.onHand[k]; }
    if (logEn && outDir) {
      const hdrs: Record<string, string> = {
        "orders_log.csv": "day,tick,seq,etype,src,dst,sku,qty", "allocations.csv": "day,wh,sku,store,req,fill",
        "disruptions.csv": "day,tick,fac,state,dur", "credit_calls.csv": "day,acc,cash,req_margin,status",
        "risk_actions.csv": "day,tick,acc,action,amount", "transshipments.csv": "day,tick,src,dst,sku,qty",
        "shocks_log.csv": "day,tick,cls,store,mag", "accounts_series.csv": "day,acc,cash,oh_tot,it_tot,status",
      };
      for (const [fn, h] of Object.entries(hdrs)) {
        const fd = openSync(`${outDir}/${fn}`, "w"); writeSync(fd, h + "\n"); this.fds[fn] = fd;
      }
    }
  }
  logLine(fn: string, line: string): void { if (this.fds[fn] !== undefined) writeSync(this.fds[fn], line + "\n"); }
  closeLogs(): void { for (const fd of Object.values(this.fds)) closeSync(fd); this.fds = {}; }
}
