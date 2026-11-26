import { readFileSync } from "node:fs";
import { FNV_OFFSET, MASK64, Xoshiro256, fnv1aU32Pair, splitmix64 } from "./rng.ts";

export const BP = 10000;
export const MAX_DAYS = 200;
export const ST_S = 0;
export const ST_E = 1;
export const ST_IA = 2;
export const ST_IS = 3;
export const ST_H = 4;
export const ST_ICU = 5;
export const ST_R = 6;
export const ST_D = 7;

export type Config = Record<string, number>;

export interface RunParams {
  n: number;
  days: number;
  beta_bp: number;
  compliance_ovr: number;
  test_ovr: number;
  vax_day_ovr: number;
  vax_cov_ovr: number;
  school_ovr: number;
  lock_ovr: number;
  var_day_ovr: number;
  mask_ovr: number;
  contact_scale_bp: number;
  vax_perfect: number;
  interv_on: number;
}

export function defaultRP(n: number, days: number, beta_bp: number, interv_on: number): RunParams {
  return {
    n, days, beta_bp,
    compliance_ovr: -1, test_ovr: -1, vax_day_ovr: -1, vax_cov_ovr: -1,
    school_ovr: -1, lock_ovr: -1, var_day_ovr: -1, mask_ovr: -1,
    contact_scale_bp: 10000, vax_perfect: 0, interv_on,
  };
}

export function loadConfig(path: string): Config {
  const c: Config = {};
  for (const line of readFileSync(path, "utf8").split("\n")) {
    const idx = line.indexOf("=");
    if (idx < 0) continue;
    c[line.slice(0, idx).trim()] = Number.parseInt(line.slice(idx + 1).trim(), 10);
  }
  return c;
}

export interface Person {
  age_grp: number;
  hh_id: number;
  sw_id: number;
  status: number;
  variant: number;
  state_day: number;
  dur_days: number;
  inf_day: number;
  vax_doses: number;
  vax_day: number;
  compliance: number;
  quar_until: number;
  super_flag: number;
  offspring: number;
  ever_inf: number;
  hosp_rej: number;
}

export class Sim {
  cfg: Config;
  rp: RunParams;
  n: number;
  hosp_cap: number;
  icu_cap: number;
  hosp_occ = 0;
  icu_occ = 0;
  pop: Person[];
  n_hh: number;
  n_sw: number;
  hh_off: Int32Array;
  hh_mem: Int32Array;
  sw_off: Int32Array;
  sw_mem: Int32Array;
  seen_day: Int32Array;
  rng: Xoshiro256[];
  cand_tgt: number[] = [];
  cand_src: number[] = [];
  cand_lay: number[] = [];
  cand_var: number[] = [];
  daily_cnt: number[][];
  daily_inc_age: number[][];
  daily_deaths_age: number[][];
  daily_hosp: number[][];
  daily_tests: number[][];
  daily_quar: number[];
  daily_vax: number[][];
  daily_var: number[][];
  layer_inf = [0, 0, 0];
  age_inf = [0, 0, 0];
  age_tot = [0, 0, 0];
  branches = new Array<number>(21).fill(0);
  inv_violations = 0;
  total_trans = 0;
  total_reinf = 0;
  event_hash = FNV_OFFSET;
  log_enabled: boolean;
  trans: [number, number, number, number, number, number][] = [];
  progs: [number, number, number, number][] = [];

  constructor(cfg: Config, rp: RunParams, logEn: boolean) {
    this.cfg = { ...cfg };
    this.rp = { ...rp };
    this.n = rp.n;
    this.hosp_cap = Math.max(Math.floor((this.n * cfg.hosp_beds_per_10k) / BP), 1);
    this.icu_cap = Math.max(Math.floor((this.n * cfg.icu_beds_per_10k) / BP), 1);
    this.log_enabled = logEn;

    this.rng = [];
    for (let i = 0; i < 5; i++) {
      const x = (BigInt(cfg.master_seed) & MASK64) ^ fnv1aU32Pair(100, i);
      const [, sd] = splitmix64(x);
      this.rng.push(new Xoshiro256(sd));
    }

    const n = this.n;
    this.pop = new Array<Person>(n);
    this.seen_day = new Int32Array(n).fill(-1);
    for (let i = 0; i < n; i++) {
      this.pop[i] = {
        age_grp: 0, hh_id: 0, sw_id: -1, status: 0, variant: 0,
        state_day: 0, dur_days: 0, inf_day: -1, vax_doses: 0, vax_day: -1,
        compliance: 0, quar_until: -1, super_flag: 0, offspring: 0, ever_inf: 0, hosp_rej: 0,
      };
    }

    const hhCnt: number[] = [];
    let pid = 0;
    let hh = 0;
    while (pid < n) {
      const sz = Math.min(this.rng[0].randInt(1, 6), n - pid);
      hhCnt.push(sz);
      for (let k = 0; k < sz; k++) {
        const p = this.pop[pid];
        p.hh_id = hh;
        const u = this.rng[0].randInt(1, 100);
        const ag = u <= 22 ? 0 : u <= 82 ? 1 : 2;
        p.age_grp = ag;
        this.age_tot[ag]++;
        p.compliance = this.rng[0].randInt(4000, 9500);
        p.super_flag = this.rng[0].randInt(1, 10) === 1 ? 1 : 0;
        pid++;
      }
      hh++;
    }
    this.n_hh = hh;
    this.hh_off = new Int32Array(hh + 1);
    for (let i = 0; i < hh; i++) this.hh_off[i + 1] = this.hh_off[i] + hhCnt[i];
    this.hh_mem = new Int32Array(n);
    const cur = new Int32Array(hh);
    for (let i = 0; i < n; i++) {
      const h = this.pop[i].hh_id;
      this.hh_mem[this.hh_off[h] + cur[h]++] = i;
    }

    const nSchools = Math.floor(this.age_tot[0] / 25) + 1;
    const nWorks = Math.floor(this.age_tot[1] / 15) + 1;
    this.n_sw = nSchools + nWorks;
    const swCnt = new Int32Array(this.n_sw);
    for (let i = 0; i < n; i++) {
      const ag = this.pop[i].age_grp;
      if (ag === 0) {
        const sw = this.rng[0].randInt(0, nSchools - 1);
        this.pop[i].sw_id = sw;
        swCnt[sw]++;
      } else if (ag === 1) {
        const sw = nSchools + this.rng[0].randInt(0, nWorks - 1);
        this.pop[i].sw_id = sw;
        swCnt[sw]++;
      }
    }
    this.sw_off = new Int32Array(this.n_sw + 1);
    for (let i = 0; i < this.n_sw; i++) this.sw_off[i + 1] = this.sw_off[i] + swCnt[i];
    this.sw_mem = new Int32Array(n);
    const swCur = new Int32Array(this.n_sw);
    for (let i = 0; i < n; i++) {
      const sw = this.pop[i].sw_id;
      if (sw >= 0) this.sw_mem[this.sw_off[sw] + swCur[sw]++] = i;
    }

    this.daily_cnt = Array.from({ length: MAX_DAYS }, () => new Array<number>(8).fill(0));
    this.daily_inc_age = Array.from({ length: MAX_DAYS }, () => [0, 0, 0]);
    this.daily_deaths_age = Array.from({ length: MAX_DAYS }, () => [0, 0, 0]);
    this.daily_hosp = Array.from({ length: MAX_DAYS }, () => [0, 0, 0, 0]);
    this.daily_tests = Array.from({ length: MAX_DAYS }, () => [0, 0]);
    this.daily_quar = new Array<number>(MAX_DAYS).fill(0);
    this.daily_vax = Array.from({ length: MAX_DAYS }, () => [0, 0]);
    this.daily_var = Array.from({ length: MAX_DAYS }, () => [0, 0]);
  }
}
