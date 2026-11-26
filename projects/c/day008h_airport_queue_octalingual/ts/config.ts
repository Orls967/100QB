import * as fs from "node:fs";

export const SLOTS_PER_DAY = 8640;
export const SLOTS_PER_HOUR = 360;
export const N_STN = 6;
export const MAX_SRV = 16;
export const MAX_Q = 4096;
export const MAX_POOL = 8192;
export const AUDIT_DAYS = 3;

export const STN_CHECKIN = 0;
export const STN_SECURITY = 1;
export const STN_IMIGRASI = 2;
export const STN_GATE = 3;
export const STN_BELT = 4;
export const STN_CUSTOMS = 5;

export const CLS_ECON = 0;
export const CLS_BIZ = 1;
export const CLS_TRANSIT = 2;

export const DEST_DOM = 0;
export const DEST_INTL = 1;

export const SRV_IDLE = 0;
export const SRV_BUSY = 1;
export const SRV_DOWN = 2;

export interface Config {
  kv: Record<string, number>;
  srvCnt: number[];
  muBp: number[];
}

export interface RunParams {
  days: number;
  arrScaleBp: number;
  breakProbBp: number;
  repairProbBp: number;
  beltCap: number;
  gateCloseLead: number;
  srvOvr: number[];
  muOvr: number[];
  priorityOn: number;
  breakdownOn: number;
  rosterOn: number;
  customsOn: number;
  permuteSrv: number;
  timescaleDiv2: number;
  auditOn: number;
}

export function defaultRp(c: Config, days: number): RunParams {
  return {
    days,
    arrScaleBp: c.kv["arr_scale_bp"],
    breakProbBp: c.kv["break_prob_bp"],
    repairProbBp: c.kv["repair_prob_bp"],
    beltCap: c.kv["belt_cap"],
    gateCloseLead: c.kv["gate_close_lead"],
    srvOvr: Array(N_STN).fill(-1),
    muOvr: Array(N_STN).fill(-1),
    priorityOn: 1,
    breakdownOn: 1,
    rosterOn: 1,
    customsOn: 0,
    permuteSrv: 0,
    timescaleDiv2: 0,
    auditOn: 0,
  };
}

export function cloneConfig(c: Config, seedOff: number): Config {
  return {
    kv: { ...c.kv, master_seed: c.kv["master_seed"] + seedOff },
    srvCnt: [...c.srvCnt],
    muBp: [...c.muBp],
  };
}

export function loadConfig(path: string): Config {
  const kv: Record<string, number> = {};
  const lines = fs.readFileSync(path, "utf8").split(/\r?\n/);
  for (const line of lines) {
    const idx = line.indexOf("=");
    if (idx < 0) continue;
    const k = line.slice(0, idx).trim();
    const v = Number.parseInt(line.slice(idx + 1).trim(), 10);
    kv[k] = v;
  }
  return {
    kv,
    srvCnt: [
      kv["checkin_servers"],
      kv["security_servers"],
      kv["imigrasi_servers"],
      kv["gate_servers"],
      kv["belt_servers"],
      kv["customs_servers"],
    ],
    muBp: [
      kv["checkin_mu_bp"],
      kv["security_mu_bp"],
      kv["imigrasi_mu_bp"],
      kv["gate_mu_bp"],
      kv["belt_mu_bp"],
      kv["customs_mu_bp"],
    ],
  };
}
