import { MAX_POOL } from "./config.ts";

export interface Agent {
  active: number;
  cls: number;
  dest: number;
  flightHr: number;
  nBags: number;
  curStn: number;
  extraSec: number;
  isBag: number;
  arrSlot: number;
  qEnterSlot: number;
  depSlot: number;
}

export function emptyAgent(active = 0): Agent {
  return {
    active,
    cls: 0,
    dest: 0,
    flightHr: 0,
    nBags: 0,
    curStn: 0,
    extraSec: 0,
    isBag: 0,
    arrSlot: 0,
    qEnterSlot: 0,
    depSlot: 0,
  };
}

export interface AgentPoolHolder {
  pool: Agent[];
  freeStk: number[];
}

export function poolAlloc(s: AgentPoolHolder): number {
  const idx = s.freeStk.pop();
  if (idx === undefined) return -1;
  s.pool[idx] = emptyAgent(1);
  return idx;
}

export function poolFree(s: AgentPoolHolder, idx: number): void {
  if (idx < 0 || idx >= MAX_POOL || s.pool[idx].active === 0) return;
  s.pool[idx].active = 0;
  s.freeStk.push(idx);
}
