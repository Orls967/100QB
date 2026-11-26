import { CLS_BIZ, MAX_Q, MAX_SRV } from "./config.ts";
import type { Agent } from "./passenger.ts";

export interface Server {
  status: number;
  enabled: number;
  agentIdx: number;
  extraSec: number;
}

export class Station {
  q: Int32Array;
  head = 0;
  tail = 0;
  len = 0;
  cap: number;
  hasPrio: number;
  cMax: number;
  cActive: number;
  busyCnt = 0;
  srv: Server[];
  arrCnt = 0;
  depCnt = 0;
  blockCnt = 0;
  abanCnt = 0;
  sumL = 0;
  sumW = 0;
  waitCnt = 0;
  busySlots = 0;

  constructor(cMax: number, cap: number, hasPrio: number) {
    this.q = new Int32Array(MAX_Q);
    this.cap = cap;
    this.hasPrio = hasPrio;
    this.cMax = cMax;
    this.cActive = cMax;
    this.srv = Array.from({ length: MAX_SRV }, (_, i) => ({
      status: 0,
      enabled: i < cMax ? 1 : 0,
      agentIdx: -1,
      extraSec: 0,
    }));
  }

  pushTail(idx: number): void {
    if (this.len >= MAX_Q) return;
    this.q[this.tail] = idx;
    this.tail = (this.tail + 1) % MAX_Q;
    this.len++;
  }

  pushHead(idx: number): void {
    if (this.len >= MAX_Q) return;
    this.head = (this.head - 1 + MAX_Q) % MAX_Q;
    this.q[this.head] = idx;
    this.len++;
  }

  removeAt(pos: number): number {
    const idx = this.q[(this.head + pos) % MAX_Q];
    for (let i = pos; i < this.len - 1; i++) {
      this.q[(this.head + i) % MAX_Q] = this.q[(this.head + i + 1) % MAX_Q];
    }
    this.tail = (this.tail - 1 + MAX_Q) % MAX_Q;
    this.len--;
    return idx;
  }

  popNext(pool: Agent[], priorityOn: number, prioBr: number): [number, number] {
    if (this.len === 0) return [-1, prioBr];
    if (priorityOn !== 0 && this.hasPrio !== 0) {
      for (let i = 0; i < this.len; i++) {
        if (pool[this.q[(this.head + i) % MAX_Q]].cls === CLS_BIZ) {
          if (i > 0) return [this.removeAt(i), prioBr + 1];
          break;
        }
      }
    }
    const res = this.q[this.head];
    this.head = (this.head + 1) % MAX_Q;
    this.len--;
    return [res, prioBr];
  }
}
