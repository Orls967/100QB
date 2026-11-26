import {
  AUDIT_DAYS,
  CLS_BIZ,
  CLS_ECON,
  CLS_TRANSIT,
  DEST_DOM,
  DEST_INTL,
  MAX_POOL,
  MAX_Q,
  MAX_SRV,
  N_STN,
  SLOTS_PER_DAY,
  SLOTS_PER_HOUR,
  SRV_BUSY,
  SRV_DOWN,
  SRV_IDLE,
  STN_BELT,
  STN_CHECKIN,
  STN_CUSTOMS,
  STN_GATE,
  STN_IMIGRASI,
  STN_SECURITY,
} from "./config.ts";
import type { Config, RunParams } from "./config.ts";
import { emptyAgent, poolAlloc, poolFree } from "./passenger.ts";
import type { Agent } from "./passenger.ts";
import { BP, FNV_OFFSET, MASK64, Xoshiro256, fnv1aI64, fnv1aU32Pair, splitmix64 } from "./rng.ts";
import { Station } from "./station.ts";

export const BR_ARR_ACC = 0;
export const BR_ARR_REJ = 1;
export const BR_SVC_DONE = 2;
export const BR_SEC_EXTRA = 3;
export const BR_SRV_BREAK = 4;
export const BR_SRV_REPAIR = 5;
export const BR_ABANDON = 6;
export const BR_MISSED_FLT = 7;
export const BR_REBOOK = 8;
export const BR_BAG_MISS = 9;
export const BR_BELT_BLOCK = 10;
export const BR_ROSTER_DROP = 11;
export const BR_ROUTE_INTL = 12;
export const BR_ROUTE_TRANSIT = 13;
export const BR_BIZ_PRIO = 14;

export class Sim {
  cfg: Config;
  rp: RunParams;
  rng: Xoshiro256[];
  pool: Agent[];
  freeStk: number[];
  stn: Station[];
  arrProf: number[];
  roster: number[][];
  branches = new Array<number>(15).fill(0);
  totalArr = 0;
  totalDone = 0;
  totalLost = 0;
  totalMissed = 0;
  totalRebook = 0;
  totalAban = 0;
  totalBags = 0;
  totalBagsMiss = 0;
  totalBeltBlock = 0;
  totalBreakdowns = 0;
  sysSumL = 0;
  sysSumSojourn = 0;
  invViolations = 0;
  eventHash = FNV_OFFSET;
  clsWaitSum = [0, 0, 0];
  clsWaitCnt = [0, 0, 0];
  clsSojournSum = [0, 0, 0];
  clsSojournCnt = [0, 0, 0];
  stnWaitSum = new Array<number>(N_STN).fill(0);
  stnWaitCnt = new Array<number>(N_STN).fill(0);
  hrArr = new Array<number>(24).fill(0);
  hrUtilSum: number[][] = Array.from({ length: 24 }, () => new Array<number>(N_STN).fill(0));
  waitSamples: number[] = [];
  sojournSamples: number[] = [];
  auditSlots: number;
  qSeries: Int32Array;
  srvSeries: Int32Array;
  stnEv: Array<Array<[number, number, number]>>;

  constructor(c: Config, rp: RunParams) {
    this.cfg = c;
    this.rp = rp;
    const ms = BigInt(c.kv["master_seed"]) & MASK64;
    this.rng = Array.from({ length: 8 }, (_, i) => {
      const [, sd] = splitmix64(ms ^ fnv1aU32Pair(200, i));
      return new Xoshiro256(sd);
    });
    this.pool = Array.from({ length: MAX_POOL }, () => emptyAgent(0));
    this.freeStk = Array.from({ length: MAX_POOL }, (_, i) => MAX_POOL - 1 - i);
    const caps = [1000, 1000, 1000, 1000, rp.beltCap, 1000];
    const prios = [1, 1, 0, 0, 0, 0];
    this.stn = Array.from({ length: N_STN }, (_, i) => {
      let cm = rp.srvOvr[i] > 0 ? rp.srvOvr[i] : c.srvCnt[i];
      if (cm > MAX_SRV) cm = MAX_SRV;
      return new Station(cm, caps[i], prios[i]);
    });
    this.arrProf = [180, 140, 120, 150, 260, 420, 680, 860, 920, 820, 700, 660, 680, 740, 820, 880, 940, 900, 760, 620, 480, 360, 280, 210];
    this.roster = Array.from({ length: 24 }, (_, h) =>
      Array.from({ length: N_STN }, (_, st) => {
        const cm = this.stn[st].cMax;
        return h >= 6 && h <= 20 ? cm : cm > 1 ? cm - 1 : 1;
      })
    );
    this.auditSlots = rp.auditOn !== 0 ? Math.min(rp.days, AUDIT_DAYS) * SLOTS_PER_DAY : 0;
    this.qSeries = new Int32Array(this.auditSlots * N_STN);
    this.srvSeries = new Int32Array(this.auditSlots * 3);
    this.stnEv = Array.from({ length: N_STN }, () => []);
  }

  effBp(base: number): number {
    return this.rp.timescaleDiv2 !== 0 ? Math.trunc(base / 2) : base;
  }

  enqueueAgent(stnId: number, aidx: number, slot: number, toHead: boolean): void {
    const st = this.stn[stnId];
    if (st.len + st.busyCnt >= st.cap) {
      this.branches[BR_ARR_REJ]++;
      st.blockCnt++;
      if (stnId === STN_BELT) this.totalBeltBlock++;
      if (this.pool[aidx].isBag === 0) {
        this.totalLost++;
        this.sysSumSojourn += slot - this.pool[aidx].arrSlot;
      }
      poolFree(this, aidx);
      return;
    }
    this.pool[aidx].curStn = stnId;
    this.pool[aidx].qEnterSlot = slot;
    st.arrCnt++;
    if (toHead) st.pushHead(aidx);
    else st.pushTail(aidx);
    if (this.rp.auditOn !== 0 && slot < this.auditSlots) {
      this.stnEv[stnId].push([slot, stnId, 1]);
    }
  }
}

export function simStepSlot(s: Sim, slot: number): void {
  const hr = Math.trunc((slot % SLOTS_PER_DAY) / SLOTS_PER_HOUR);
  const comp: Array<[number, number]> = [];

  let arrP = Math.trunc((s.arrProf[hr] * s.rp.arrScaleBp) / BP);
  if (s.rp.timescaleDiv2 !== 0) arrP = Math.trunc(arrP / 2);
  if (s.rng[0].bernoulli(arrP) !== 0) {
    const aidx = poolAlloc(s);
    if (aidx >= 0) {
      const ag = s.pool[aidx];
      const uCls = s.rng[0].randInt(1, BP);
      ag.cls = uCls <= s.cfg.kv["biz_prob_bp"] ? CLS_BIZ : uCls <= s.cfg.kv["biz_prob_bp"] + s.cfg.kv["transit_prob_bp"] ? CLS_TRANSIT : CLS_ECON;
      ag.dest = s.rng[0].bernoulli(s.cfg.kv["intl_prob_bp"]) !== 0 ? DEST_INTL : DEST_DOM;
      ag.nBags = s.rng[0].randInt(0, 2);
      ag.arrSlot = slot;
      ag.flightHr = hr;
      ag.depSlot = slot + s.rng[0].randInt(110, 220);
      s.totalArr++;
      s.hrArr[hr]++;
      s.branches[BR_ARR_ACC]++;
      s.eventHash = fnv1aI64(s.eventHash, slot * 100 + ag.cls * 10 + ag.dest);
      let startStn = STN_CHECKIN;
      if (ag.cls === CLS_TRANSIT) {
        s.branches[BR_ROUTE_TRANSIT]++;
        startStn = STN_SECURITY;
      }
      if (ag.nBags > 0 && ag.cls !== CLS_TRANSIT) {
        for (let b = 0; b < ag.nBags; b++) {
          const bidx = poolAlloc(s);
          if (bidx >= 0) {
            const bg = s.pool[bidx];
            bg.isBag = 1;
            bg.cls = ag.cls;
            bg.dest = ag.dest;
            bg.arrSlot = slot;
            bg.depSlot = ag.depSlot;
            s.totalBags++;
            s.enqueueAgent(STN_BELT, bidx, slot, false);
          }
        }
      }
      s.enqueueAgent(startStn, aidx, slot, false);
    }
  }

  for (let st = 0; st < N_STN; st++) {
    const sn = s.stn[st];
    const baseMu = s.rp.muOvr[st] > 0 ? s.rp.muOvr[st] : s.cfg.muBp[st];
    for (let idxK = 0; idxK < sn.cMax; idxK++) {
      const k = s.rp.permuteSrv !== 0 ? (idxK + 1) % sn.cMax : idxK;
      const sv = sn.srv[k];
      if (sv.status !== SRV_BUSY || sv.agentIdx < 0) continue;
      const mu = sv.extraSec !== 0 ? s.cfg.kv["security_extra_mu_bp"] : baseMu;
      if (s.rng[1].bernoulli(s.effBp(mu)) !== 0) {
        comp.push([sv.agentIdx, st]);
        sv.agentIdx = -1;
        sv.status = SRV_IDLE;
        sv.extraSec = 0;
        sn.busyCnt--;
        sn.depCnt++;
        s.branches[BR_SVC_DONE]++;
        if (s.rp.auditOn !== 0 && slot < s.auditSlots) s.stnEv[st].push([slot, st, -1]);
      }
    }
  }

  for (const [aidx, st] of comp) {
    const ag = s.pool[aidx];
    if (ag.isBag !== 0) {
      if (s.rng[4].bernoulli(s.cfg.kv["bag_miss_bp"]) !== 0) {
        s.totalBagsMiss++;
        s.branches[BR_BAG_MISS]++;
      }
      poolFree(s, aidx);
      continue;
    }
    let nextSt = -1;
    if (st === STN_CHECKIN) nextSt = STN_SECURITY;
    else if (st === STN_SECURITY) {
      if (ag.dest === DEST_INTL) {
        nextSt = STN_IMIGRASI;
        s.branches[BR_ROUTE_INTL]++;
      } else nextSt = STN_GATE;
    } else if (st === STN_IMIGRASI) nextSt = STN_GATE;
    else if (st === STN_GATE) {
      const lead = s.rp.timescaleDiv2 !== 0 ? s.rp.gateCloseLead * 2 : s.rp.gateCloseLead;
      if (slot >= ag.depSlot - lead) {
        s.totalMissed++;
        s.branches[BR_MISSED_FLT]++;
        if (s.rng[2].bernoulli(s.cfg.kv["rebook_prob_bp"]) !== 0) {
          s.totalRebook++;
          s.branches[BR_REBOOK]++;
        }
        s.totalLost++;
        s.sysSumSojourn += slot - ag.arrSlot;
        poolFree(s, aidx);
        continue;
      }
      if (s.rp.customsOn !== 0 && ag.dest === DEST_INTL) nextSt = STN_CUSTOMS;
    }
    if (nextSt >= 0) {
      s.enqueueAgent(nextSt, aidx, slot, false);
    } else {
      s.totalDone++;
      const soj = slot - ag.arrSlot;
      s.sysSumSojourn += soj;
      s.clsSojournSum[ag.cls] += soj;
      s.clsSojournCnt[ag.cls]++;
      if (s.sojournSamples.length < 50000) s.sojournSamples.push(soj);
      poolFree(s, aidx);
    }
  }

  for (let st = 0; st < N_STN; st++) {
    const sn = s.stn[st];
    for (let idxK = 0; idxK < sn.cMax; idxK++) {
      const k = s.rp.permuteSrv !== 0 ? (idxK + 1) % sn.cMax : idxK;
      const sv = sn.srv[k];
      if (sv.enabled === 0 || sv.status !== SRV_IDLE || sn.len <= 0) continue;
      const [aidx, nb] = sn.popNext(s.pool, s.rp.priorityOn, s.branches[BR_BIZ_PRIO]);
      s.branches[BR_BIZ_PRIO] = nb;
      if (aidx >= 0) {
        const w = slot - s.pool[aidx].qEnterSlot;
        sn.sumW += w;
        sn.waitCnt++;
        s.stnWaitSum[st] += w;
        s.stnWaitCnt[st]++;
        if (s.pool[aidx].isBag === 0) {
          const cl = s.pool[aidx].cls;
          s.clsWaitSum[cl] += w;
          s.clsWaitCnt[cl]++;
          if (s.waitSamples.length < 50000) s.waitSamples.push(w);
        }
        const extra = st === STN_SECURITY && s.rng[5].bernoulli(s.cfg.kv["security_extra_prob_bp"]) !== 0 ? 1 : 0;
        if (extra !== 0) s.branches[BR_SEC_EXTRA]++;
        sv.status = SRV_BUSY;
        sv.agentIdx = aidx;
        sv.extraSec = extra;
        sn.busyCnt++;
      }
    }
  }

  const pat = s.rp.timescaleDiv2 !== 0 ? s.cfg.kv["patience_slots"] * 2 : s.cfg.kv["patience_slots"];
  for (let st = 0; st < N_STN; st++) {
    if (st === STN_BELT) continue;
    const sn = s.stn[st];
    let pos = 0;
    while (pos < sn.len) {
      const aidx = sn.q[(sn.head + pos) % MAX_Q];
      if (slot - s.pool[aidx].qEnterSlot >= pat) {
        sn.removeAt(pos);
        sn.abanCnt++;
        sn.depCnt++;
        s.totalAban++;
        s.totalLost++;
        s.branches[BR_ABANDON]++;
        s.sysSumSojourn += slot - s.pool[aidx].arrSlot;
        if (s.rp.auditOn !== 0 && slot < s.auditSlots) s.stnEv[st].push([slot, st, -1]);
        poolFree(s, aidx);
      } else {
        pos++;
      }
    }
  }

  if (s.rp.breakdownOn !== 0 && s.rp.breakProbBp > 0) {
    const bp = s.effBp(s.rp.breakProbBp);
    const rpb = s.effBp(s.rp.repairProbBp);
    for (let st = 0; st < N_STN; st++) {
      const sn = s.stn[st];
      for (let idxK = 0; idxK < sn.cMax; idxK++) {
        const k = s.rp.permuteSrv !== 0 ? (idxK + 1) % sn.cMax : idxK;
        const sv = sn.srv[k];
        if (sv.status === SRV_DOWN) {
          if (s.rng[3].bernoulli(rpb) !== 0) {
            sv.status = SRV_IDLE;
            s.branches[BR_SRV_REPAIR]++;
          }
        } else if (sv.enabled !== 0 && s.rng[3].bernoulli(bp) !== 0) {
          if (sv.status === SRV_BUSY && sv.agentIdx >= 0) {
            sn.pushHead(sv.agentIdx);
            sv.agentIdx = -1;
            sn.busyCnt--;
          }
          sv.status = SRV_DOWN;
          s.totalBreakdowns++;
          s.branches[BR_SRV_BREAK]++;
        }
      }
    }
  }

  if (s.rp.rosterOn !== 0 && slot % SLOTS_PER_HOUR === 0) {
    for (let st = 0; st < N_STN; st++) {
      const sn = s.stn[st];
      const target = s.roster[hr][st];
      if (target < sn.cActive) s.branches[BR_ROSTER_DROP]++;
      sn.cActive = target;
      for (let idxK = 0; idxK < sn.cMax; idxK++) {
        const k = s.rp.permuteSrv !== 0 ? (idxK + 1) % sn.cMax : idxK;
        sn.srv[k].enabled = idxK < target ? 1 : 0;
      }
    }
  }

  let ti = 0;
  let tb = 0;
  let td = 0;
  for (let st = 0; st < N_STN; st++) {
    const sn = s.stn[st];
    const inSt = sn.len + sn.busyCnt;
    sn.sumL += inSt;
    sn.busySlots += sn.busyCnt;
    if (st !== STN_BELT) s.sysSumL += inSt;
    s.hrUtilSum[hr][st] += sn.busyCnt;
    if (sn.busyCnt > sn.cMax || inSt > sn.cap) s.invViolations++;
    if (s.rp.auditOn !== 0 && slot < s.auditSlots) s.qSeries[slot * N_STN + st] = inSt;
    for (let k = 0; k < sn.cMax; k++) {
      const sv = sn.srv[k];
      if (sv.status === SRV_IDLE && sv.enabled !== 0) ti++;
      else if (sv.status === SRV_BUSY) tb++;
      else if (sv.status === SRV_DOWN) td++;
    }
  }
  if (s.rp.auditOn !== 0 && slot < s.auditSlots) {
    s.srvSeries[slot * 3] = ti;
    s.srvSeries[slot * 3 + 1] = tb;
    s.srvSeries[slot * 3 + 2] = td;
  }
}
