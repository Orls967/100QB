import { runAllocAuction } from "./allocation.ts";
import { simEodClearing, simRiskCheck } from "./clearing.ts";
import type { Sim } from "./entities.ts";
import { ACC_CUST, ACC_FEE, ACC_SUPP, CARRIER_BASE_LT, carrierRate, checkInvariants, NUM_ACCOUNTS, NUM_CARRIERS, NUM_STORES, SKU_PROD, SKU_RETAIL, SKU_WHOLESALE, STORE_WH, transferCash } from "./inventory.ts";
import { fnv1aI64, type SimEvent } from "./rng.ts";

function doShip(s: Sim, src: number, dst: number, sku: number, qty: number, day: number, tick: number, tAbs: number, bin: number): void {
  if (qty <= 0 || s.acc[src].onHand[sku] <= 0) return;
  if (qty > s.acc[src].onHand[sku]) qty = s.acc[src].onHand[sku];
  const whCap = 20000 * s.capScale, cCap = 12000 * s.capScale;
  if (dst >= 10 && dst < 14 && s.acc[dst].onHand.reduce((a, b) => a + b, 0) + qty > whCap) { s.branches[4]++; return; }
  let cid = (src + dst + sku) % NUM_CARRIERS;
  if (s.carrierDayVol[cid] + qty > cCap) {
    cid = (cid + 1) % NUM_CARRIERS;
    if (s.carrierDayVol[cid] + qty > cCap) { s.branches[5]++; return; }
  }
  const rate = carrierRate(s.carrierDayVol[cid]);
  s.carrierTierCnt[rate === 200 ? 0 : (rate === 160 ? 1 : 2)]++;
  s.carrierDayVol[cid] += qty;
  s.carrierVolTot += qty;
  if (bin >= 0 && bin < 100) s.binOrders[sku][bin] += qty;
  const goods = qty * SKU_WHOLESALE[sku], freight = qty * rate;
  s.acc[src].onHand[sku] -= qty; s.acc[dst].inTransit[sku] += qty;
  transferCash(s, dst, src, goods, -1, day, tick);
  transferCash(s, dst, 17 + cid, freight, -1, day, tick);
  if (day < 6) s.dailyCostEchelon[day][2] += freight;
  const delay = s.rng[17 + cid].randInt(0, 25);
  if (delay > 18) s.branches[15]++;
  const lt = CARRIER_BASE_LT[cid] + delay;
  if (s.leadtimePairs.length < 256) s.leadtimePairs.push([lt, qty]);
  if (s.logEnabled) s.ships.push([day, tick, src, dst, cid, sku, qty, goods, freight, 1]);
  s.heap.push({ time: tAbs + lt, seq: s.nextSeq++, etype: 2, src, dst, sku, qty, aux: cid });
}

function hashEv(s: Sim, ev: SimEvent): void {
  for (const v of [ev.time, ev.etype, ev.src, ev.dst, ev.sku, ev.qty]) s.eventHash = fnv1aI64(s.eventHash, v);
}

export function simRunDays(s: Sim, days: number, evPerDay: number, safetyOvr: number, disruptOvr: number, transshipOn: number): void {
  const tDay = s.cfg.ticks_per_day, cs = s.capScale; let evProc = 0;
  const chkStep = Math.max(Math.floor(evPerDay / 8), 1), shkStep = Math.max(Math.floor(evPerDay / 250), 1);
  for (let d = 0; d < days; d++) {
    const dayBase = d * tDay; s.carrierDayVol = [0, 0, 0];
    s.carrierCapTot += NUM_CARRIERS * 12000 * cs;
    for (let f = 14; f < 17; f++) {
      if (s.acc[f].plantState === 1) { s.acc[f].disruptRem--; if (s.acc[f].disruptRem <= 0) s.acc[f].plantState = 2; }
      else if (s.acc[f].plantState === 2) s.acc[f].plantState = 0;
    }
    if (d > 0) {
      for (let wh = 10; wh < 13; wh++) {
        for (let sku = 0; sku < 4; sku++) {
          const sts: number[] = [], reqs: number[] = [], seqs: number[] = [], prio: number[] = [];
          const sLvl = (safetyOvr > 0 ? safetyOvr : 450) * cs;
          for (let st = 0; st < NUM_STORES; st++) {
            if (STORE_WH[st] !== wh) continue;
            const inv = s.acc[st].onHand[sku] + s.acc[st].inTransit[sku];
            if (inv < sLvl) { sts.push(st); reqs.push(sLvl - inv); seqs.push(s.nextSeq++); prio.push(st % 3 === 0 ? 2 : 1); }
          }
          if (sts.length > 0) {
            const fills = runAllocAuction(s, wh, sku, sts, reqs, seqs, prio);
            const b0 = Math.floor((d * 100) / days);
            for (let i = 0; i < sts.length; i++) {
              s.logLine("allocations.csv", `${d},${wh},${sku},${sts[i]},${reqs[i]},${fills[i]}`);
              if (fills[i] > 0) doShip(s, wh, sts[i], sku, fills[i], d, 0, dayBase, b0);
            }
          }
        }
      }
    }
    for (let e = 0; e < evPerDay; e++) {
      const tAbs = dayBase + Math.floor((e * tDay) / evPerDay);
      const st = s.rng[NUM_ACCOUNTS].randInt(0, 9), sku = s.rng[NUM_ACCOUNTS].randInt(0, 3);
      const qty = s.rng[NUM_ACCOUNTS].randInt(0, 9) === 0 ? s.rng[NUM_ACCOUNTS].randInt(12, 24) : s.rng[NUM_ACCOUNTS].randInt(1, 5);
      s.heap.push({ time: tAbs, seq: s.nextSeq++, etype: 0, src: ACC_CUST, dst: st, sku, qty, aux: 0 });
      if (e % 20 === 0) { const f = Math.floor(e / 20) % 3; s.heap.push({ time: tAbs + 1, seq: s.nextSeq++, etype: 3, src: 14 + f, dst: 10 + f, sku: f, qty: 130, aux: 0 }); }
      if (e % 40 === 0) { const kq = s.rng[13].randInt(15, 35); s.heap.push({ time: tAbs + 2, seq: s.nextSeq++, etype: 4, src: 13, dst: 13, sku: 3, qty: kq, aux: 0 }); }
      if (e % 60 === 0 && transshipOn !== 0) { const tq = s.rng[13].randInt(30, 50); s.heap.push({ time: tAbs + 3, seq: s.nextSeq++, etype: 5, src: 10, dst: 13, sku: 0, qty: tq, aux: 0 }); }
      if (e % chkStep === 0) {
        const dp = disruptOvr >= 0 ? disruptOvr : s.cfg.disrupt_prob_pct;
        if (s.rng[NUM_ACCOUNTS + 1].randInt(1, 100) <= dp) { const f = s.rng[NUM_ACCOUNTS + 1].randInt(0, 2); s.heap.push({ time: tAbs + 4, seq: s.nextSeq++, etype: 7, src: 14 + f, dst: 14 + f, sku: f, qty: 1, aux: 0 }); }
        s.heap.push({ time: tAbs + 6, seq: s.nextSeq++, etype: 8, src: 0, dst: 0, sku: 0, qty: 0, aux: 0 });
      }
      if (e % shkStep === 0) {
        if (s.rng[NUM_ACCOUNTS + 1].randInt(1, 100) <= s.cfg.shock_prob_pct) { const stS = s.rng[NUM_ACCOUNTS + 1].randInt(0, 9), cls = s.rng[NUM_ACCOUNTS + 1].randInt(0, 2); s.heap.push({ time: tAbs + 5, seq: s.nextSeq++, etype: 6, src: ACC_CUST, dst: stS, sku: 0, qty: (cls + 1) * 10, aux: cls }); }
      }
    }
    const dayEnd = (d + 1) * tDay;
    while (!s.heap.empty() && s.heap.peekTime() < dayEnd) {
      const ev = s.heap.pop(), tick = ev.time - dayBase;
      const bin = Math.max(0, Math.min(99, Math.floor(((d * tDay + tick) * 100) / (days * tDay))));
      hashEv(s, ev); evProc++;
      if (evProc <= 5000) s.logLine("orders_log.csv", `${d},${tick},${ev.seq},${ev.etype},${ev.src},${ev.dst},${ev.sku},${ev.qty}`);
      if (ev.etype === 0) {
        const st = ev.dst, sku = ev.sku, qty = ev.qty;
        s.demEventsTot++; s.demUnitsTot += qty;
        s.binDemand[sku][bin] += qty; s.binStoreDem[st][bin] += qty;
        if (d < 6) { s.dailyDemand[d][sku] += qty; s.dailyFillStore[st][0] += qty; }
        const al = s.cfg.ewma_alpha_milli;
        s.acc[st].ewmaMilli[sku] = Math.floor((al * qty * 1000 + (1000 - al) * s.acc[st].ewmaMilli[sku]) / 1000);
        if (s.acc[st].status === 3) s.branches[2]++;
        else if (s.acc[st].onHand[sku] >= qty) {
          s.acc[st].onHand[sku] -= qty; s.totalSold[sku] += qty; s.demUnitsFilled += qty; const rev = qty * SKU_RETAIL[sku];
          transferCash(s, ACC_CUST, st, rev, -1, d, tick);
          if (s.logEnabled) s.sales.push([d, tick, st, sku, qty, rev]);
          s.branches[0]++; if (d < 6) s.dailyFillStore[st][1] += qty;
        } else {
          const part = s.acc[st].onHand[sku];
          if (part > 0) {
            s.acc[st].onHand[sku] = 0; s.totalSold[sku] += part; s.demUnitsFilled += part; const rev = part * SKU_RETAIL[sku];
            transferCash(s, ACC_CUST, st, rev, -1, d, tick);
            if (s.logEnabled) s.sales.push([d, tick, st, sku, part, rev]);
            if (d < 6) s.dailyFillStore[st][1] += part;
          }
          s.acc[st].backorder[sku] += qty - part; s.demBoEvents++; s.branches[1]++;
        }
        const capS = (safetyOvr > 0 ? safetyOvr : 450) * cs;
        const sLvl = Math.floor((capS * 40) / 100);
        const invPos = s.acc[st].onHand[sku] + s.acc[st].inTransit[sku] - s.acc[st].backorder[sku];
        if (invPos < sLvl && s.acc[st].status === 0) {
          s.branches[3]++;
          const srcWh = (sku === 3 || st === 9) ? 13 : (10 + sku);
          doShip(s, srcWh, st, sku, capS - invPos, d, tick, ev.time, bin);
        }
      } else if (ev.etype === 2) {
        if (s.acc[ev.dst].inTransit[ev.sku] >= ev.qty) { s.acc[ev.dst].inTransit[ev.sku] -= ev.qty; s.acc[ev.dst].onHand[ev.sku] += ev.qty; }
        if (s.logEnabled) s.ships.push([d, tick, ev.src, ev.dst, ev.aux, ev.sku, ev.qty, 0, 0, 2]);
      } else if (ev.etype === 3) {
        const qty = s.acc[ev.src].plantState === 1 ? 0 : (s.acc[ev.src].plantState === 2 ? Math.floor(ev.qty / 2) : ev.qty);
        if (qty > 0) {
          const cost = qty * SKU_PROD[ev.sku]; s.acc[ev.src].onHand[ev.sku] += qty; s.totalProduced[ev.sku] += qty;
          transferCash(s, ev.src, ACC_SUPP, cost, -1, d, tick); if (d < 6) s.dailyCostEchelon[d][0] += cost;
          if (s.logEnabled) s.prods.push([d, tick, ev.src, ev.sku, qty, cost]);
          doShip(s, ev.src, ev.dst, ev.sku, qty, d, tick, ev.time, bin);
        }
      } else if (ev.etype === 4) {
        const m = Math.min(s.acc[13].onHand[0], s.acc[13].onHand[1], s.acc[13].onHand[2]);
        const mx = Math.max(s.acc[13].onHand[0], s.acc[13].onHand[1], s.acc[13].onHand[2]);
        const imb = mx - m; if (s.kitImbSeries.length < 64) s.kitImbSeries.push(imb);
        let q = ev.qty; if (m >= q) s.branches[8]++; else { q = m > 0 ? m : 0; s.branches[9]++; }
        if (q > 0) {
          for (let k = 0; k < 3; k++) s.acc[13].onHand[k] -= q;
          s.acc[13].onHand[3] += q; s.totalKits += q; const fee = q * s.cfg.kit_fee;
          transferCash(s, 13, ACC_FEE, fee, -1, d, tick); if (d < 6) s.dailyCostEchelon[d][1] += fee;
          if (s.logEnabled) s.kits.push([d, tick, 13, q, fee, imb]);
        }
      } else if (ev.etype === 5) {
        let minK = 0; for (let k = 1; k < 3; k++) if (s.acc[13].onHand[k] < s.acc[13].onHand[minK]) minK = k;
        let bestW = 10; for (let w = 11; w < 13; w++) if (s.acc[w].onHand[minK] > s.acc[bestW].onHand[minK]) bestW = w;
        s.branches[10]++; s.nTransships++;
        doShip(s, bestW, 13, minK, ev.qty, d, tick, ev.time, bin);
        s.logLine("transshipments.csv", `${d},${tick},${bestW},13,${minK},${ev.qty}`);
      } else if (ev.etype === 7) {
        s.acc[ev.src].plantState = 1; s.acc[ev.src].disruptRem = 1; s.branches[11]++;
        s.logLine("disruptions.csv", `${d},${tick},${ev.src},1,1`);
      } else if (ev.etype === 6) {
        const cls = ev.aux, mag = (cls + 1) * 15 * cs;
        s.branches[12]++; s.shockCnt[cls]++; s.shockResp[cls] += mag;
        if (s.shockEvs.length < 128) s.shockEvs.push([bin, cls, ev.dst]);
        s.logLine("shocks_log.csv", `${d},${tick},${cls},${ev.dst},${mag}`);
        s.heap.push({ time: ev.time + 1, seq: s.nextSeq++, etype: 0, src: ACC_CUST, dst: ev.dst, sku: 0, qty: Math.min(mag, 5000), aux: 0 });
      } else if (ev.etype === 8) simRiskCheck(s, d, tick);
      if (evProc % s.cfg.k_chk === 0) checkInvariants(s);
    }
    simEodClearing(s, d); checkInvariants(s);
  }
}
