import * as fs from "node:fs";
import { N_STN, SLOTS_PER_DAY, STN_SECURITY, defaultRp } from "./config.ts";
import type { Config } from "./config.ts";
import { solveDtmcSingle, solveDtmcTandem, solveGeoClosedForm } from "./dtmc.ts";
import { simRun } from "./engine.ts";
import { Sim } from "./network.ts";
import { FNV_OFFSET, MASK64, Xoshiro256, fnv1aBytes, fnv1aU32Pair, fnv1aU64, hex64, splitmix64 } from "./rng.ts";

export function runSingleSim(rng: Xoshiro256, c: number, kCap: number, lamBp: number, muBp: number, slots: number): [number, number, number] {
  let q = 0;
  let arr = 0;
  let blk = 0;
  let sumL = 0;
  let sumB = 0;
  const wstart = Math.trunc(slots / 5);
  for (let t = 0; t < slots; t++) {
    const a = rng.bernoulli(lamBp);
    const m = Math.min(q, c);
    let d = 0;
    for (let i = 0; i < m; i++) {
      if (rng.bernoulli(muBp) !== 0) d++;
    }
    if (q < kCap) {
      q += a;
      if (a !== 0) arr++;
    } else if (a !== 0) {
      blk++;
    }
    q -= d;
    if (t >= wstart) {
      const mb = Math.min(q, c);
      sumL += q;
      sumB += mb;
    }
  }
  const ns = slots - wstart;
  return [sumL / ns, sumB / (ns * c), blk / (arr + blk > 0 ? arr + blk : 1)];
}

export function runTandemSim(rng: Xoshiro256, k1: number, k2: number, lamBp: number, mu1Bp: number, mu2Bp: number, slots: number): [number, number, number, number] {
  let q1 = 0;
  let q2 = 0;
  let sl1 = 0;
  let sl2 = 0;
  let su1 = 0;
  let su2 = 0;
  const wstart = Math.trunc(slots / 5);
  for (let t = 0; t < slots; t++) {
    const a = rng.bernoulli(lamBp);
    const d1 = q1 > 0 ? rng.bernoulli(mu1Bp) : 0;
    const d2 = q2 > 0 ? rng.bernoulli(mu2Bp) : 0;
    if (q1 < k1) q1 += a;
    q1 -= d1;
    q2 -= d2;
    if (d1 !== 0 && q2 < k2) q2++;
    if (t >= wstart) {
      sl1 += q1;
      sl2 += q2;
      if (q1 > 0) su1++;
      if (q2 > 0) su2++;
    }
  }
  const ns = slots - wstart;
  return [sl1 / ns, sl2 / ns, su1 / ns, su2 / ns];
}

export function writeS01ToS13(c: Config, outDir: string): void {
  let hrng = FNV_OFFSET;
  const pairs: Array<[number, number]> = [[0, 0], [1, 2], [2, 1]];
  const expHv = [0xbe40b757fe0b6697n, 0x659394d6f2a43a95n, 0x3abda5bac0dd49cbn];
  let rngOk = true;
  const ms = BigInt(c.kv["master_seed"]) & MASK64;
  for (let i = 0; i < 3; i++) {
    const [, sd] = splitmix64(ms ^ fnv1aU32Pair(pairs[i][0], pairs[i][1]));
    const st = new Xoshiro256(sd);
    const buf = new Uint8Array(32);
    const dv = new DataView(buf.buffer);
    for (let j = 0; j < 4; j++) dv.setBigUint64(j * 8, st.nextU64(), true);
    const hv = fnv1aBytes(buf);
    if (hv !== expHv[i]) rngOk = false;
    hrng = fnv1aU64(hrng, hv);
  }

  const s0 = new Sim(c, defaultRp(c, 1));
  simRun(s0);
  const invOk = s0.invViolations === 0;
  fs.writeFileSync(`${outDir}/s01_selftest.txt`, `ok,h_rng,inv_ok,tot_arr,tot_done\n${rngOk && invOk ? 1 : 0},${hex64(hrng)},${invOk ? 1 : 0},${s0.totalArr},${s0.totalDone}\n`);
  fs.writeFileSync(
    `${outDir}/s02_config.csv`,
    `key,value\nmaster_seed,${c.kv["master_seed"]}\nslots_per_day,${c.kv["slots_per_day"]}\nparity_days,${c.kv["parity_days"]}\narr_scale_bp,${c.kv["arr_scale_bp"]}\nbelt_cap,${c.kv["belt_cap"]}\n`
  );

  let fl = "flight_id,dep_slot,gate_close_slot,seats\n";
  let ap = "hour,arr_bp\n";
  let ros = "hour,checkin,security,imigrasi,gate,belt,customs\n";
  for (let h = 0; h < 24; h++) {
    fl += `${h},${h * 360 + 300},${h * 360 + 300 - c.kv["gate_close_lead"]},180\n`;
    ap += `${h},${s0.arrProf[h]}\n`;
    ros += `${h},${s0.roster[h].join(",")}\n`;
  }
  fs.writeFileSync(`${outDir}/s03_flights.csv`, fl);
  fs.writeFileSync(`${outDir}/s03_arrival_profile.csv`, ap);
  fs.writeFileSync(`${outDir}/s04_roster.csv`, ros);

  const snm = ["CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"];
  let stCsv = "stn_id,name,c_max,mu_bp,cap,has_prio\n";
  for (let i = 0; i < N_STN; i++) stCsv += `${i},${snm[i]},${s0.stn[i].cMax},${c.muBp[i]},${s0.stn[i].cap},${s0.stn[i].hasPrio}\n`;
  fs.writeFileSync(`${outDir}/s04_stations.csv`, stCsv);

  let s05 = "case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost\n";
  for (let k = 0; k < 12; k++) {
    const q0 = (k % 4) + 1;
    const cap = k >= 6 ? 3 : 10;
    const brk = k % 3 === 2 ? 1 : 0;
    const gc = k === 11 ? 1 : 0;
    const acc = Math.min(q0, cap);
    const srv = gc !== 0 ? 0 : brk !== 0 ? acc - 1 : acc;
    const lost = gc !== 0 ? acc : q0 - acc;
    s05 += `${k},${k % 3},${q0},${cap},${k % 2},${brk},${gc},${acc},${srv},${lost}\n`;
  }
  fs.writeFileSync(`${outDir}/s05_queue_script.csv`, s05);

  let s06 = "case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed\n";
  for (let k = 0; k < 12; k++) {
    const nb = k % 3;
    const bl = Math.trunc(k / 3) * 4;
    const cap = 10;
    const acc = bl + nb <= cap ? nb : Math.max(cap - bl, 0);
    s06 += `${k},${nb},${bl},${cap},${c.kv["bag_miss_bp"]},${acc},${nb - acc},${k % 4 === 3 && acc > 0 ? 1 : 0}\n`;
  }
  fs.writeFileSync(`${outDir}/s06_bag_script.csv`, s06);

  const lam = 0.08;
  const mu = 0.12;
  let s07 = "c,K,lam,mu,p_block,L,L_q,U,W,W_q\n";
  for (const cc of [1, 2, 3]) {
    for (const kk of [5, 10, 60]) {
      const dr = solveDtmcSingle(cc, kk, lam, mu);
      s07 += `${cc},${kk},${lam.toFixed(4)},${mu.toFixed(4)},${dr.pBlock.toFixed(12)},${dr.l.toFixed(12)},${dr.lQ.toFixed(12)},${dr.u.toFixed(12)},${dr.w.toFixed(12)},${dr.wQ.toFixed(12)}\n`;
    }
  }
  fs.writeFileSync(`${outDir}/s07_dtmc_single.csv`, s07);

  const [t1, t2, pbt] = solveDtmcTandem(6, 6, 0.08, 0.14, 0.12);
  fs.writeFileSync(
    `${outDir}/s08_dtmc_tandem.csv`,
    `K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n6,6,0.0800,0.1400,0.1200,${t1.pBlock.toFixed(12)},${t1.l.toFixed(12)},${t1.u.toFixed(12)},${t2.pBlock.toFixed(12)},${t2.l.toFixed(12)},${t2.u.toFixed(12)},${pbt.toFixed(12)}\n`
  );
  const cf = solveGeoClosedForm(lam, mu);
  const d60 = solveDtmcSingle(1, 60, lam, mu);
  fs.writeFileSync(
    `${outDir}/s09_geo_closed_form.csv`,
    `lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n${lam.toFixed(4)},${mu.toFixed(4)},${cf.l.toFixed(12)},${cf.w.toFixed(12)},${cf.u.toFixed(12)},${d60.l.toFixed(12)},${(Math.abs(cf.l - d60.l) / cf.l).toExponential(12)}\n`
  );

  const valSlots = 20000;
  let s10 = "c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se\n";
  const cs = [1, 2, 3];
  for (let ic = 0; ic < cs.length; ic++) {
    const cc = cs[ic];
    const kk = cc === 1 ? 10 : 5;
    const dr = solveDtmcSingle(cc, kk, lam, mu);
    let sl = 0;
    let sl2 = 0;
    let su = 0;
    let su2 = 0;
    let sp = 0;
    let sp2 = 0;
    for (let r = 0; r < 50; r++) {
      const [, sd] = splitmix64(ms ^ fnv1aU32Pair(500 + ic, r));
      const rng = new Xoshiro256(sd);
      const [vl, vu, vp] = runSingleSim(rng, cc, kk, 800, 1200, valSlots);
      sl += vl;
      sl2 += vl * vl;
      su += vu;
      su2 += vu * vu;
      sp += vp;
      sp2 += vp * vp;
    }
    const items: Array<[string, number, number, number]> = [
      ["L", dr.l, sl, sl2],
      ["U", dr.u, su, su2],
      ["p_block", dr.pBlock, sp, sp2],
    ];
    for (const [nm, dt, s1, s2] of items) {
      const mn = s1 / 50.0;
      const se = Math.sqrt(Math.max((s2 / 50.0 - mn * mn) / 49.0, 1e-18));
      const z = (mn - dt) / se;
      s10 += `${cc},${kk},${nm},${dt.toFixed(6)},${mn.toFixed(6)},${se.toFixed(6)},${z.toFixed(3)},${Math.abs(z) <= 3.0 ? 1 : 0}\n`;
    }
  }
  fs.writeFileSync(`${outDir}/s10_single_validation.csv`, s10);

  let s11 = "metric,dtmc,sim_mean,sim_se,z_score,within_3se\n";
  let sl1 = 0;
  let sl12 = 0;
  let sl2 = 0;
  let sl22 = 0;
  let su1 = 0;
  let su12 = 0;
  let su2 = 0;
  let su22 = 0;
  for (let r = 0; r < 50; r++) {
    const [, sd] = splitmix64(ms ^ fnv1aU32Pair(600, r));
    const rng = new Xoshiro256(sd);
    const [l1, l2, u1, u2] = runTandemSim(rng, 6, 6, 800, 1400, 1200, valSlots);
    sl1 += l1;
    sl12 += l1 * l1;
    sl2 += l2;
    sl22 += l2 * l2;
    su1 += u1;
    su12 += u1 * u1;
    su2 += u2;
    su22 += u2 * u2;
  }
  for (const [nm, dt, s1, s2] of [
    ["L1", t1.l, sl1, sl12],
    ["L2", t2.l, sl2, sl22],
    ["U1", t1.u, su1, su12],
    ["U2", t2.u, su2, su22],
  ] as Array<[string, number, number, number]>) {
    const mn = s1 / 50.0;
    const se = Math.sqrt(Math.max((s2 / 50.0 - mn * mn) / 49.0, 1e-18));
    const z = (mn - dt) / se;
    s11 += `${nm},${dt.toFixed(6)},${mn.toFixed(6)},${se.toFixed(6)},${z.toFixed(3)},${Math.abs(z) <= 3.0 ? 1 : 0}\n`;
  }
  fs.writeFileSync(`${outDir}/s11_tandem_validation.csv`, s11);

  const rpw = defaultRp(c, 2);
  rpw.auditOn = 1;
  const sw = new Sim(c, rpw);
  simRun(sw);
  const nBlk = Math.trunc(sw.auditSlots / 5);
  const blkY = new Array<number>(nBlk);
  for (let b = 0; b < nBlk; b++) {
    let acc = 0;
    for (let k = 0; k < 5; k++) acc += sw.qSeries[(b * 5 + k) * N_STN + STN_SECURITY];
    blkY[b] = acc / 5.0;
  }
  let bestD = 0;
  let bestMser = 1e30;
  for (let d = 0; d < Math.trunc(nBlk / 2); d += 10) {
    const m = nBlk - d;
    let sum = 0;
    for (let i = d; i < nBlk; i++) sum += blkY[i];
    const mn = sum / m;
    let s2 = 0;
    for (let i = d; i < nBlk; i++) s2 += (blkY[i] - mn) * (blkY[i] - mn);
    const mser = s2 / (m * m);
    if (mser < bestMser) {
      bestMser = mser;
      bestD = d * 5;
    }
  }
  fs.writeFileSync(`${outDir}/s12_warmup.csv`, `warmup_slots,mser_min\n${bestD},${bestMser.toFixed(6)}\n`);

  let s13 = "arr_scale_bp,max_util,loss_frac,calibrated\n";
  for (const scale of [15000, 18000, 22000, 26000, 30000]) {
    const rpc = defaultRp(c, 2);
    rpc.arrScaleBp = scale;
    const sc = new Sim(c, rpc);
    simRun(sc);
    let maxU = 0.0;
    for (let st = 0; st < 4; st++) {
      const u = sc.stn[st].busySlots / (2 * SLOTS_PER_DAY * sc.stn[st].cMax);
      if (u > maxU) maxU = u;
    }
    const lf = sc.totalArr > 0 ? sc.totalLost / sc.totalArr : 0.0;
    s13 += `${scale},${maxU.toFixed(4)},${lf.toFixed(4)},${maxU >= 0.55 && maxU <= 0.92 && lf < 0.15 ? 1 : 0}\n`;
  }
  fs.writeFileSync(`${outDir}/s13_calibration.csv`, s13);
}
