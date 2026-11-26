import { writeFileSync } from "node:fs";
import { runAllocAuction } from "./allocation.ts";
import { simRunDays } from "./engine.ts";
import { type Config, Sim } from "./entities.ts";
import { ACC_CUST, ACC_FEE, ACC_SUPP, CARRIER_BASE_LT, checkInvariants, NUM_ACCOUNTS, STORE_WH } from "./inventory.ts";
import { fnv1aBytes, fnv1aI64, fnv1aU32Pair, splitmix64, Xoshiro256 } from "./rng.ts";

const SKU_L = ["A", "B", "C", "D"];
const wf = (d: string, n: string, txt: string): void => writeFileSync(`${d}/${n}`, txt);
const hex64 = (h: bigint): string => h.toString(16).padStart(16, "0");

export function writeS01ToS09(c: Config, outDir: string): void {
  let hrng = 0xcbf29ce484222325n;
  const pairs: [number, number][] = [[0, 0], [1, 2], [2, 1]];
  const expHv = [0xbe40b757fe0b6697n, 0x659394d6f2a43a95n, 0x3abda5bac0dd49cbn];
  let rngOk = true;
  for (let i = 0; i < 3; i++) {
    const [sd] = splitmix64(BigInt(c.master_seed) ^ fnv1aU32Pair(pairs[i][0], pairs[i][1]));
    const st = new Xoshiro256(sd), buf = new Uint8Array(32), dv = new DataView(buf.buffer);
    for (let j = 0; j < 4; j++) dv.setBigUint64(j * 8, st.nextU64(), true);
    const hv = fnv1aBytes(buf);
    if (hv !== expHv[i]) rngOk = false;
    hrng = fnv1aI64(hrng, hv);
  }
  const s = new Sim(c, false, "");
  const fills = runAllocAuction(s, 10, 0, [0, 1, 2], [1000, 1500, 1000], [1, 2, 3], [2, 1, 1]);
  let halloc = 0xcbf29ce484222325n;
  for (const fl of fills) halloc = fnv1aI64(halloc, fl);
  const invOk = checkInvariants(s), ok = rngOk && invOk === 1 && fills.reduce((a, b) => a + b, 0) === 1600 ? 1 : 0;
  wf(outDir, "s01_selftest.txt", `ok,h_rng,h_alloc,inv_ok\n${ok},${hex64(hrng)},${hex64(halloc)},${invOk}\n`);
  wf(outDir, "s02_config.csv", `key,value\nmaster_seed,${c.master_seed}\nnum_days,${c.num_days}\nparity_days,${c.parity_days}\nparity_events_per_day,${c.parity_events_per_day}\nwh_cap,${c.wh_cap}\nkit_fee,${c.kit_fee}\n`);

  const ent = ["id,type,name", ...Array.from({ length: 10 }, (_, i) => `${i},STORE,S${i}`), ...Array.from({ length: 4 }, (_, i) => `${10 + i},WAREHOUSE,W${i}`), ...Array.from({ length: 3 }, (_, i) => `${14 + i},FACTORY,F${i}`), ...Array.from({ length: 3 }, (_, i) => `${17 + i},CARRIER,C${i}`), "20,CUSTOMER,ExtCustomer", "21,SUPPLIER,RawSupplier", "22,FEEPOOL,FeePool"];
  wf(outDir, "s03_entities.csv", ent.join("\n") + "\n");
  wf(outDir, "s03_accounts_init.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n" + s.acc.map((a, i) => `${i},${a.cash},${a.onHand.join(",")},${a.creditLimit},${a.status}`).join("\n") + "\n");
  wf(outDir, "s03_carriers.csv", "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n" + [0, 1, 2].map((i) => `${i},${CARRIER_BASE_LT[i]},${c.carrier_cap_daily},200,160,120`).join("\n") + "\n");
  const edges = ["src,dst,edge_type", ...[0, 1, 2].map((i) => `${14 + i},${10 + i},FAC_WH`), ...[0, 1, 2].map((i) => `13,${10 + i},KIT_TS`), ...Array.from({ length: 10 }, (_, i) => `${STORE_WH[i]},${i},WH_STORE`)];
  wf(outDir, "s03_network_edges.csv", edges.join("\n") + "\n");

  const pol = ["case_id,inv_pos,s_lvl,S_lvl,order_qty,branch"];
  for (let k = 0; k < 60; k++) {
    const inv = ((k * 37) % 1000) - 100, oq = inv < 250 && k % 11 !== 0 ? 800 - inv : 0;
    pol.push(`${k},${inv},250,800,${oq},${k % 11 === 0 ? "SUSPENDED" : (oq > 0 ? "REORDER" : "NO_OP")}`);
  }
  wf(outDir, "s04_policy_script.csv", pol.join("\n") + "\n");

  const al = ["case_id,avail,f0,f1,f2,f3"];
  for (let k = 0; k < 10; k++) {
    const sa = new Sim(c, false, ""); sa.acc[10].onHand[0] = 200 + k * 150;
    const fl = runAllocAuction(sa, 10, 0, [0, 1, 2, 3], [400, 300 + k * 20, 500, 250], [10 + k, 5, 20 - k, 8], [2, 1, 2, 1]);
    al.push(`${k},${sa.acc[10].onHand[0]},${fl.join(",")}`);
  }
  wf(outDir, "s05_alloc_script.csv", al.join("\n") + "\n");

  const sess = ["day,session,start_tick,end_tick"];
  for (let d = 0; d < c.num_days; d++) sess.push(`${d},PRE_DAY,0,0`, `${d},ALLOCATION_AUCTION,0,1`, `${d},CONTINUOUS_OPS,0,${c.ticks_per_day}`, `${d},EOD_CLEARING,${c.ticks_per_day},${c.ticks_per_day}`);
  wf(outDir, "s06_sessions.csv", sess.join("\n") + "\n");
  const [sd0] = splitmix64(BigInt(c.master_seed) ^ fnv1aU32Pair(2, 23)); const r0 = new Xoshiro256(sd0);
  const sh = ["idx,time,seq,etype,dst,sku,qty"];
  for (let i = 0; i < 200; i++) sh.push(`${i},${Math.floor((i * c.ticks_per_day) / 8000)},${i + 1},0,${r0.randInt(0, 9)},${r0.randInt(0, 3)},${r0.randInt(1, 14)}`);
  wf(outDir, "s06_schedule_head.csv", sh.join("\n") + "\n");
  wf(outDir, "s07_shock_schedule.csv", "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n");

  const sw = new Sim(c, false, ""); simRunDays(sw, 1, 2000, 0, 0, 1);
  wf(outDir, "s08_warmup_summary.csv", `days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,${sw.totalSold.join(",")},${sw.invViolations === 0 ? 1 : 0}\n`);
  for (let k = 0; k < 4; k++) wf(outDir, `s08_warmup_sku_${SKU_L[k]}.csv`, `sku,sold,produced,end_snap\n${SKU_L[k]},${sw.totalSold[k]},${sw.totalProduced[k]},${sw.dayInvSnap[0][k]}\n`);
  wf(outDir, "s09_day1_summary.csv", `sold_tot,kits,alloc_cases,inv_violations\n${sw.totalSold.reduce((a, b) => a + b, 0)},${sw.totalKits},${sw.allocCases},${sw.invViolations}\n`);
  const calFr = sw.demUnitsTot > 0 ? sw.demUnitsFilled / sw.demUnitsTot : 1.0;
  const calBo = sw.demEventsTot > 0 ? sw.demBoEvents / sw.demEventsTot : 0.0;
  const calCu = sw.carrierCapTot > 0 ? sw.carrierVolTot / sw.carrierCapTot : 0.0;
  const calOk = (calFr >= 0.85 && calFr <= 0.99 && calBo < 0.20 && calCu >= 0.30 && calCu <= 0.90) ? 1 : 0;
  wf(outDir, "s09_calibration.csv", `fill_rate,bo_event_frac,carrier_util,calibrated\n${calFr.toFixed(6)},${calBo.toFixed(6)},${calCu.toFixed(6)},${calOk}\n`);
}

export function writeS11ToS38(s: Sim, c: Config, outDir: string): void {
  for (let k = 0; k < 4; k++) {
    wf(outDir, `sales_${SKU_L[k]}.csv`, "day,tick,store_id,sku,qty,rev_cents\n" + s.sales.filter((r) => r[3] === k).map((r) => r.join(",")).join("\n") + "\n");
    wf(outDir, `shipments_${SKU_L[k]}.csv`, "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status\n" + s.ships.filter((r) => r[5] === k).map((r) => r.join(",")).join("\n") + "\n");
    wf(outDir, `inventory_snap_${SKU_L[k]}.csv`, "day,sku,on_hand_tot\n" + Array.from({ length: c.num_days }, (_, d) => `${d},${SKU_L[k]},${s.dayInvSnap[d][k]}`).join("\n") + "\n");
  }
  for (let i = 0; i < 10; i++) wf(outDir, `store_daily_${i}.csv`, `store,cash,oh_tot,req_tot,fill_tot\n${i},${s.acc[i].cash},${s.acc[i].onHand.reduce((a, b) => a + b, 0)},${s.dailyFillStore[i][0]},${s.dailyFillStore[i][1]}\n`);
  for (let i = 0; i < 4; i++) { const a = s.acc[10 + i]; wf(outDir, `warehouse_daily_${i}.csv`, `wh,cash,oh_A,oh_B,oh_C,oh_D\n${i},${a.cash},${a.onHand.join(",")}\n`); }
  for (let i = 0; i < 3; i++) wf(outDir, `factory_daily_${i}.csv`, `fac,cash,produced,state\n${i},${s.acc[14 + i].cash},${s.totalProduced[i]},${s.acc[14 + i].plantState}\n`);
  wf(outDir, "production.csv", "day,tick,factory_id,sku,qty,cost_cents\n" + s.prods.map((r) => r.join(",")).join("\n") + "\n");
  wf(outDir, "kit_assembly.csv", "day,tick,wh_id,qty,fee_cents,imbalance\n" + s.kits.map((r) => r.join(",")).join("\n") + "\n");
  wf(outDir, "fees_ledger.csv", "day,tick,from_id,to_id,amount,kind\n" + s.fees.map((r) => r.join(",")).join("\n") + "\n");
  wf(outDir, "entities_final.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n" + s.acc.map((a, i) => `${i},${a.cash},${a.onHand.join(",")},${a.inTransit.join(",")},${a.status}`).join("\n") + "\n");

  wf(outDir, "s11_invariants.csv", `violations,ok\n${s.invViolations},${s.invViolations === 0 ? 1 : 0}\n`);
  wf(outDir, "s12_clearing.csv", `days,fee_pool_cash\n${c.num_days},${s.acc[ACC_FEE].cash}\n`);
  wf(outDir, "s12_credit.csv", `credit_holds,liquidations\n${s.branches[13]},${s.branches[14]}\n`);
  wf(outDir, "s13_risk.csv", `credit_holds,liquidations\n${s.branches[13]},${s.branches[14]}\n`);
  wf(outDir, "s14_factories.csv", `fac,produced\n0,${s.totalProduced[0]}\n1,${s.totalProduced[1]}\n2,${s.totalProduced[2]}\n`);
  wf(outDir, "s15_carrier_fees.csv", `tier,count\n0,${s.carrierTierCnt[0]}\n1,${s.carrierTierCnt[1]}\n2,${s.carrierTierCnt[2]}\n`);
  wf(outDir, "s16_kits.csv", `total_kits,full,partial\n${s.totalKits},${s.branches[8]},${s.branches[9]}\n`);
  wf(outDir, "s16_kit_imbalance_series.csv", "idx,imbalance\n" + s.kitImbSeries.map((v, i) => `${i},${v}`).join("\n") + "\n");
  wf(outDir, "s17_allocations.csv", `cases,nonzero,total_req,total_fill\n${s.allocCases},${s.allocNonzero},${s.allocTotalReq},${s.allocTotalFill}\n`);
  wf(outDir, "s18_cost_entities.csv", "acc,cash_delta\n" + Array.from({ length: 17 }, (_, i) => `${i},${s.acc[i].cash - (i < 10 ? 2000000000 : 5000000000) * s.capScale}`).join("\n") + "\n");
  wf(outDir, "s19_cost_types.csv", "echelon,cost_cents\n" + [0, 1, 2, 3].map((e) => `${e},${Array.from({ length: c.num_days }, (_, d) => s.dailyCostEchelon[d][e]).reduce((a, b) => a + b, 0)}`).join("\n") + "\n");
  wf(outDir, "s19_cost_daily.csv", "day,prod,kit,freight,hold_pen\n" + Array.from({ length: c.num_days }, (_, d) => `${d},${s.dailyCostEchelon[d].join(",")}`).join("\n") + "\n");
  for (let k = 0; k < 4; k++) wf(outDir, `s20_demand_${SKU_L[k]}.csv`, "day,demand\n" + Array.from({ length: c.num_days }, (_, d) => `${d},${s.dailyDemand[d][k]}`).join("\n") + "\n");
  const bwRows: string[] = [];
  for (let k = 0; k < 4; k++) {
    const md = s.binDemand[k].reduce((a, b) => a + b, 0) / 100, mo = s.binOrders[k].reduce((a, b) => a + b, 0) / 100;
    let m2 = 0, m4 = 0, vo = 0, cov1 = 0;
    for (let b = 0; b < 100; b++) {
      const d = s.binDemand[k][b] - md, o = s.binOrders[k][b] - mo;
      m2 += d * d; m4 += d * d * d * d; vo += o * o;
      if (b > 0) cov1 += (s.binDemand[k][b] - md) * (s.binDemand[k][b - 1] - md);
    }
    const kurt = m2 > 0 ? (m4 / 100) / ((m2 / 100) ** 2) : 0.0, acf1 = m2 > 0 ? cov1 / m2 : 0.0, bwr = m2 > 0 ? vo / m2 : 1.0;
    bwRows.push(`${SKU_L[k]},${kurt.toFixed(6)},${acf1.toFixed(6)},${bwr.toFixed(6)}`);
  }
  wf(outDir, "s21_bullwhip.csv", "sku,kurtosis,acf_lag1,bullwhip_ratio\n" + bwRows.join("\n") + "\n");
  wf(outDir, "s22_service.csv", "store,req,fill,fill_rate\n" + Array.from({ length: 10 }, (_, i) => `${i},${s.dailyFillStore[i][0]},${s.dailyFillStore[i][1]},${(s.dailyFillStore[i][0] > 0 ? s.dailyFillStore[i][1] / s.dailyFillStore[i][0] : 1.0).toFixed(6)}`).join("\n") + "\n");
  wf(outDir, "s23_leadtime_impact.csv", "idx,leadtime,qty\n" + s.leadtimePairs.map((p, i) => `${i},${p[0]},${p[1]}`).join("\n") + "\n");
  const imbM = s.kitImbSeries.length ? s.kitImbSeries.reduce((a, b) => a + b, 0) / s.kitImbSeries.length : 0.0;
  wf(outDir, "s24_kit_coupling.csv", `mean_imbalance,total_kits\n${imbM.toFixed(4)},${s.totalKits}\n`);
  const shkRows: string[] = [];
  for (let cl = 0; cl < 3; cl++) {
    const valsR: number[] = [];
    for (const [b, cev, st] of s.shockEvs) {
      if (cev !== cl) continue;
      let pre = 0, post = 0, pc = 0, qc = 0;
      for (const db of [-2, -1]) if (b + db >= 0) { pre += s.binStoreDem[st][b + db]; pc++; }
      for (const db of [0, 1]) if (b + db < 100) { post += s.binStoreDem[st][b + db]; qc++; }
      if (pc > 0 && qc > 0 && pre > 0) valsR.push((post / qc - pre / pc) / (pre / pc));
    }
    const cnt = valsR.length, mn = cnt > 0 ? valsR.reduce((a, b) => a + b, 0) / cnt : 0.0;
    const se = cnt > 1 ? Math.sqrt(Math.max(0, (valsR.reduce((a, b) => a + b * b, 0) / cnt - mn * mn) / (cnt - 1))) : 0.0;
    shkRows.push(`${cl},${cnt},${mn.toFixed(6)},${se.toFixed(6)}`);
  }
  wf(outDir, "s25_shock_event_study.csv", "cls,count,mean_rel_dev,se_rel_dev\n" + shkRows.join("\n") + "\n");
  wf(outDir, "s26_disruptions.csv", `disrupt_hits\n${s.branches[11]}\n`);
  wf(outDir, "s27_fees.csv", `n_fee_entries,fee_pool_cash\n${s.fees.length},${s.acc[ACC_FEE].cash}\n`);
  wf(outDir, "s28_capacity_rejects.csv", `wh_cap_reject,carrier_cap_reject\n${s.branches[4]},${s.branches[5]}\n`);

  const sw1 = ["S_lvl,rep,fill_rate,avg_oh"];
  for (const sv of [80, 220, 350]) for (let rp = 0; rp < 5; rp++) {
    const ss = new Sim({ ...c, master_seed: c.master_seed + rp + 1 }, false, ""); simRunDays(ss, 2, 2000, sv, -1, 1);
    const oh = ss.acc.slice(0, 10).reduce((a, x) => a + x.onHand[0], 0);
    const fr = ss.demUnitsTot > 0 ? ss.demUnitsFilled / ss.demUnitsTot : 1.0;
    sw1.push(`${sv},${rp},${fr.toFixed(6)},${(oh / 10).toFixed(2)}`);
  }
  wf(outDir, "s29_sweep_safety_stock.csv", sw1.join("\n") + "\n");
  const sw2 = ["disrupt_pct,rep,produced,sold"];
  for (const dp of [0, 4, 12]) for (let rp = 0; rp < 5; rp++) {
    const ss = new Sim({ ...c, master_seed: c.master_seed + rp + 10 }, false, ""); simRunDays(ss, 2, 2000, 0, dp, 1);
    sw2.push(`${dp},${rp},${ss.totalProduced.reduce((a, b) => a + b, 0)},${ss.totalSold.reduce((a, b) => a + b, 0)}`);
  }
  wf(outDir, "s30_sweep_shock_disrupt.csv", sw2.join("\n") + "\n");
  const sw3 = ["transship_on,rep,mean_imbalance,n_transships"];
  for (const ts of [0, 1]) for (let rp = 0; rp < 5; rp++) {
    const ss = new Sim({ ...c, master_seed: c.master_seed + rp + 20 }, false, ""); simRunDays(ss, 2, 2000, 0, -1, ts);
    const im = ss.kitImbSeries.length ? ss.kitImbSeries.reduce((a, b) => a + b, 0) / ss.kitImbSeries.length : 0.0;
    sw3.push(`${ts},${rp},${im.toFixed(4)},${ss.nTransships}`);
  }
  wf(outDir, "s31_sweep_transship.csv", sw3.join("\n") + "\n");

  const bn = ["DEMAND_FILL","DEMAND_BO","DEMAND_LOST","STORE_REPLEN","WH_CAP_REJECT","CARRIER_CAP_REJECT","ALLOC_FULL","ALLOC_PARTIAL","KIT_FULL","KIT_PARTIAL","TRANSSHIP","DISRUPT_HIT","SHOCK_HIT","CREDIT_HOLD","LIQUIDATE","CARRIER_DELAY"];
  wf(outDir, "s32_branches.csv", "branch,count\n" + bn.map((nm, i) => `${nm},${s.branches[i]}`).join("\n") + "\n");

  const sp1 = new Sim(c, false, "", c.parity_events_per_day); simRunDays(sp1, c.parity_days, c.parity_events_per_day, 0, -1, 1);
  const sp2 = new Sim(c, false, "", c.parity_events_per_day); simRunDays(sp2, c.parity_days, c.parity_events_per_day, 0, -1, 1);
  let h1 = 0xcbf29ce484222325n, h2 = 0xcbf29ce484222325n;
  for (let i = 0; i < NUM_ACCOUNTS; i++) {
    h1 = fnv1aI64(h1, sp1.acc[i].cash); h2 = fnv1aI64(h2, sp2.acc[i].cash);
    for (let k = 0; k < 4; k++) { h1 = fnv1aI64(h1, sp1.acc[i].onHand[k]); h2 = fnv1aI64(h2, sp2.acc[i].onHand[k]); }
  }
  wf(outDir, "s33_replay.csv", `ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n${hex64(sp1.eventHash)},${hex64(sp2.eventHash)},${hex64(h1)},${hex64(h2)},${sp1.eventHash === sp2.eventHash && h1 === h2 ? 1 : 0}\n`);
  wf(outDir, "s34_parity.csv", `parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n${c.parity_days},${c.parity_events_per_day},${hex64(sp1.eventHash)},${hex64(h1)},${sp1.totalSold.reduce((a, b) => a + b, 0)},${sp1.totalKits},${sp1.invViolations}\n`);

  const vals: number[] = [];
  for (let rp = 0; rp < 5; rp++) { const sb = new Sim(c, false, ""); const t0 = process.hrtime.bigint(); simRunDays(sb, 1, 8000, 0, -1, 1); vals.push(Number(process.hrtime.bigint() - t0) / 8000); }
  vals.sort((a, b) => a - b); const med = vals[2], dev = vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b);
  wf(outDir, "s35_bench.csv", `events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n8000,${med.toFixed(4)},${dev[2].toFixed(4)},${(1e9 / med).toFixed(2)}\n`);
  wf(outDir, "s36_sizes.csv", `n_sales,n_ships,n_prods,n_kits,n_fees\n${s.sales.length},${s.ships.length},${s.prods.length},${s.kits.length},${s.fees.length}\n`);

  const raw = s.sales.slice(0, 4096).map((r) => r[4]), nq = raw.length, rle: number[] = [];
  for (let i = 0; i < nq; ) { const dv = i === 0 ? raw[0] : raw[i] - raw[i - 1]; let run = 1; while (i + run < nq && raw[i + run] - raw[i + run - 1] === dv) run++; rle.push(dv, run); i += run; }
  const dec: number[] = [];
  for (let j = 0; j < rle.length; j += 2) for (let r = 0; r < rle[j + 1]; r++) dec.push(dec.length ? dec[dec.length - 1] + rle[j] : rle[j]);
  wf(outDir, "s37_compress.csv", `n_items,rle_ints,ratio,roundtrip_ok\n${nq},${rle.length},${(rle.length ? nq / rle.length : 1).toFixed(4)},${dec.length === nq && dec.every((v, i) => v === raw[i]) ? 1 : 0}\n`);

  const si = new Sim(c, false, "", s.capScale * 2000), lc = si.acc.map((a) => a.cash), loh = si.acc.map((a) => [...a.onHand]), lit = Array.from({ length: NUM_ACCOUNTS }, () => [0, 0, 0, 0]);
  for (const r of s.sales) { lc[ACC_CUST] -= r[5]; lc[r[2]] += r[5]; loh[r[2]][r[3]] -= r[4]; }
  for (const r of s.ships) {
    if (r[9] === 1) { loh[r[2]][r[5]] -= r[6]; lit[r[3]][r[5]] += r[6]; lc[r[3]] -= r[7] + r[8]; lc[r[2]] += r[7]; lc[17 + r[4]] += r[8]; }
    else { lit[r[3]][r[5]] -= r[6]; loh[r[3]][r[5]] += r[6]; }
  }
  for (const r of s.prods) { lc[r[2]] -= r[5]; lc[ACC_SUPP] += r[5]; loh[r[2]][r[3]] += r[4]; }
  for (const r of s.kits) { lc[r[2]] -= r[4]; lc[ACC_FEE] += r[4]; for (let k = 0; k < 3; k++) loh[r[2]][k] -= r[3]; loh[r[2]][3] += r[3]; }
  for (const r of s.fees) { lc[r[2]] -= r[4]; lc[r[3]] += r[4]; }
  let mcd = 0, mid = 0;
  for (let i = 0; i < NUM_ACCOUNTS; i++) { mcd = Math.max(mcd, Math.abs(lc[i] - s.acc[i].cash)); for (let k = 0; k < 4; k++) mid = Math.max(mid, Math.abs(loh[i][k] - s.acc[i].onHand[k]) + Math.abs(lit[i][k] - s.acc[i].inTransit[k])); }
  wf(outDir, "s38_ledger_replay.csv", `max_cash_diff,max_inv_diff,match\n${mcd},${mid},${mcd === 0 && mid === 0 ? 1 : 0}\n`);
}
