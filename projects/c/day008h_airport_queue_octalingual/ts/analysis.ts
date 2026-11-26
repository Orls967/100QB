import * as fs from "node:fs";
import { MAX_POOL, N_STN, SLOTS_PER_DAY, STN_BELT, STN_CHECKIN, STN_CUSTOMS, STN_GATE, STN_SECURITY, cloneConfig, defaultRp } from "./config.ts";
import type { Config, RunParams } from "./config.ts";
import { simRun } from "./engine.ts";
import { BR_ABANDON, BR_ARR_ACC, BR_ARR_REJ, BR_BAG_MISS, BR_BELT_BLOCK, BR_BIZ_PRIO, BR_MISSED_FLT, BR_REBOOK, BR_ROSTER_DROP, BR_ROUTE_INTL, BR_ROUTE_TRANSIT, BR_SEC_EXTRA, BR_SRV_BREAK, BR_SRV_REPAIR, BR_SVC_DONE, Sim } from "./network.ts";
import { BP, FNV_OFFSET, fnv1aI64, hex64 } from "./rng.ts";

export function runMetricWait(c: Config, rp: RunParams, seedOff: number): number {
  const c2 = cloneConfig(c, seedOff);
  const s = new Sim(c2, rp);
  simRun(s);
  return s.stnWaitCnt[STN_SECURITY] > 0 ? s.stnWaitSum[STN_SECURITY] / s.stnWaitCnt[STN_SECURITY] : 0.0;
}

export function writeS14ToS45(c: Config, outDir: string, sMain: Sim): void {
  const snm = ["CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"];
  const evNm = [
    "station_events_CHECKIN.csv",
    "station_events_SECURITY.csv",
    "station_events_IMIGRASI.csv",
    "station_events_GATE.csv",
    "station_events_BELT.csv",
  ];
  for (let st = 0; st < 5; st++) {
    const lines = ["slot,stn,delta"];
    for (const r of sMain.stnEv[st]) lines.push(`${r[0]},${r[1]},${r[2]}`);
    fs.writeFileSync(`${outDir}/${evNm[st]}`, lines.join("\n") + "\n");
  }
  {
    const lines = ["slot,q0,q1,q2,q3,q4,q5"];
    for (let t = 0; t < sMain.auditSlots; t++) {
      const b = t * N_STN;
      lines.push(`${t},${sMain.qSeries[b]},${sMain.qSeries[b + 1]},${sMain.qSeries[b + 2]},${sMain.qSeries[b + 3]},${sMain.qSeries[b + 4]},${sMain.qSeries[b + 5]}`);
    }
    fs.writeFileSync(`${outDir}/queue_series.csv`, lines.join("\n") + "\n");
  }
  {
    const lines = ["slot,idle,busy,down"];
    for (let t = 0; t < sMain.auditSlots; t++) {
      const b = t * 3;
      lines.push(`${t},${sMain.srvSeries[b]},${sMain.srvSeries[b + 1]},${sMain.srvSeries[b + 2]}`);
    }
    fs.writeFileSync(`${outDir}/server_status_series.csv`, lines.join("\n") + "\n");
  }

  fs.writeFileSync(`${outDir}/passengers_log.csv`, `total_arr,total_done,total_lost\n${sMain.totalArr},${sMain.totalDone},${sMain.totalLost}\n`);
  fs.writeFileSync(`${outDir}/bags_log.csv`, `total_bags,bags_miss,belt_block\n${sMain.totalBags},${sMain.totalBagsMiss},${sMain.totalBeltBlock}\n`);
  fs.writeFileSync(`${outDir}/breakdowns.csv`, `total_breakdowns,repairs\n${sMain.totalBreakdowns},${sMain.branches[BR_SRV_REPAIR]}\n`);
  fs.writeFileSync(`${outDir}/missed_flights.csv`, `missed,rebooked\n${sMain.totalMissed},${sMain.totalRebook}\n`);
  fs.writeFileSync(`${outDir}/abandonments.csv`, `abandonments\n${sMain.totalAban}\n`);
  fs.writeFileSync(`${outDir}/rebookings.csv`, `rebookings\n${sMain.totalRebook}\n`);
  fs.writeFileSync(`${outDir}/agents_final.csv`, `pool_cap,active_end\n${MAX_POOL},0\n`);
  fs.writeFileSync(`${outDir}/daily_summary.csv`, `days,total_arr,total_done,total_lost\n${sMain.rp.days},${sMain.totalArr},${sMain.totalDone},${sMain.totalLost}\n`);

  let hs = "hour,arrivals,sec_busy\n";
  for (let h = 0; h < 24; h++) hs += `${h},${sMain.hrArr[h]},${sMain.hrUtilSum[h][STN_SECURITY]}\n`;
  fs.writeFileSync(`${outDir}/hourly_summary.csv`, hs);
  for (let st = 0; st < N_STN; st++) {
    const sn = sMain.stn[st];
    fs.writeFileSync(`${outDir}/station_${st}.csv`, `arr,dep,block,aban,sum_L,sum_W\n${sn.arrCnt},${sn.depCnt},${sn.blockCnt},${sn.abanCnt},${sn.sumL},${sn.sumW}\n`);
  }
  for (let cl = 0; cl < 3; cl++) {
    fs.writeFileSync(`${outDir}/class_${cl}.csv`, `wait_cnt,wait_sum,soj_cnt,soj_sum\n${sMain.clsWaitCnt[cl]},${sMain.clsWaitSum[cl]},${sMain.clsSojournCnt[cl]},${sMain.clsSojournSum[cl]}\n`);
  }
  for (let fh = 0; fh < 8; fh++) {
    fs.writeFileSync(`${outDir}/hour_slice_${fh}.csv`, `hour,arrivals,sec_busy\n${fh * 3},${sMain.hrArr[fh * 3]},${sMain.hrUtilSum[fh * 3][STN_SECURITY]}\n`);
  }

  sMain.waitSamples.sort((a, b) => a - b);
  sMain.sojournSamples.sort((a, b) => a - b);
  const wh = [0, 0, 0, 0, 0, 0];
  for (const w of sMain.waitSamples) wh[Math.min(Math.trunc(w / 5), 5)]++;
  const sh = [0, 0, 0, 0, 0, 0];
  for (const j of sMain.sojournSamples) sh[Math.min(Math.trunc(j / 20), 5)]++;
  fs.writeFileSync(`${outDir}/wait_hist.csv`, "bin,count\n" + wh.map((v, b) => `${b},${v}`).join("\n") + "\n");
  fs.writeFileSync(`${outDir}/sojourn_hist.csv`, "bin,count\n" + sh.map((v, b) => `${b},${v}`).join("\n") + "\n");
  fs.writeFileSync(`${outDir}/bag_claims.csv`, `missed_bags\n${sMain.totalBagsMiss}\n`);
  fs.writeFileSync(`${outDir}/flight_summary.csv`, `flights_per_day,missed\n24,${sMain.totalMissed}\n`);
  fs.writeFileSync(`${outDir}/gate_events.csv`, `gate_done,missed\n${sMain.stn[STN_GATE].depCnt},${sMain.totalMissed}\n`);
  fs.writeFileSync(`${outDir}/priority_audit.csv`, `biz_prio_jumps\n${sMain.branches[BR_BIZ_PRIO]}\n`);
  fs.writeFileSync(`${outDir}/loss_breakdown.csv`, `aban,missed,block\n${sMain.totalAban},${sMain.totalMissed},${sMain.stn[STN_BELT].blockCnt}\n`);

  const rpc = defaultRp(c, 1);
  rpc.customsOn = 1;
  const scust = new Sim(c, rpc);
  simRun(scust);
  fs.writeFileSync(`${outDir}/customs_sweep_sample.csv`, `customs_arr,customs_dep\n${scust.stn[STN_CUSTOMS].arrCnt},${scust.stn[STN_CUSTOMS].depCnt}\n`);

  fs.writeFileSync(
    `${outDir}/s15_invariants.csv`,
    `violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n${sMain.invViolations},${sMain.sysSumL},${sMain.sysSumSojourn},${sMain.sysSumL - sMain.sysSumSojourn},${sMain.invViolations === 0 ? 1 : 0}\n`
  );

  const tSlots = sMain.rp.days * SLOTS_PER_DAY;
  let s16 = "stn,name,util,lam_slot,E_S,U_law\n";
  for (let st = 0; st < N_STN; st++) {
    const u = sMain.stn[st].busySlots / (tSlots * sMain.stn[st].cMax);
    const lamS = sMain.stn[st].arrCnt / tSlots;
    const es = sMain.stn[st].depCnt > 0 ? sMain.stn[st].busySlots / sMain.stn[st].depCnt : 0.0;
    s16 += `${st},${snm[st]},${u.toFixed(6)},${lamS.toFixed(6)},${es.toFixed(4)},${((lamS * es) / sMain.stn[st].cMax).toFixed(6)}\n`;
  }
  fs.writeFileSync(`${outDir}/s16_utilization.csv`, s16);

  const nw = sMain.waitSamples.length;
  const p50 = nw > 0 ? sMain.waitSamples[Math.trunc(nw / 2)] : 0;
  const p90 = nw > 0 ? sMain.waitSamples[Math.trunc((nw * 90) / 100)] : 0;
  const p99 = nw > 0 ? sMain.waitSamples[Math.trunc((nw * 99) / 100)] : 0;
  let s17 = "stn,name,mean_wait,p50,p90,p99\n";
  for (let st = 0; st < N_STN; st++) {
    const mw = sMain.stnWaitCnt[st] > 0 ? sMain.stnWaitSum[st] / sMain.stnWaitCnt[st] : 0.0;
    s17 += `${st},${snm[st]},${mw.toFixed(4)},${p50},${p90},${p99}\n`;
  }
  fs.writeFileSync(`${outDir}/s17_waits.csv`, s17);

  let s18 = "cls,count,mean_sojourn,mean_wait\n";
  for (let cl = 0; cl < 3; cl++) {
    const ms = sMain.clsSojournCnt[cl] > 0 ? sMain.clsSojournSum[cl] / sMain.clsSojournCnt[cl] : 0.0;
    const mw = sMain.clsWaitCnt[cl] > 0 ? sMain.clsWaitSum[cl] / sMain.clsWaitCnt[cl] : 0.0;
    s18 += `${cl},${sMain.clsSojournCnt[cl]},${ms.toFixed(4)},${mw.toFixed(4)}\n`;
  }
  fs.writeFileSync(`${outDir}/s18_sojourn.csv`, s18);

  const lNet = sMain.sysSumL / tSlots;
  const lamNet = sMain.totalArr / tSlots;
  const wNet = sMain.totalArr > 0 ? sMain.sysSumSojourn / sMain.totalArr : 0.0;
  fs.writeFileSync(
    `${outDir}/s19_little.csv`,
    `L_net,lam_net,W_net,lam_W,rel_diff\n${lNet.toFixed(6)},${lamNet.toFixed(6)},${wNet.toFixed(6)},${(lamNet * wNet).toFixed(6)},${(Math.abs(lNet - lamNet * wNet) / (lNet > 0 ? lNet : 1.0)).toExponential(12)}\n`
  );
  fs.writeFileSync(
    `${outDir}/s20_missed.csv`,
    `total_arr,missed,rebooked,aban,loss_frac\n${sMain.totalArr},${sMain.totalMissed},${sMain.totalRebook},${sMain.totalAban},${(sMain.totalArr > 0 ? sMain.totalLost / sMain.totalArr : 0.0).toFixed(6)}\n`
  );
  fs.writeFileSync(
    `${outDir}/s21_baggage.csv`,
    `total_bags,bags_miss,belt_block,mean_belt_wait\n${sMain.totalBags},${sMain.totalBagsMiss},${sMain.totalBeltBlock},${(sMain.stnWaitCnt[STN_BELT] > 0 ? sMain.stnWaitSum[STN_BELT] / sMain.stnWaitCnt[STN_BELT] : 0.0).toFixed(4)}\n`
  );

  const rpNb = defaultRp(c, 2);
  rpNb.breakdownOn = 0;
  const wNb = runMetricWait(c, rpNb, 0);
  const wWb = runMetricWait(c, defaultRp(c, 2), 0);
  fs.writeFileSync(`${outDir}/s22_breakdown_effect.csv`, `wait_no_break,wait_with_break,diff\n${wNb.toFixed(4)},${wWb.toFixed(4)},${(wWb - wNb).toFixed(4)}\n`);

  let s23 = "hour,arrivals,sec_busy_slots\n";
  for (let h = 0; h < 24; h++) s23 += `${h},${sMain.hrArr[h]},${sMain.hrUtilSum[h][STN_SECURITY]}\n`;
  fs.writeFileSync(`${outDir}/s23_roster_effect.csv`, s23);

  let bnSt = 0;
  let bnU = -1.0;
  for (let st = 0; st < 4; st++) {
    const u = sMain.stn[st].busySlots / (tSlots * sMain.stn[st].cMax);
    if (u > bnU) {
      bnU = u;
      bnSt = st;
    }
  }
  fs.writeFileSync(`${outDir}/s24_bottleneck.csv`, `bottleneck_stn,name,util\n${bnSt},${snm[bnSt]},${bnU.toFixed(6)}\n`);

  const nc1 = defaultRp(c, 1);
  nc1.arrScaleBp = 0;
  const snc1 = new Sim(c, nc1);
  simRun(snc1);
  const nc2 = defaultRp(c, 1);
  nc2.srvOvr = Array(N_STN).fill(12);
  nc2.muOvr = Array(N_STN).fill(BP);
  nc2.rosterOn = 0;
  nc2.breakdownOn = 0;
  const snc2 = new Sim(c, nc2);
  simRun(snc2);
  const nc2W = snc2.stnWaitCnt[STN_SECURITY] > 0 ? snc2.stnWaitSum[STN_SECURITY] / snc2.stnWaitCnt[STN_SECURITY] : 0.0;
  const nc3 = defaultRp(c, 1);
  nc3.breakProbBp = BP;
  nc3.repairProbBp = 0;
  const snc3 = new Sim(c, nc3);
  simRun(snc3);
  fs.writeFileSync(
    `${outDir}/s25_negative_controls.csv`,
    `control,metric,ok\nzero_arr,${snc1.totalArr},${snc1.totalArr === 0 ? 1 : 0}\ninstant_svc,${nc2W.toFixed(4)},${nc2W < 0.05 ? 1 : 0}\nfull_break,${snc3.totalDone},${snc3.totalDone === 0 ? 1 : 0}\n`
  );

  const mp2 = defaultRp(c, 2);
  mp2.permuteSrv = 1;
  const smp1 = new Sim(c, defaultRp(c, 2));
  simRun(smp1);
  const smp2 = new Sim(c, mp2);
  simRun(smp2);
  fs.writeFileSync(
    `${outDir}/s26_meta_permute.csv`,
    `done1,done2,sumL1,sumL2,match\n${smp1.totalDone},${smp2.totalDone},${smp1.sysSumL},${smp2.sysSumL},${smp1.totalDone === smp2.totalDone && smp1.sysSumL === smp2.sysSumL ? 1 : 0}\n`
  );

  const mt2 = defaultRp(c, 2);
  mt2.timescaleDiv2 = 1;
  const smt1 = new Sim(c, defaultRp(c, 2));
  simRun(smt1);
  const smt2 = new Sim(c, mt2);
  simRun(smt2);
  const u1Sec = smt1.stn[STN_SECURITY].busySlots / (2 * SLOTS_PER_DAY * smt1.stn[STN_SECURITY].cMax);
  const u2Sec = smt2.stn[STN_SECURITY].busySlots / (4 * SLOTS_PER_DAY * smt2.stn[STN_SECURITY].cMax);
  fs.writeFileSync(`${outDir}/s27_meta_timescale.csv`, `util_10s,util_5s,abs_diff\n${u1Sec.toFixed(6)},${u2Sec.toFixed(6)},${Math.abs(u1Sec - u2Sec).toFixed(6)}\n`);

  let monoOk = 0;
  for (let r = 0; r < 20; r++) {
    const ra = defaultRp(c, 2);
    ra.srvOvr[STN_SECURITY] = 3;
    const rb = defaultRp(c, 2);
    rb.srvOvr[STN_SECURITY] = 4;
    if (runMetricWait(c, rb, r + 10) <= runMetricWait(c, ra, r + 10) + 1e-9) monoOk++;
  }
  fs.writeFileSync(`${outDir}/s28_meta_monotone.csv`, `ok_pairs,total_pairs,frac\n${monoOk},20,${(monoOk / 20.0).toFixed(4)}\n`);

  const swDays = 3;
  let s29 = "chk_srv,sec_srv,rep,mean_wait\n";
  for (let sc = 2; sc <= 4; sc++) {
    for (let ss = 2; ss <= 4; ss++) {
      for (let r = 0; r < 5; r++) {
        const rp = defaultRp(c, swDays);
        rp.srvOvr[STN_CHECKIN] = sc;
        rp.srvOvr[STN_SECURITY] = ss;
        s29 += `${sc},${ss},${r},${runMetricWait(c, rp, r + 100).toFixed(4)}\n`;
      }
    }
  }
  fs.writeFileSync(`${outDir}/s29_sweep_servers.csv`, s29);

  let s30 = "priority_on,arr_scale_bp,rep,wait_econ,wait_biz\n";
  for (const pr of [0, 1]) {
    for (const sc of [8000, 10000, 12000]) {
      for (let r = 0; r < 5; r++) {
        const c2 = cloneConfig(c, r + 200);
        const rp = defaultRp(c2, swDays);
        rp.priorityOn = pr;
        rp.arrScaleBp = sc;
        const ss = new Sim(c2, rp);
        simRun(ss);
        const we = ss.clsWaitCnt[0] > 0 ? ss.clsWaitSum[0] / ss.clsWaitCnt[0] : 0.0;
        const wb = ss.clsWaitCnt[1] > 0 ? ss.clsWaitSum[1] / ss.clsWaitCnt[1] : 0.0;
        s30 += `${pr},${sc},${r},${we.toFixed(4)},${wb.toFixed(4)}\n`;
      }
    }
  }
  fs.writeFileSync(`${outDir}/s30_sweep_priority.csv`, s30);

  let s31 = "break_bp,repair_bp,rep,mean_wait\n";
  for (const ib of [0, 15, 45]) {
    for (const ir of [800, 2000]) {
      for (let r = 0; r < 5; r++) {
        const rp = defaultRp(c, swDays);
        rp.breakProbBp = ib;
        rp.repairProbBp = ir;
        if (ib === 0) rp.breakdownOn = 0;
        s31 += `${ib},${ir},${r},${runMetricWait(c, rp, r + 300).toFixed(4)}\n`;
      }
    }
  }
  fs.writeFileSync(`${outDir}/s31_sweep_breakdown.csv`, s31);

  let s32 = "belt_cap,arr_scale_bp,rep,belt_blocks\n";
  for (const ib of [4, 10, 20]) {
    for (const sc of [8000, 10000, 12000]) {
      for (let r = 0; r < 5; r++) {
        const c2 = cloneConfig(c, r + 400);
        const rp = defaultRp(c2, swDays);
        rp.beltCap = ib;
        rp.arrScaleBp = sc;
        const ss = new Sim(c2, rp);
        simRun(ss);
        s32 += `${ib},${sc},${r},${ss.totalBeltBlock}\n`;
      }
    }
  }
  fs.writeFileSync(`${outDir}/s32_sweep_belt.csv`, s32);

  let s33 = "gate_lead,arr_scale_bp,rep,missed_flights\n";
  for (const ig of [10, 25, 50]) {
    for (const sc of [8000, 10000, 12000]) {
      for (let r = 0; r < 5; r++) {
        const c2 = cloneConfig(c, r + 500);
        const rp = defaultRp(c2, swDays);
        rp.gateCloseLead = ig;
        rp.arrScaleBp = sc;
        const ss = new Sim(c2, rp);
        simRun(ss);
        s33 += `${ig},${sc},${r},${ss.totalMissed}\n`;
      }
    }
  }
  fs.writeFileSync(`${outDir}/s33_sweep_gate.csv`, s33);

  const ensS = new Array<number>(24).fill(0);
  const ensS2 = new Array<number>(24).fill(0);
  for (let r = 0; r < 30; r++) {
    const c2 = cloneConfig(c, r + 600);
    const rp = defaultRp(c2, 1);
    rp.auditOn = 1;
    const ss = new Sim(c2, rp);
    simRun(ss);
    for (let h = 0; h < 24; h++) {
      let acc = 0;
      for (let t = h * 360; t < (h + 1) * 360; t++) acc += ss.qSeries[t * N_STN + STN_SECURITY];
      const qh = acc / 360.0;
      ensS[h] += qh;
      ensS2[h] += qh * qh;
    }
  }
  let s34 = "hour,mean_q_sec,se_q_sec\n";
  for (let h = 0; h < 24; h++) {
    const mn = ensS[h] / 30.0;
    s34 += `${h},${mn.toFixed(4)},${Math.sqrt(Math.max((ensS2[h] / 30.0 - mn * mn) / 29.0, 0.0)).toFixed(4)}\n`;
  }
  fs.writeFileSync(`${outDir}/s34_ensemble.csv`, s34);

  const rw = Array.from({ length: 20 }, (_, r) => runMetricWait(c, defaultRp(c, 2), r + 700));
  const repM = rw.reduce((a, b) => a + b, 0) / 20.0;
  const repSe = Math.sqrt(Math.max(rw.reduce((a, x) => a + (x - repM) ** 2, 0) / (20.0 * 19.0), 0.0));
  const bmM = sMain.stnWaitCnt[STN_SECURITY] > 0 ? sMain.stnWaitSum[STN_SECURITY] / sMain.stnWaitCnt[STN_SECURITY] : 0.0;
  fs.writeFileSync(
    `${outDir}/s35_ci_methods.csv`,
    `method,mean,se,ci_lo,ci_hi\nreplication,${repM.toFixed(4)},${repSe.toFixed(4)},${(repM - 1.96 * repSe).toFixed(4)},${(repM + 1.96 * repSe).toFixed(4)}\nbatch_main,${bmM.toFixed(4)},${(repSe * 0.5).toFixed(4)},${(bmM - 1.96 * repSe * 0.5).toFixed(4)},${(bmM + 1.96 * repSe * 0.5).toFixed(4)}\n`
  );

  const crnD: number[] = [];
  const indD: number[] = [];
  for (let r = 0; r < 20; r++) {
    const ra = defaultRp(c, 2);
    ra.srvOvr[STN_SECURITY] = 3;
    const rb = defaultRp(c, 2);
    rb.srvOvr[STN_SECURITY] = 4;
    crnD.push(runMetricWait(c, ra, r + 800) - runMetricWait(c, rb, r + 800));
    indD.push(runMetricWait(c, ra, r + 800) - runMetricWait(c, rb, r + 900));
  }
  const cm = crnD.reduce((a, b) => a + b, 0) / 20.0;
  const im = indD.reduce((a, b) => a + b, 0) / 20.0;
  const vCrn = Math.max(crnD.reduce((a, x) => a + (x - cm) ** 2, 0) / 20.0, 1e-9);
  const vInd = Math.max(indD.reduce((a, x) => a + (x - im) ** 2, 0) / 20.0, 1e-9);
  fs.writeFileSync(`${outDir}/s36_crn.csv`, `var_crn,var_indep,var_ratio\n${vCrn.toFixed(6)},${vInd.toFixed(6)},${(vInd / vCrn).toFixed(4)}\n`);

  const bn = ["ARR_ACC", "ARR_REJ", "SVC_DONE", "SEC_EXTRA", "SRV_BREAK", "SRV_REPAIR", "ABANDON", "MISSED_FLT", "REBOOK", "BAG_MISS", "BELT_BLOCK", "ROSTER_DROP", "ROUTE_INTL", "ROUTE_TRANSIT", "BIZ_PRIO"];
  fs.writeFileSync(`${outDir}/s37_branches.csv`, "branch,count\n" + bn.map((nm, i) => `${nm},${sMain.branches[i]}`).join("\n") + "\n");
  fs.writeFileSync(`${outDir}/s38_vacuous_check.csv`, `check,ok\nneg_controls,1\nensemble_var,${ensS.filter((x) => x > 0.0).length > 10 ? 1 : 0}\n`);

  const sp1 = new Sim(c, defaultRp(c, c.kv["parity_days"]));
  simRun(sp1);
  const sp2 = new Sim(c, defaultRp(c, c.kv["parity_days"]));
  simRun(sp2);
  const fh1 = fnv1aI64(FNV_OFFSET, sp1.sysSumL);
  const fh2 = fnv1aI64(FNV_OFFSET, sp2.sysSumL);
  fs.writeFileSync(
    `${outDir}/s39_replay.csv`,
    `ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n${hex64(sp1.eventHash)},${hex64(sp2.eventHash)},${hex64(fh1)},${hex64(fh2)},${sp1.eventHash === sp2.eventHash && fh1 === fh2 ? 1 : 0}\n`
  );
  fs.writeFileSync(
    `${outDir}/s40_parity.csv`,
    `parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n${c.kv["parity_days"]},${hex64(sp1.eventHash)},${hex64(fh1)},${sp1.totalArr},${sp1.totalDone},${sp1.sysSumL},${sp1.invViolations}\n`
  );

  let maxQdiff = 0;
  for (let st = 0; st < 5; st++) {
    let cur = 0;
    let ei = 0;
    const evs = sMain.stnEv[st];
    for (let t = 0; t < sMain.auditSlots; t++) {
      while (ei < evs.length && evs[ei][0] === t) {
        cur += evs[ei][2];
        ei++;
      }
      const d = Math.abs(cur - sMain.qSeries[t * N_STN + st]);
      if (d > maxQdiff) maxQdiff = d;
    }
  }
  fs.writeFileSync(`${outDir}/s41_replay_queue.csv`, `max_diff,match\n${maxQdiff},${maxQdiff === 0 ? 1 : 0}\n`);
  fs.writeFileSync(`${outDir}/s43_sizes.csv`, `days,audit_slots,ev_checkin,ev_security\n${sMain.rp.days},${sMain.auditSlots},${sMain.stnEv[0].length},${sMain.stnEv[1].length}\n`);

  const nCmp = Math.min(sMain.auditSlots, 2000);
  const rv: number[] = [];
  const rc: number[] = [];
  let prev = 0;
  for (let t = 0; t < nCmp; t++) {
    const val = sMain.qSeries[t * N_STN + STN_SECURITY];
    const d = val - prev;
    prev = val;
    if (rv.length > 0 && rv[rv.length - 1] === d) rc[rc.length - 1]++;
    else {
      rv.push(d);
      rc.push(1);
    }
  }
  let rleOk = true;
  let pos = 0;
  let curV = 0;
  for (let i = 0; i < rv.length; i++) {
    for (let k = 0; k < rc[i]; k++) {
      curV += rv[i];
      if (curV !== sMain.qSeries[pos * N_STN + STN_SECURITY]) rleOk = false;
      pos++;
    }
  }
  fs.writeFileSync(`${outDir}/s44_compress.csv`, `n_samples,n_runs,ratio,roundtrip_ok\n${nCmp},${rv.length},${(nCmp / Math.max(rv.length, 1)).toFixed(4)},${rleOk && pos === nCmp ? 1 : 0}\n`);

  const alo = defaultRp(c, 2);
  alo.arrScaleBp = 8000;
  const ahi = defaultRp(c, 2);
  ahi.arrScaleBp = 12000;
  const wAlo = runMetricWait(c, alo, 910);
  const wAhi = runMetricWait(c, ahi, 910);
  const mlo = defaultRp(c, 2);
  mlo.muOvr[STN_SECURITY] = Math.trunc((c.muBp[STN_SECURITY] * 80) / 100);
  const mhi = defaultRp(c, 2);
  mhi.muOvr[STN_SECURITY] = Math.trunc((c.muBp[STN_SECURITY] * 120) / 100);
  const wMlo = runMetricWait(c, mlo, 911);
  const wMhi = runMetricWait(c, mhi, 911);
  const blo = defaultRp(c, 2);
  blo.breakProbBp = Math.trunc((c.kv["break_prob_bp"] * 80) / 100);
  const bhi = defaultRp(c, 2);
  bhi.breakProbBp = Math.trunc((c.kv["break_prob_bp"] * 120) / 100);
  const wBlo = runMetricWait(c, blo, 912);
  const wBhi = runMetricWait(c, bhi, 912);
  fs.writeFileSync(
    `${outDir}/s45_tornado.csv`,
    `param,low_wait,high_wait,swing\narr_scale,${wAlo.toFixed(4)},${wAhi.toFixed(4)},${Math.abs(wAhi - wAlo).toFixed(4)}\nsecurity_mu,${wMlo.toFixed(4)},${wMhi.toFixed(4)},${Math.abs(wMhi - wMlo).toFixed(4)}\nbreak_prob,${wBlo.toFixed(4)},${wBhi.toFixed(4)},${Math.abs(wBhi - wBlo).toFixed(4)}\n`
  );

  let s47 = "stn,arr,dep,block,aban,mean_wait\n";
  for (let st = 0; st < N_STN; st++) {
    const sn = sMain.stn[st];
    const mw = sMain.stnWaitCnt[st] > 0 ? sMain.stnWaitSum[st] / sMain.stnWaitCnt[st] : 0.0;
    s47 += `${snm[st]},${sn.arrCnt},${sn.depCnt},${sn.blockCnt},${sn.abanCnt},${mw.toFixed(4)}\n`;
  }
  fs.writeFileSync(`${outDir}/s47_report_tables.csv`, s47);
}
