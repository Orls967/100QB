import { writeFileSync } from "node:fs";
import { analyticalR0, simRun, solveFinalSize } from "./engine.ts";
import { BP, ST_D, Sim, defaultRP } from "./population.ts";
import type { Config } from "./population.ts";
import { FNV_OFFSET, MASK64, Xoshiro256, fnv1aBytes, fnv1aI64, fnv1aU32Pair, fnv1aU64, hex64, splitmix64 } from "./rng.ts";

function netDegrees(s: Sim): [number, number, number] {
  let shh = 0, ssw = 0, scom = 0;
  for (const p of s.pop) {
    const h = p.hh_id;
    shh += s.hh_off[h + 1] - s.hh_off[h] - 1;
    ssw += p.age_grp === 0 ? 6.0 : p.age_grp === 1 ? 5.0 : 0.0;
    scom += p.super_flag !== 0 ? 12.0 : 3.0;
  }
  const nf = s.n;
  return [shh / nf, ssw / nf, scom / nf];
}

export function writeS01ToS09(c: Config, outDir: string): void {
  let hrng = FNV_OFFSET;
  const pairs: [number, number][] = [[0, 0], [1, 2], [2, 1]];
  const expHv = [0xbe40b757fe0b6697n, 0x659394d6f2a43a95n, 0x3abda5bac0dd49cbn];
  let rngOk = true;
  for (let i = 0; i < 3; i++) {
    const x = (BigInt(c.master_seed) & MASK64) ^ fnv1aU32Pair(pairs[i][0], pairs[i][1]);
    const [, sd] = splitmix64(x);
    const st = new Xoshiro256(sd);
    const buf = new Uint8Array(32);
    const dv = new DataView(buf.buffer);
    for (let j = 0; j < 4; j++) dv.setBigUint64(j * 8, st.nextU64(), true);
    const hv = fnv1aBytes(buf);
    if (hv !== expHv[i]) rngOk = false;
    hrng = fnv1aU64(hrng, hv);
  }
  const st0 = new Sim(c, defaultRP(200, 10, 450, 0), false);
  simRun(st0);
  const invOk = st0.inv_violations === 0;
  writeFileSync(`${outDir}/s01_selftest.txt`, `ok,h_rng,inv_ok,tot_trans\n${rngOk && invOk ? 1 : 0},${hex64(hrng)},${invOk ? 1 : 0},${st0.total_trans}\n`);
  writeFileSync(`${outDir}/s02_config.csv`, `key,value\nmaster_seed,${c.master_seed}\nparity_N,${c.parity_N}\nparity_days,${c.parity_days}\nmain_days,${c.main_days}\nbeta_bp_init,${c.beta_bp_init}\n`);

  const sp = new Sim(c, defaultRP(1000, 1, 450, 0), false);
  writeFileSync(`${outDir}/s03_population.csv`, `age_grp,count\n0,${sp.age_tot[0]}\n1,${sp.age_tot[1]}\n2,${sp.age_tot[2]}\n`);
  const hhLines = ["hh_id,size"];
  for (let h = 0; h < sp.n_hh && h < 50; h++) hhLines.push(`${h},${sp.hh_off[h + 1] - sp.hh_off[h]}`);
  writeFileSync(`${outDir}/s03_households.csv`, hhLines.join("\n") + "\n");
  const swLines = ["sw_id,size"];
  for (let w = 0; w < sp.n_sw && w < 50; w++) swLines.push(`${w},${sp.sw_off[w + 1] - sp.sw_off[w]}`);
  writeFileSync(`${outDir}/s03_schools_work.csv`, swLines.join("\n") + "\n");
  const [chh, csw, ccom] = netDegrees(sp);
  writeFileSync(`${outDir}/s04_network.csv`, `n_hh,n_sw,deg_hh,deg_sw,deg_comm\n${sp.n_hh},${sp.n_sw},${chh.toFixed(4)},${csw.toFixed(4)},${ccom.toFixed(4)}\n`);

  const cLines = ["case_id,vax,mask,quar,eff_bp"];
  for (let k = 0; k < 12; k++) {
    const vax = k % 3, mask = Math.floor(k / 3) % 2, quar = Math.floor(k / 6) % 2;
    let p = quar !== 0 ? 0 : 450;
    if (mask !== 0) p = Math.floor((p * (BP - c.mask_eff_bp)) / BP);
    if (vax === 1) p = Math.floor((p * (BP - c.vax_eff1_bp)) / BP);
    else if (vax === 2) p = Math.floor((p * (BP - c.vax_eff2_bp)) / BP);
    cLines.push(`${k},${vax},${mask},${quar},${p}`);
  }
  writeFileSync(`${outDir}/s05_contact_script.csv`, cLines.join("\n") + "\n");

  const pLines = ["case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp"];
  for (let k = 0; k < 12; k++) {
    const ag = k % 3, rej = Math.floor(k / 3) % 2;
    const ph = ag === 0 ? 120 : ag === 1 ? 450 : 1800;
    let pd = ag === 0 ? 100 : ag === 1 ? 500 : 2000;
    if (rej !== 0) pd *= 2;
    pLines.push(`${k},${ag},${rej},${ph},${pd}`);
  }
  writeFileSync(`${outDir}/s06_progression_script.csv`, pLines.join("\n") + "\n");
  writeFileSync(`${outDir}/s07_interventions.csv`, `name,start_day,param_bp\nschool_close,${c.school_close_day},10000\nlockdown,${c.lockdown_day},5000\nmask,${c.mask_day},${c.mask_eff_bp}\ntest,${c.test_day},${c.test_prob_bp}\nvax,${c.vax_day},${c.vax_daily_bp}\nvariant,${c.variant_day},${c.variant_mult_bp}\n`);

  const calLines = ["beta_bp,r0_analytic,attack_rate,calibrated"];
  for (const bv of [300, 380, 450, 520, 600]) {
    const sc = new Sim(c, defaultRP(1000, 120, bv, 0), false);
    simRun(sc);
    const [ch, cs2, cc2] = netDegrees(sc);
    const r0 = analyticalR0(c, bv, ch, cs2, cc2);
    const ar = sc.total_trans / sc.n;
    calLines.push(`${bv},${r0.toFixed(4)},${ar.toFixed(4)},${r0 >= 1.5 && r0 <= 3.5 && ar >= 0.20 && ar <= 0.80 ? 1 : 0}`);
  }
  writeFileSync(`${outDir}/s08_calibration.csv`, calLines.join("\n") + "\n");

  const s9 = new Sim(c, defaultRP(500, 30, 450, 0), false);
  simRun(s9);
  writeFileSync(`${outDir}/s09_baseline_small.csv`, `N,days,total_trans,deaths,inv_ok\n500,30,${s9.total_trans},${s9.daily_cnt[29][ST_D]},${s9.inv_violations === 0 ? 1 : 0}\n`);
}

function writeDailyCSV(s: Sim, outDir: string, fname: string): void {
  const lines = ["day,S,E,Ia,Is,H,ICU,R,D"];
  for (let d = 0; d < s.rp.days; d++) lines.push(`${d},${s.daily_cnt[d].join(",")}`);
  writeFileSync(`${outDir}/${fname}`, lines.join("\n") + "\n");
}

function writeTransCSV(s: Sim, outDir: string, fname: string): void {
  const lines = ["day,src_id,dst_id,layer,variant,gen_int"];
  for (const r of s.trans) lines.push(r.join(","));
  writeFileSync(`${outDir}/${fname}`, lines.join("\n") + "\n");
}

export function writeS10ToS37(sNo: Sim, sMain: Sim, c: Config, outDir: string, betaCal: number): void {
  const dCnt = sMain.rp.days;
  writeDailyCSV(sNo, outDir, "daily_counts_nointerv.csv");
  writeTransCSV(sNo, outDir, "transmissions_log_nointerv.csv");
  writeDailyCSV(sMain, outDir, "daily_counts.csv");
  writeTransCSV(sMain, outDir, "transmissions_log.csv");

  const wLines = (fname: string, hdr: string, rows: string[]) => writeFileSync(`${outDir}/${fname}`, hdr + "\n" + (rows.length ? rows.join("\n") + "\n" : ""));
  wLines("incidence_by_age.csv", "day,age_0_17,age_18_59,age_60p", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_inc_age[d].join(",")}`));
  wLines("hospital_daily.csv", "day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_hosp[d].join(",")}`));
  wLines("deaths_by_age.csv", "day,d_0_17,d_18_59,d_60p", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_deaths_age[d].join(",")}`));
  wLines("tests.csv", "day,tests,positives", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_tests[d][0]},${sMain.daily_tests[d][1]}`));
  wLines("quarantine.csv", "day,quarantined", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_quar[d]}`));
  wLines("vaccination.csv", "day,dose1,dose2", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_vax[d][0]},${sMain.daily_vax[d][1]}`));
  wLines("variant_share.csv", "day,v0_inf,v1_inf", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_var[d][0]},${sMain.daily_var[d][1]}`));
  wLines("generation_interval.csv", "idx,gen_int", sMain.trans.slice(0, 2000).map((r, i) => `${i},${r[5]}`));
  wLines("agents_final.csv", "id,age_grp,status,vax_doses,offspring,ever_inf", sMain.pop.slice(0, 2000).map((p, i) => `${i},${p.age_grp},${p.status},${p.vax_doses},${p.offspring},${p.ever_inf}`));
  wLines("progression_log.csv", "day,pid,from_st,to_st", sMain.progs.map((r) => r.join(",")));

  for (let a = 0; a < 3; a++) wLines(`daily_age_${a}.csv`, "day,inc,deaths", Array.from({ length: dCnt }, (_, d) => `${d},${sMain.daily_inc_age[d][a]},${sMain.daily_deaths_age[d][a]}`));
  ["layer_hh.csv", "layer_sw.csv", "layer_comm.csv"].forEach((nm, l) => writeFileSync(`${outDir}/${nm}`, `layer,infections\n${l},${sMain.layer_inf[l]}\n`));
  wLines("hospital_cap_nointerv.csv", "day,hosp_occ,icu_occ", Array.from({ length: dCnt }, (_, d) => `${d},${sNo.daily_hosp[d][0]},${sNo.daily_hosp[d][1]}`));
  wLines("deaths_nointerv.csv", "day,deaths", Array.from({ length: dCnt }, (_, d) => `${d},${sNo.daily_cnt[d][ST_D]}`));
  const ohBins = [0, 0, 0, 0, 0, 0];
  for (const p of sNo.pop) if (p.ever_inf !== 0) ohBins[Math.min(p.offspring, 5)]++;
  wLines("offspring_hist.csv", "bin,count", ohBins.map((v, b) => `${b},${v}`));
  wLines("rt_nointerv.csv", "day,inc", Array.from({ length: dCnt }, (_, d) => `${d},${sNo.daily_inc_age[d][0] + sNo.daily_inc_age[d][1] + sNo.daily_inc_age[d][2]}`));

  const totViol = sNo.inv_violations + sMain.inv_violations;
  writeFileSync(`${outDir}/s12_invariants.csv`, `violations,ok\n${totViol},${totViol === 0 ? 1 : 0}\n`);
  wLines("s13_rt.csv", "day,rt", Array.from({ length: dCnt }, (_, d) => {
    let num = 0, den = 0;
    for (let k = 0; k < 5; k++) if (d >= k) num += sMain.daily_inc_age[d - k][0] + sMain.daily_inc_age[d - k][1] + sMain.daily_inc_age[d - k][2];
    for (let k = 4; k < 9; k++) if (d >= k) den += sMain.daily_inc_age[d - k][0] + sMain.daily_inc_age[d - k][1] + sMain.daily_inc_age[d - k][2];
    return `${d},${(den > 0 ? num / den : 0.0).toFixed(4)}`;
  }));

  const [chh, csw, ccom] = netDegrees(sNo);
  const r0Ana = analyticalR0(c, betaCal, chh, csw, ccom);
  const sumInc = (s: Sim, d: number) => s.daily_inc_age[d][0] + s.daily_inc_age[d][1] + s.daily_inc_age[d][2];
  const iEarly = (sumInc(sNo, 3) + sumInc(sNo, 4) + 1) / 2.0;
  let iMid = 0;
  for (let d = 10; d < Math.min(16, dCnt); d++) iMid += sumInc(sNo, d);
  iMid /= 6.0;
  const rGrow = iMid > iEarly ? Math.log(iMid / iEarly) / 9.0 : 0.0;
  writeFileSync(`${outDir}/s14_growth.csv`, `growth_r,doubling_days\n${rGrow.toFixed(6)},${(rGrow > 0 ? Math.log(2.0) / rGrow : 0.0).toFixed(4)}\n`);
  const r0Emp = 1.0 + rGrow * 4.5;
  writeFileSync(`${outDir}/s15_r0.csv`, `r0_analytic,r0_empirical,ratio\n${r0Ana.toFixed(4)},${r0Emp.toFixed(4)},${(r0Ana > 0 ? r0Emp / r0Ana : 0.0).toFixed(4)}\n`);
  const zAna = solveFinalSize(r0Ana), herd = r0Ana > 1.0 ? 1.0 - 1.0 / r0Ana : 0.0;
  writeFileSync(`${outDir}/s16_final_size.csv`, `z_analytic,herd_thresh,attack_nointerv,attack_interv\n${zAna.toFixed(4)},${herd.toFixed(4)},${(sNo.total_trans / sNo.n).toFixed(4)},${(sMain.total_trans / sMain.n).toFixed(4)}\n`);
  writeFileSync(`${outDir}/s17_attack.csv`, `category,count,frac\nHH,${sMain.layer_inf[0]},${(sMain.layer_inf[0] / sMain.n).toFixed(4)}\nSW,${sMain.layer_inf[1]},${(sMain.layer_inf[1] / sMain.n).toFixed(4)}\nCOMM,${sMain.layer_inf[2]},${(sMain.layer_inf[2] / sMain.n).toFixed(4)}\nAGE0,${sMain.age_inf[0]},${(sMain.age_inf[0] / Math.max(sMain.age_tot[0], 1)).toFixed(4)}\nAGE1,${sMain.age_inf[1]},${(sMain.age_inf[1] / Math.max(sMain.age_tot[1], 1)).toFixed(4)}\nAGE2,${sMain.age_inf[2]},${(sMain.age_inf[2] / Math.max(sMain.age_tot[2], 1)).toFixed(4)}\n`);

  let m1 = 0, m2 = 0, nInf = 0;
  for (const p of sNo.pop) if (p.ever_inf !== 0) { const o = p.offspring; m1 += o; m2 += o * o; nInf++; }
  if (nInf > 0) { m1 /= nInf; m2 /= nInf; }
  const varO = m2 - m1 * m1;
  const kDisp = varO > m1 && m1 > 0 ? (m1 * m1) / (varO - m1) : 10.0;
  writeFileSync(`${outDir}/s18_offspring.csv`, `n_inf,mean_offspring,var_offspring,k_disp\n${nInf},${m1.toFixed(4)},${varO.toFixed(4)},${kDisp.toFixed(4)}\n`);

  let pkHNo = 0, pkINo = 0, pkHM = 0, pkIM = 0;
  for (let d = 0; d < dCnt; d++) {
    if (sNo.daily_hosp[d][0] > pkHNo) pkHNo = sNo.daily_hosp[d][0];
    if (sNo.daily_hosp[d][1] > pkINo) pkINo = sNo.daily_hosp[d][1];
    if (sMain.daily_hosp[d][0] > pkHM) pkHM = sMain.daily_hosp[d][0];
    if (sMain.daily_hosp[d][1] > pkIM) pkIM = sMain.daily_hosp[d][1];
  }
  writeFileSync(`${outDir}/s19_hospital.csv`, `hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n${sMain.hosp_cap},${sMain.icu_cap},${pkHNo},${pkINo},${sNo.branches[9] + sNo.branches[11]},${pkHM},${pkIM},${sMain.branches[9] + sMain.branches[11]}\n`);
  wLines("s20_ifr.csv", "age_grp,infections,deaths,ifr", [0, 1, 2].map((a) => {
    const dt = sMain.daily_deaths_age.slice(0, dCnt).reduce((acc, r) => acc + r[a], 0);
    return `${a},${sMain.age_inf[a]},${dt},${(sMain.age_inf[a] > 0 ? dt / sMain.age_inf[a] : 0.0).toFixed(6)}`;
  }));
  const tTot = sMain.daily_tests.slice(0, dCnt).reduce((a, r) => a + r[0], 0);
  const tPos = sMain.daily_tests.slice(0, dCnt).reduce((a, r) => a + r[1], 0);
  writeFileSync(`${outDir}/s21_testing.csv`, `tests,positives,quar_traces\n${tTot},${tPos},${sMain.branches[17]}\n`);

  const vn = [0, 0, 0], vi = [0, 0, 0];
  for (const p of sMain.pop) { const vd = Math.min(p.vax_doses, 2); vn[vd]++; if (p.ever_inf !== 0) vi[vd]++; }
  wLines("s22_vaccine.csv", "doses,n_people,n_inf,attack_rate", [0, 1, 2].map((d) => `${d},${vn[d]},${vi[d]},${(vn[d] > 0 ? vi[d] / vn[d] : 0.0).toFixed(6)}`));
  const v0 = sMain.daily_var.slice(0, dCnt).reduce((a, r) => a + r[0], 0);
  const v1 = sMain.daily_var.slice(0, dCnt).reduce((a, r) => a + r[1], 0);
  writeFileSync(`${outDir}/s23_variant.csv`, `v0_total,v1_total,v1_share\n${v0},${v1},${(v0 + v1 > 0 ? v1 / (v0 + v1) : 0.0).toFixed(6)}\n`);
  writeFileSync(`${outDir}/s24_waning.csv`, `waning_events,reinfections\n${sMain.branches[14]},${sMain.total_reinf}\n`);

  const r1 = defaultRP(500, 30, betaCal, 0); r1.contact_scale_bp = 0; const nc1 = new Sim(c, r1, false); simRun(nc1);
  const r2 = defaultRP(500, 30, 0, 0); r2.contact_scale_bp = -1; const nc2 = new Sim(c, r2, false); simRun(nc2);
  const r3 = defaultRP(500, 30, betaCal, 1); r3.vax_perfect = 1; const nc3 = new Sim(c, r3, false); simRun(nc3);
  writeFileSync(`${outDir}/s25_negative_controls.csv`, `control,infections,ok\nzero_contact,${nc1.total_trans},${nc1.total_trans === 0 ? 1 : 0}\nzero_seed,${nc2.total_trans},${nc2.total_trans === 0 ? 1 : 0}\nperfect_vax,${nc3.total_trans},${nc3.total_trans === 0 ? 1 : 0}\n`);

  const nsw = 1000, dsw = 90;
  const s26: string[] = []; for (const cv of [2000, 6000, 9500]) for (const tv of [1000, 5000, 9000]) for (let rp = 0; rp < 3; rp++) { const c2 = { ...c, master_seed: c.master_seed + rp + 1 }; const r = defaultRP(nsw, dsw, betaCal, 1); r.compliance_ovr = cv; r.test_ovr = tv; const ss = new Sim(c2, r, false); simRun(ss); s26.push(`${cv},${tv},${rp},${(ss.total_trans / nsw).toFixed(6)}`); }
  wLines("s26_sweep_comp_test.csv", "comp_bp,test_bp,rep,attack_rate", s26);
  const s27: string[] = []; for (const vd of [10, 45, 80]) for (const vc of [2000, 5500, 9000]) for (let rp = 0; rp < 3; rp++) { const c2 = { ...c, master_seed: c.master_seed + rp + 10 }; const r = defaultRP(nsw, dsw, betaCal, 1); r.vax_day_ovr = vd; r.vax_cov_ovr = vc; const ss = new Sim(c2, r, false); simRun(ss); s27.push(`${vd},${vc},${rp},${(ss.total_trans / nsw).toFixed(6)}`); }
  wLines("s27_sweep_vax.csv", "vax_day,vax_cov_bp,rep,attack_rate", s27);
  const s28: string[] = []; for (const sc of [15, 200]) for (const lk of [15, 45, 200]) for (let rp = 0; rp < 3; rp++) { const c2 = { ...c, master_seed: c.master_seed + rp + 20 }; const r = defaultRP(nsw, dsw, betaCal, 1); r.school_ovr = sc; r.lock_ovr = lk; const ss = new Sim(c2, r, false); simRun(ss); s28.push(`${sc},${lk},${rp},${(ss.total_trans / nsw).toFixed(6)}`); }
  wLines("s28_sweep_school_lock.csv", "school_day,lock_day,rep,attack_rate", s28);
  const s29: string[] = []; for (const vr of [15, 45, 200]) for (const mk of [10, 40, 200]) for (let rp = 0; rp < 3; rp++) { const c2 = { ...c, master_seed: c.master_seed + rp + 30 }; const r = defaultRP(nsw, dsw, betaCal, 1); r.var_day_ovr = vr; r.mask_ovr = mk; const ss = new Sim(c2, r, false); simRun(ss); s29.push(`${vr},${mk},${rp},${(ss.total_trans / nsw).toFixed(6)}`); }
  wLines("s29_sweep_variant_mask.csv", "var_day,mask_day,rep,attack_rate", s29);

  const ensS = new Array<number>(90).fill(0), ensS2 = new Array<number>(90).fill(0);
  for (let rp = 0; rp < 20; rp++) {
    const c2 = { ...c, master_seed: c.master_seed + rp + 100 };
    const ss = new Sim(c2, defaultRP(nsw, dsw, betaCal, 1), false); simRun(ss);
    for (let d = 0; d < 90; d++) { const inc = sumInc(ss, d); ensS[d] += inc; ensS2[d] += inc * inc; }
  }
  wLines("s30_ensemble.csv", "day,mean_inc,se_inc", Array.from({ length: 90 }, (_, d) => {
    const mn = ensS[d] / 20.0;
    return `${d},${mn.toFixed(4)},${Math.sqrt(Math.max((ensS2[d] / 20.0 - mn * mn) / 19.0, 0.0)).toFixed(4)}`;
  }));

  const bn = ["INF_HH", "INF_SW", "INF_COMM", "FAIL_VAX", "FAIL_MASK", "FAIL_QUAR", "PROG_IA", "PROG_IS", "HOSP_ADM", "HOSP_REJ", "ICU_ADM", "ICU_REJ", "DEATH", "RECOVER", "WANING", "TEST_POS", "TEST_NEG", "QUAR_TRACE", "VAX_DOSE", "VAX_OUT", "VAR_IMPORT"];
  wLines("s31_branches.csv", "branch,count", bn.map((nm, i) => `${nm},${sMain.branches[i]}`));
  writeFileSync(`${outDir}/s32_vacuous_check.csv`, `check,ok\nneg_controls,${nc1.total_trans === 0 && nc2.total_trans === 0 && nc3.total_trans === 0 ? 1 : 0}\nensemble_var,${ensS.filter((x) => x > 0).length > 5 ? 1 : 0}\n`);

  const sp1 = new Sim(c, defaultRP(c.parity_N, c.parity_days, betaCal, 1), false); simRun(sp1);
  const sp2 = new Sim(c, defaultRP(c.parity_N, c.parity_days, betaCal, 1), false); simRun(sp2);
  let h1 = FNV_OFFSET, h2 = FNV_OFFSET;
  for (let i = 0; i < sp1.n; i++) {
    h1 = fnv1aI64(h1, sp1.pop[i].status); h2 = fnv1aI64(h2, sp2.pop[i].status);
    h1 = fnv1aI64(h1, sp1.pop[i].offspring); h2 = fnv1aI64(h2, sp2.pop[i].offspring);
  }
  writeFileSync(`${outDir}/s33_replay.csv`, `ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n${hex64(sp1.event_hash)},${hex64(sp2.event_hash)},${hex64(h1)},${hex64(h2)},${sp1.event_hash === sp2.event_hash && h1 === h2 ? 1 : 0}\n`);
  writeFileSync(`${outDir}/s34_parity.csv`, `parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n${c.parity_N},${c.parity_days},${hex64(sp1.event_hash)},${hex64(h1)},${sp1.total_trans},${sp1.daily_cnt[c.parity_days - 1][ST_D]},${sp1.inv_violations}\n`);
  writeFileSync(`${outDir}/s36_sizes.csv`, `N,days,n_trans,n_progs\n${sMain.n},${sMain.rp.days},${sMain.trans.length},${sMain.progs.length}\n`);

  const cntR = [sMain.n, 0, 0, 0, 0, 0, 0, 0];
  let pi = 0, maxDiff = 0;
  for (let d = 0; d < dCnt; d++) {
    while (pi < sMain.progs.length && sMain.progs[pi][0] === d) { cntR[sMain.progs[pi][2]]--; cntR[sMain.progs[pi][3]]++; pi++; }
    for (let k = 0; k < 8; k++) maxDiff = Math.max(maxDiff, Math.abs(cntR[k] - sMain.daily_cnt[d][k]));
  }
  writeFileSync(`${outDir}/s37_replay_counts.csv`, `max_diff,match\n${maxDiff},${maxDiff === 0 ? 1 : 0}\n`);
}
