import { BP, ST_D, ST_E, ST_H, ST_IA, ST_ICU, ST_IS, ST_R, ST_S, Sim } from "./population.ts";
import { fnv1aI64 } from "./rng.ts";

export function simSetState(s: Sim, pid: number, newSt: number, day: number, dur: number): void {
  const oldSt = s.pop[pid].status;
  s.pop[pid].status = newSt;
  s.pop[pid].state_day = day;
  s.pop[pid].dur_days = dur;
  if (s.log_enabled) s.progs.push([day, pid, oldSt, newSt]);
}

export function simDiseaseStep(s: Sim, day: number): void {
  const n = s.n;
  for (let i = 0; i < n; i++) {
    const p = s.pop[i];
    const elapsed = day - p.state_day;
    const st = p.status;
    if (st === ST_E && elapsed >= p.dur_days) {
      const dur = s.rng[2].randInt(s.cfg.inf_days_min, s.cfg.inf_days_max);
      if (s.rng[2].randInt(1, BP) <= s.cfg.asymp_prob_bp) {
        s.branches[6]++;
        simSetState(s, i, ST_IA, day, dur);
      } else {
        s.branches[7]++;
        simSetState(s, i, ST_IS, day, dur);
      }
    } else if ((st === ST_IA || st === ST_IS) && elapsed >= p.dur_days) {
      if (st === ST_IS) {
        let pHosp = p.age_grp === 0 ? 120 : p.age_grp === 1 ? 450 : 1800;
        if (p.vax_doses > 0) pHosp = Math.floor((pHosp * (BP - s.cfg.vax_hosp_bp)) / BP);
        if (s.rng[4].randInt(1, BP) <= pHosp) {
          if (s.hosp_occ < s.hosp_cap) {
            s.hosp_occ++;
            p.hosp_rej = 0;
            s.branches[8]++;
          } else {
            p.hosp_rej = 1;
            s.branches[9]++;
          }
          const dur = s.rng[4].randInt(4, 8);
          simSetState(s, i, ST_H, day, dur);
          continue;
        }
      }
      s.branches[13]++;
      simSetState(s, i, ST_R, day, s.cfg.waning_min_days);
    } else if (st === ST_H && elapsed >= p.dur_days) {
      if (p.hosp_rej === 0 && s.hosp_occ > 0) s.hosp_occ--;
      const pIcu = p.age_grp === 2 ? 3500 : 1800;
      if (s.rng[4].randInt(1, BP) <= pIcu) {
        if (s.icu_occ < s.icu_cap) {
          s.icu_occ++;
          p.hosp_rej = 0;
          s.branches[10]++;
        } else {
          p.hosp_rej = 1;
          s.branches[11]++;
        }
        const dur = s.rng[4].randInt(4, 7);
        simSetState(s, i, ST_ICU, day, dur);
      } else {
        let pDie = p.age_grp === 0 ? 100 : p.age_grp === 1 ? 500 : 2000;
        if (p.hosp_rej !== 0) pDie *= 2;
        if (s.rng[4].randInt(1, BP) <= pDie) {
          s.branches[12]++;
          s.daily_deaths_age[day][p.age_grp]++;
          simSetState(s, i, ST_D, day, 0);
        } else {
          s.branches[13]++;
          simSetState(s, i, ST_R, day, s.cfg.waning_min_days);
        }
      }
    } else if (st === ST_ICU && elapsed >= p.dur_days) {
      if (p.hosp_rej === 0 && s.icu_occ > 0) s.icu_occ--;
      let pDie = p.age_grp === 0 ? 1500 : p.age_grp === 1 ? 2800 : 4800;
      if (p.hosp_rej !== 0) pDie *= 2;
      if (pDie > BP) pDie = BP;
      if (s.rng[4].randInt(1, BP) <= pDie) {
        s.branches[12]++;
        s.daily_deaths_age[day][p.age_grp]++;
        simSetState(s, i, ST_D, day, 0);
      } else {
        s.branches[13]++;
        simSetState(s, i, ST_R, day, s.cfg.waning_min_days);
      }
    } else if (st === ST_R && elapsed >= p.dur_days) {
      if (s.rng[2].randInt(1, BP) <= s.cfg.waning_bp) {
        s.branches[14]++;
        simSetState(s, i, ST_S, day, 0);
      }
    }
  }

  for (let c = 0; c < s.cand_tgt.length; c++) {
    const tgt = s.cand_tgt[c];
    const src = s.cand_src[c];
    const lay = s.cand_lay[c];
    const vari = s.cand_var[c];
    const pt = s.pop[tgt];
    if (pt.status !== ST_S) continue;
    if (pt.ever_inf !== 0) s.total_reinf++;
    pt.ever_inf = 1;
    pt.variant = vari;
    pt.inf_day = day;
    s.pop[src].offspring++;
    s.layer_inf[lay]++;
    s.branches[lay]++;
    const ag = pt.age_grp;
    s.daily_inc_age[day][ag]++;
    s.age_inf[ag]++;
    s.daily_var[day][vari]++;
    s.total_trans++;
    const gi = Math.max(day - s.pop[src].inf_day, 1);
    if (s.log_enabled) s.trans.push([day, src, tgt, lay, vari, gi]);
    s.event_hash = fnv1aI64(s.event_hash, day * 1000000 + tgt);
    simSetState(s, tgt, ST_E, day, s.cfg.latent_days);
  }
}
