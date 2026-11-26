import { simContactStep } from "./contact.ts";
import { simDiseaseStep, simSetState } from "./disease.ts";
import { simIntervStep } from "./interventions.ts";
import { BP, ST_D, ST_IS, ST_S, Sim } from "./population.ts";
import type { Config } from "./population.ts";

export function analyticalR0(c: Config, betaBp: number, cHH: number, cSW: number, cComm: number): number {
  const beta = betaBp / BP;
  const dInf = 0.5 * (c.inf_days_min + c.inf_days_max);
  const pIA = c.asymp_prob_bp / BP;
  const rIA = c.asymp_rel_bp / BP;
  return beta * (cHH + cSW + cComm) * dInf * ((1.0 - pIA) + pIA * rIA);
}

export function solveFinalSize(r0: number): number {
  if (r0 <= 1.0) return 0.0;
  let z = 0.5;
  for (let i = 0; i < 80; i++) z = 1.0 - Math.exp(-r0 * z);
  return z;
}

export function checkInvariants(s: Sim, day: number): boolean {
  const cnt = [0, 0, 0, 0, 0, 0, 0, 0];
  for (let i = 0; i < s.pop.length; i++) cnt[s.pop[i].status]++;
  let sum = 0;
  for (let k = 0; k < 8; k++) {
    sum += cnt[k];
    s.daily_cnt[day][k] = cnt[k];
  }
  let ok = sum === s.n;
  if (day > 0 && s.daily_cnt[day][ST_D] < s.daily_cnt[day - 1][ST_D]) ok = false;
  if (s.hosp_occ > s.hosp_cap || s.icu_occ > s.icu_cap) ok = false;
  if (s.log_enabled && s.trans.length !== s.total_trans) ok = false;
  if (!ok) s.inv_violations++;
  return ok;
}

export function simRun(s: Sim): void {
  const n = s.n;
  let nSeed = s.cfg.init_infected;
  if (s.rp.beta_bp === 0 && s.rp.contact_scale_bp < 0) nSeed = 0;
  if (s.rp.vax_perfect !== 0) {
    for (let i = 0; i < n; i++) {
      s.pop[i].vax_doses = 2;
      s.pop[i].vax_day = 0;
    }
  }
  if (nSeed > 0) {
    let seeded = 0;
    for (let k = 0; k < n && seeded < nSeed; k++) {
      const pid = s.rng[2].randInt(0, n - 1);
      if (s.pop[pid].status === ST_S) {
        s.pop[pid].ever_inf = 1;
        s.pop[pid].inf_day = 0;
        simSetState(s, pid, ST_IS, 0, 4);
        seeded++;
      }
    }
  }
  for (let d = 0; d < s.rp.days; d++) {
    simIntervStep(s, d);
    simContactStep(s, d);
    simDiseaseStep(s, d);
    s.daily_hosp[d][0] = s.hosp_occ;
    s.daily_hosp[d][1] = s.icu_occ;
    s.daily_hosp[d][2] = s.branches[9];
    s.daily_hosp[d][3] = s.branches[11];
    let q = 0;
    for (let i = 0; i < n; i++) if (s.pop[i].quar_until >= d) q++;
    s.daily_quar[d] = q;
    checkInvariants(s, d);
  }
}
