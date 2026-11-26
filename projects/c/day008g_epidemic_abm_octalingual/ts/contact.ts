import { BP, ST_IA, ST_IS, ST_S, Sim } from "./population.ts";
import type { Person } from "./population.ts";

function vaxEffBp(s: Sim, p: Person, day: number): number {
  if (p.vax_doses <= 0 || p.vax_day < 0) return 0;
  if (s.rp.vax_perfect !== 0) return BP;
  const base = p.vax_doses === 1 ? s.cfg.vax_eff1_bp : s.cfg.vax_eff2_bp;
  const elapsed = day - p.vax_day;
  return elapsed <= 30 ? base : Math.max(base - (elapsed - 30) * 35, 2000);
}

function tryInfect(s: Sim, src: number, tgt: number, lay: number, day: number, betaEff: number): void {
  const pt = s.pop[tgt];
  if (pt.status !== ST_S || s.seen_day[tgt] === day) return;
  if (pt.quar_until >= day) {
    s.branches[5]++;
    return;
  }
  let pBp = betaEff;
  if (s.rp.interv_on !== 0) {
    const mDay = s.rp.mask_ovr >= 0 ? s.rp.mask_ovr : s.cfg.mask_day;
    if (day >= mDay) {
      const red = Math.floor((s.cfg.mask_eff_bp * s.pop[src].compliance) / BP);
      const nextBp = Math.floor((pBp * (BP - red)) / BP);
      if (nextBp < pBp) s.branches[4]++;
      pBp = nextBp;
    }
  }
  const ve = vaxEffBp(s, pt, day);
  if (ve > 0) {
    pBp = Math.floor((pBp * (BP - ve)) / BP);
    if (pBp <= 0) {
      s.branches[3]++;
      return;
    }
  }
  if (s.rng[1].randInt(1, BP) <= pBp) {
    s.seen_day[tgt] = day;
    s.cand_tgt.push(tgt);
    s.cand_src.push(src);
    s.cand_lay.push(lay);
    s.cand_var.push(s.pop[src].variant);
  } else if (ve > 0) {
    s.branches[3]++;
  }
}

export function simContactStep(s: Sim, day: number): void {
  s.cand_tgt.length = 0;
  s.cand_src.length = 0;
  s.cand_lay.length = 0;
  s.cand_var.length = 0;
  if (s.rp.contact_scale_bp === 0) return;

  const n = s.n;
  const scDay = s.rp.school_ovr >= 0 ? s.rp.school_ovr : s.cfg.school_close_day;
  const lkDay = s.rp.lock_ovr >= 0 ? s.rp.lock_ovr : s.cfg.lockdown_day;
  const schoolClosed = s.rp.interv_on !== 0 && day >= scDay;
  const locked = s.rp.interv_on !== 0 && day >= lkDay;

  for (let i = 0; i < n; i++) {
    const pi = s.pop[i];
    const st = pi.status;
    if (st !== ST_IA && st !== ST_IS) continue;

    let baseBp = s.rp.beta_bp;
    if (st === ST_IA) baseBp = Math.floor((baseBp * s.cfg.asymp_rel_bp) / BP);
    if (pi.variant === 1) baseBp = Math.floor((baseBp * s.cfg.variant_mult_bp) / BP);
    if (baseBp > BP) baseBp = BP;

    const h = pi.hh_id;
    for (let idx = s.hh_off[h]; idx < s.hh_off[h + 1]; idx++) {
      const j = s.hh_mem[idx];
      if (j !== i) tryInfect(s, i, j, 0, day, baseBp);
    }
    if (pi.quar_until >= day) {
      s.branches[5]++;
      continue;
    }

    if (pi.sw_id >= 0 && !(pi.age_grp === 0 && schoolClosed)) {
      const sw = pi.sw_id;
      const lo = s.sw_off[sw];
      const hi = s.sw_off[sw + 1];
      const sz = hi - lo;
      if (sz > 1) {
        let nSwC = pi.age_grp === 0 ? 6 : 5;
        if (locked) nSwC = Math.floor((nSwC + 1) / 2);
        for (let c = 0; c < nSwC; c++) {
          const j = s.sw_mem[lo + s.rng[1].randInt(0, sz - 1)];
          if (j !== i) tryInfect(s, i, j, 1, day, baseBp);
        }
      }
    }

    let nComm = pi.super_flag !== 0 ? 12 : 3;
    if (locked) nComm = Math.floor((nComm + 1) / 3);
    for (let c = 0; c < nComm; c++) {
      const j = s.rng[1].randInt(0, n - 1);
      if (j !== i) tryInfect(s, i, j, 2, day, baseBp);
    }
  }
}
