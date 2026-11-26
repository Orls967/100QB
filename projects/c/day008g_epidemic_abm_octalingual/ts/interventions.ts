import { simSetState } from "./disease.ts";
import { BP, ST_D, ST_IS, ST_S, Sim } from "./population.ts";

export function simIntervStep(s: Sim, day: number): void {
  const n = s.n;
  const vd = s.rp.var_day_ovr >= 0 ? s.rp.var_day_ovr : s.cfg.variant_day;
  if (day === vd) {
    let imported = 0;
    for (let k = 0; k < n && imported < 5; k++) {
      const pid = s.rng[3].randInt(0, n - 1);
      if (s.pop[pid].status === ST_S) {
        s.pop[pid].variant = 1;
        s.pop[pid].ever_inf = 1;
        s.pop[pid].inf_day = day;
        simSetState(s, pid, ST_IS, day, 4);
        s.branches[20]++;
        imported++;
      }
    }
  }
  if (s.rp.interv_on === 0) return;

  if (day >= s.cfg.test_day) {
    const tp = s.rp.test_ovr >= 0 ? s.rp.test_ovr : s.cfg.test_prob_bp;
    const cp = s.rp.compliance_ovr;
    for (let i = 0; i < n; i++) {
      const p = s.pop[i];
      if (p.status === ST_IS && p.quar_until < day) {
        const effComp = cp >= 0 ? cp : p.compliance;
        if (s.rng[3].randInt(1, BP) <= Math.floor((tp * effComp) / BP)) {
          s.daily_tests[day][0]++;
          if (s.rng[3].randInt(1, BP) <= 8500) {
            s.daily_tests[day][1]++;
            s.branches[15]++;
            p.quar_until = day + 7;
            const h = p.hh_id;
            for (let idx = s.hh_off[h]; idx < s.hh_off[h + 1]; idx++) {
              const j = s.hh_mem[idx];
              if (j !== i && s.pop[j].quar_until < day && s.rng[3].randInt(1, BP) <= s.cfg.trace_prob_bp) {
                s.pop[j].quar_until = day + 5;
                s.branches[17]++;
              }
            }
          } else {
            s.branches[16]++;
          }
        }
      }
    }
  }

  const vday = s.rp.vax_day_ovr >= 0 ? s.rp.vax_day_ovr : s.cfg.vax_day;
  if (day >= vday) {
    const covCap = s.rp.vax_cov_ovr >= 0 ? Math.floor((n * s.rp.vax_cov_ovr) / BP) : Math.floor((n * 8000) / BP);
    const dosesToday = Math.max(Math.floor((n * s.cfg.vax_daily_bp) / BP), 1);
    let totalVaxed = 0;
    for (let i = 0; i < n; i++) if (s.pop[i].vax_doses > 0) totalVaxed++;
    if (totalVaxed >= covCap) {
      s.branches[19]++;
      return;
    }
    let given = 0;
    for (let grp = 2; grp >= 0 && given < dosesToday; grp--) {
      for (let i = 0; i < n && given < dosesToday; i++) {
        const p = s.pop[i];
        if (p.age_grp !== grp || p.status === ST_D) continue;
        if (p.vax_doses === 0 && totalVaxed < covCap) {
          p.vax_doses = 1;
          p.vax_day = day;
          given++;
          totalVaxed++;
          s.daily_vax[day][0]++;
          s.branches[18]++;
        } else if (p.vax_doses === 1 && day - p.vax_day >= 21) {
          p.vax_doses = 2;
          p.vax_day = day;
          given++;
          s.daily_vax[day][1]++;
          s.branches[18]++;
        }
      }
    }
    if (given === 0) s.branches[19]++;
  }
}
