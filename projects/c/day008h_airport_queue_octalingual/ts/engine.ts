import { N_STN, SLOTS_PER_DAY, STN_BELT } from "./config.ts";
import { Sim, simStepSlot } from "./network.ts";

export function simRun(s: Sim): void {
  const totalSlots = s.rp.days * (s.rp.timescaleDiv2 !== 0 ? SLOTS_PER_DAY * 2 : SLOTS_PER_DAY);
  for (let t = 0; t < totalSlots; t++) {
    simStepSlot(s, t);
  }
  let remPax = 0;
  for (const ag of s.pool) {
    if (ag.active !== 0 && ag.isBag === 0) {
      remPax++;
      s.sysSumSojourn += totalSlots - ag.arrSlot;
    }
  }
  let stnRem = 0;
  for (let st = 0; st < N_STN; st++) {
    if (st !== STN_BELT) stnRem += s.stn[st].len + s.stn[st].busyCnt;
  }
  if (s.totalArr !== s.totalDone + s.totalLost + remPax || stnRem !== remPax) {
    s.invViolations++;
  }
  if (s.sysSumL !== s.sysSumSojourn) {
    s.invViolations++;
  }
}
