import type { Sim } from "./entities.ts";
import { ACC_FEE, NUM_STORES, NUM_WAREHOUSES, SKU_WHOLESALE, transferCash } from "./inventory.ts";

export function simRiskCheck(s: Sim, day: number, tick: number): void {
  for (let i = 0; i < NUM_STORES + NUM_WAREHOUSES; i++) {
    let boVal = 0;
    for (let k = 0; k < 4; k++) boVal += s.acc[i].backorder[k] * SKU_WHOLESALE[k];
    const reqMargin = Math.floor((boVal * s.cfg.credit_margin_bps) / 10000);
    if (s.acc[i].cash < reqMargin && s.acc[i].status === 0) {
      s.acc[i].status = 1; s.branches[13]++;
      s.logLine("credit_calls.csv", `${day},${i},${s.acc[i].cash},${reqMargin},1`);
    } else if (s.acc[i].status === 1 && s.acc[i].cash >= reqMargin) {
      s.acc[i].status = 0;
    }
    if (s.acc[i].cash < 0 && s.acc[i].status !== 2) {
      s.acc[i].status = 2; s.branches[14]++;
      transferCash(s, i, ACC_FEE, 50000, 2, day, tick);
      s.logLine("risk_actions.csv", `${day},${tick},${i},LIQUIDATE,50000`);
    }
  }
}

export function simEodClearing(s: Sim, day: number): void {
  for (let i = 0; i < NUM_STORES + NUM_WAREHOUSES; i++) {
    const ohTot = s.acc[i].onHand.reduce((a, b) => a + b, 0);
    const itTot = s.acc[i].inTransit.reduce((a, b) => a + b, 0);
    const boTot = s.acc[i].backorder.reduce((a, b) => a + b, 0);
    if (day < 6) for (let k = 0; k < 4; k++) s.dayInvSnap[day][k] += s.acc[i].onHand[k];
    const hold = ohTot * s.cfg.holding_cost_per_day, pen = boTot * s.cfg.backorder_penalty;
    if (hold > 0) { transferCash(s, i, ACC_FEE, hold, 0, day, s.cfg.ticks_per_day); if (day < 6) s.dailyCostEchelon[day][3] += hold; }
    if (pen > 0) { transferCash(s, i, ACC_FEE, pen, 1, day, s.cfg.ticks_per_day); if (day < 6) s.dailyCostEchelon[day][3] += pen; }
    s.logLine("accounts_series.csv", `${day},${i},${s.acc[i].cash},${ohTot},${itTot},${s.acc[i].status}`);
  }
  simRiskCheck(s, day, s.cfg.ticks_per_day);
}
