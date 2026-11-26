import type { Sim } from "./entities.ts";

export const NUM_SKUS = 4;
export const NUM_STORES = 10;
export const NUM_WAREHOUSES = 4;
export const NUM_CARRIERS = 3;
export const NUM_ACCOUNTS = 23;
export const ACC_CUST = 20;
export const ACC_SUPP = 21;
export const ACC_FEE = 22;

export const SKU_RETAIL = [10000, 10200, 9800, 31500];
export const SKU_WHOLESALE = [8000, 8200, 7800, 26000];
export const SKU_PROD = [6000, 6200, 5800, 0];
export const STORE_WH = [10, 10, 10, 11, 11, 11, 12, 12, 12, 12];
export const CARRIER_BASE_LT = [150, 250, 200];

export const carrierRate = (v: number): number => v < 500 ? 200 : (v < 2000 ? 160 : 120);

export function transferCash(s: Sim, fromId: number, toId: number, amt: number, kind: number, day: number, tick: number): void {
  if (amt <= 0 || fromId === toId) return;
  s.acc[fromId].cash -= amt;
  s.acc[toId].cash += amt;
  if (kind >= 0) s.fees.push([day, tick, fromId, toId, amt, kind]);
}

export function checkInvariants(s: Sim): number {
  let totCash = 0;
  for (const a of s.acc) totCash += a.cash;
  if (totCash !== s.initCash) { s.invViolations++; return 0; }
  for (let k = 0; k < 4; k++) {
    let oh = 0, it = 0;
    for (const a of s.acc) {
      if (a.onHand[k] < 0 || a.inTransit[k] < 0 || a.backorder[k] < 0) { s.invViolations++; return 0; }
      oh += a.onHand[k]; it += a.inTransit[k];
    }
    const lhs = oh + it + s.totalSold[k] + (k < 3 ? s.totalKits : 0);
    const rhs = s.initStock[k] + (k < 3 ? s.totalProduced[k] : s.totalKits);
    if (lhs !== rhs) { s.invViolations++; return 0; }
  }
  return 1;
}
