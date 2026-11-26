import type { Sim } from "./entities.ts";

export function runAllocAuction(s: Sim, wh: number, sku: number, stores: number[], reqs: number[], seqs: number[], prio: number[]): number[] {
  const n = stores.length, fills = new Array(n).fill(0);
  const totReq = reqs.reduce((a, b) => a + b, 0);
  if (totReq <= 0) return fills;
  const avail = s.acc[wh].onHand[sku];
  s.allocCases++; s.allocTotalReq += totReq;
  if (avail <= 0) return fills;
  if (avail >= totReq) {
    s.allocTotalFill += totReq; s.allocNonzero++; s.branches[6]++;
    return [...reqs];
  }
  const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) =>
    prio[a] !== prio[b] ? prio[b] - prio[a] : (seqs[a] !== seqs[b] ? seqs[a] - seqs[b] : stores[a] - stores[b])
  );
  let rem = avail;
  while (rem > 0) {
    let maxShort = -1;
    for (const i of idx) if (reqs[i] - fills[i] > maxShort) maxShort = reqs[i] - fills[i];
    if (maxShort <= 0) break;
    let prog = false;
    for (const i of idx) {
      if (rem <= 0) break;
      if (reqs[i] - fills[i] === maxShort) { fills[i]++; rem--; prog = true; }
    }
    if (!prog) break;
  }
  const filled = avail - rem;
  s.allocTotalFill += filled;
  if (filled > 0) s.allocNonzero++;
  s.branches[7]++;
  return fills;
}
