#include "common.hpp"

void run_alloc_auction(Sim &s, int32_t wh, int32_t sku, const int32_t *stores, const int64_t *reqs, const int64_t *seqs, const int32_t *prio, int n, int64_t *fills) {
    int64_t avail = s.acc[wh].on_hand[sku], tot_req = 0;
    for (int i = 0; i < n; ++i) { fills[i] = 0; tot_req += reqs[i]; }
    if (tot_req <= 0) return;
    s.alloc_cases++; s.alloc_total_req += tot_req;
    if (avail <= 0) return;
    if (avail >= tot_req) {
        for (int i = 0; i < n; ++i) fills[i] = reqs[i];
        s.alloc_total_fill += tot_req; s.alloc_nonzero++; s.branches[BR_ALLOC_FULL]++;
        return;
    }
    int idx[NUM_STORES];
    for (int i = 0; i < n; ++i) idx[i] = i;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int a = idx[i], b = idx[j];
            bool sw = (prio[b] > prio[a]) || (prio[b] == prio[a] && seqs[b] < seqs[a]) || (prio[b] == prio[a] && seqs[b] == seqs[a] && stores[b] < stores[a]);
            if (sw) { idx[i] = b; idx[j] = a; }
        }
    }
    int64_t rem = avail;
    bool progress = true;
    while (rem > 0 && progress) {
        progress = false;
        int64_t max_short = -1;
        for (int ii = 0; ii < n; ++ii) {
            int64_t sh = reqs[idx[ii]] - fills[idx[ii]];
            if (sh > max_short) max_short = sh;
        }
        if (max_short <= 0) break;
        for (int ii = 0; ii < n && rem > 0; ++ii) {
            int i = idx[ii];
            if (reqs[i] - fills[i] == max_short) { fills[i]++; rem--; progress = true; }
        }
    }
    int64_t filled = avail - rem;
    s.alloc_total_fill += filled;
    if (filled > 0) s.alloc_nonzero++;
    s.branches[BR_ALLOC_PARTIAL]++;
}
