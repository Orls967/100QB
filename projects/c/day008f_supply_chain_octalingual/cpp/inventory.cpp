#include "common.hpp"

const int64_t SKU_RETAIL[NUM_SKUS]    = {10000, 10200, 9800, 31500};
const int64_t SKU_WHOLESALE[NUM_SKUS] = {8000, 8200, 7800, 26000};
const int64_t SKU_PROD[NUM_SKUS]      = {6000, 6200, 5800, 0};
const int32_t STORE_WH[NUM_STORES]    = {10, 10, 10, 11, 11, 11, 12, 12, 12, 12};
const int32_t CARRIER_BASE_LT[NUM_CARRIERS] = {120, 200, 160};

int64_t carrier_rate(int64_t day_vol) { return day_vol < 500 ? 200 : (day_vol < 2000 ? 160 : 120); }

void transfer_cash(Sim &s, int32_t from_id, int32_t to_id, int64_t amt, int32_t kind, int64_t day, int64_t tick) {
    if (amt <= 0 || from_id == to_id) return;
    s.acc[from_id].cash -= amt;
    s.acc[to_id].cash += amt;
    if (kind >= 0) s.fees.push_back({day, tick, from_id, to_id, kind, amt});
}

int check_invariants(Sim &s) {
    int64_t tot_cash = 0;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) tot_cash += s.acc[i].cash;
    if (tot_cash != s.init_cash) { s.inv_violations++; return 0; }
    for (int k = 0; k < NUM_SKUS; ++k) {
        int64_t oh = 0, it = 0;
        for (int i = 0; i < NUM_ACCOUNTS; ++i) {
            if (s.acc[i].on_hand[k] < 0 || s.acc[i].in_transit[k] < 0 || s.acc[i].backorder[k] < 0) {
                s.inv_violations++; return 0;
            }
            oh += s.acc[i].on_hand[k];
            it += s.acc[i].in_transit[k];
        }
        int64_t lhs = oh + it + s.total_sold[k] + (k < 3 ? s.total_kits : 0);
        int64_t rhs = s.init_stock[k] + (k < 3 ? s.total_produced[k] : s.total_kits);
        if (lhs != rhs) { s.inv_violations++; return 0; }
    }
    return 1;
}
