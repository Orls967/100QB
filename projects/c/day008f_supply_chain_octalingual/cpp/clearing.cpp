#include "common.hpp"

void sim_risk_check(Sim &s, int64_t day, int64_t tick) {
    for (int i = 0; i < NUM_STORES + NUM_WAREHOUSES; ++i) {
        int64_t bo_val = 0;
        for (int k = 0; k < NUM_SKUS; ++k) bo_val += s.acc[i].backorder[k] * SKU_WHOLESALE[k];
        int64_t req_margin = (bo_val * s.cfg.credit_margin_bps) / 10000;
        if (s.acc[i].cash < req_margin && s.acc[i].status == ST_ACTIVE) {
            s.acc[i].status = ST_CREDIT_HOLD; s.branches[BR_CREDIT_HOLD]++;
            if (s.log_enabled && s.f_credit) fprintf(s.f_credit, "%lld,%d,%lld,%lld,1\n", static_cast<long long>(day), i, static_cast<long long>(s.acc[i].cash), static_cast<long long>(req_margin));
        } else if (s.acc[i].status == ST_CREDIT_HOLD && s.acc[i].cash >= req_margin) {
            s.acc[i].status = ST_ACTIVE;
        }
        if (s.acc[i].cash < 0 && s.acc[i].status != ST_LIQUIDATION) {
            s.acc[i].status = ST_LIQUIDATION; s.branches[BR_LIQUIDATE]++;
            int64_t pen = 50000;
            transfer_cash(s, i, ACC_FEE, pen, 2, day, tick);
            if (s.log_enabled && s.f_risk) fprintf(s.f_risk, "%lld,%lld,%d,LIQUIDATE,%lld\n", static_cast<long long>(day), static_cast<long long>(tick), i, static_cast<long long>(pen));
        }
    }
}

void sim_eod_clearing(Sim &s, int64_t day) {
    for (int i = 0; i < NUM_STORES + NUM_WAREHOUSES; ++i) {
        int64_t oh_tot = 0, it_tot = 0, bo_tot = 0;
        for (int k = 0; k < NUM_SKUS; ++k) {
            oh_tot += s.acc[i].on_hand[k]; it_tot += s.acc[i].in_transit[k]; bo_tot += s.acc[i].backorder[k];
            if (day < 6) s.day_inv_snap[day][k] += s.acc[i].on_hand[k];
        }
        int64_t hold = oh_tot * s.cfg.holding_cost_per_day, pen = bo_tot * s.cfg.backorder_penalty;
        if (hold > 0) { transfer_cash(s, i, ACC_FEE, hold, 0, day, s.cfg.ticks_per_day); if (day < 6) s.daily_cost_echelon[day][3] += hold; }
        if (pen > 0)  { transfer_cash(s, i, ACC_FEE, pen, 1, day, s.cfg.ticks_per_day);  if (day < 6) s.daily_cost_echelon[day][3] += pen; }
        if (s.log_enabled && s.f_acc_ser) fprintf(s.f_acc_ser, "%lld,%d,%lld,%lld,%lld,%d\n", static_cast<long long>(day), i, static_cast<long long>(s.acc[i].cash), static_cast<long long>(oh_tot), static_cast<long long>(it_tot), s.acc[i].status);
    }
    sim_risk_check(s, day, s.cfg.ticks_per_day);
}
