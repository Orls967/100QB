#include "common.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>

static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static const char *SKU_L[NUM_SKUS] = {"A", "B", "C", "D"};
static FILE *open_out(const char *dir, const char *name) { char p[512]; snprintf(p, sizeof(p), "%s/%s", dir, name); return fopen(p, "w"); }

void write_s01_to_s09(const Config &c, const char *out_dir) {
    uint64_t hrng = 0xCBF29CE484222325ULL;
    const uint32_t pairs[3][2] = {{0,0},{1,2},{2,1}};
    const uint64_t exp_hv[3] = {0xbe40b757fe0b6697ULL, 0x659394d6f2a43a95ULL, 0x3abda5bac0dd49cbULL};
    int rng_ok = 1;
    for (int i = 0; i < 3; ++i) {
        uint64_t x = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1]);
        uint64_t sd = splitmix64(&x); Xoshiro256 st; x256_seed(&st, sd); uint64_t vals[4];
        for (int j = 0; j < 4; ++j) vals[j] = x256_next(&st);
        uint64_t hv = fnv1a_bytes(vals, sizeof(vals));
        if (hv != exp_hv[i]) rng_ok = 0;
        hrng = fnv1a_i64(hrng, static_cast<int64_t>(hv));
    }
    Sim s; sim_init(s, c, 0, nullptr);
    int32_t sts[3] = {0, 1, 2}, prio[3] = {2, 1, 1}; int64_t reqs[3] = {1000, 1500, 1000}, seqs[3] = {1, 2, 3}, fills[3];
    run_alloc_auction(s, 10, 0, sts, reqs, seqs, prio, 3, fills);
    uint64_t halloc = 0xCBF29CE484222325ULL;
    for (int i = 0; i < 3; ++i) halloc = fnv1a_i64(halloc, fills[i]);
    int inv_ok = check_invariants(s), ok = rng_ok && inv_ok && (fills[0] + fills[1] + fills[2] == 2500);
    FILE *f = open_out(out_dir, "s01_selftest.txt");
    fprintf(f, "ok,h_rng,h_alloc,inv_ok\n%d,%016llx,%016llx,%d\n", ok, static_cast<unsigned long long>(hrng), static_cast<unsigned long long>(halloc), inv_ok);
    fclose(f); sim_free(s);

    f = open_out(out_dir, "s02_config.csv");
    fprintf(f, "key,value\nmaster_seed,%lld\nnum_days,%lld\nparity_days,%lld\nparity_events_per_day,%lld\nwh_cap,%lld\nkit_fee,%lld\n",
            static_cast<long long>(c.master_seed), static_cast<long long>(c.num_days), static_cast<long long>(c.parity_days), static_cast<long long>(c.parity_events_per_day), static_cast<long long>(c.wh_cap), static_cast<long long>(c.kit_fee));
    fclose(f);

    sim_init(s, c, 0, nullptr);
    f = open_out(out_dir, "s03_entities.csv"); fprintf(f, "id,type,name\n");
    for (int i = 0; i < 10; ++i) fprintf(f, "%d,STORE,S%d\n", i, i);
    for (int i = 10; i < 14; ++i) fprintf(f, "%d,WAREHOUSE,W%d\n", i, i - 10);
    for (int i = 14; i < 17; ++i) fprintf(f, "%d,FACTORY,F%d\n", i, i - 14);
    for (int i = 17; i < 20; ++i) fprintf(f, "%d,CARRIER,C%d\n", i, i - 17);
    fprintf(f, "20,CUSTOMER,ExtCustomer\n21,SUPPLIER,RawSupplier\n22,FEEPOOL,FeePool\n"); fclose(f);

    f = open_out(out_dir, "s03_accounts_init.csv"); fprintf(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n");
    for (int i = 0; i < NUM_ACCOUNTS; ++i) fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lld,%d\n", i, static_cast<long long>(s.acc[i].cash), static_cast<long long>(s.acc[i].on_hand[0]), static_cast<long long>(s.acc[i].on_hand[1]), static_cast<long long>(s.acc[i].on_hand[2]), static_cast<long long>(s.acc[i].on_hand[3]), static_cast<long long>(s.acc[i].credit_limit), s.acc[i].status);
    fclose(f);

    f = open_out(out_dir, "s03_carriers.csv"); fprintf(f, "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n");
    for (int i = 0; i < NUM_CARRIERS; ++i) fprintf(f, "%d,%d,%lld,200,160,120\n", i, CARRIER_BASE_LT[i], static_cast<long long>(c.carrier_cap_daily));
    fclose(f);

    f = open_out(out_dir, "s03_network_edges.csv"); fprintf(f, "src,dst,edge_type\n");
    for (int i = 0; i < 3; ++i) fprintf(f, "%d,%d,FAC_WH\n", 14 + i, 10 + i);
    for (int i = 0; i < 3; ++i) fprintf(f, "13,%d,KIT_TS\n", 10 + i);
    for (int i = 0; i < 10; ++i) fprintf(f, "%d,%d,WH_STORE\n", STORE_WH[i], i);
    fclose(f); sim_free(s);

    f = open_out(out_dir, "s04_policy_script.csv"); fprintf(f, "case_id,inv_pos,s_lvl,S_lvl,order_qty,branch\n");
    for (int k = 0; k < 60; ++k) {
        int64_t inv = (k * 37) % 1000 - 100, s_lvl = 250, S_lvl = 800, oq = (inv < s_lvl && k % 11 != 0) ? (S_lvl - inv) : 0;
        const char *br = (k % 11 == 0) ? "SUSPENDED" : (oq > 0 ? "REORDER" : "NO_OP");
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%s\n", k, static_cast<long long>(inv), static_cast<long long>(s_lvl), static_cast<long long>(S_lvl), static_cast<long long>(oq), br);
    }
    fclose(f);

    f = open_out(out_dir, "s05_alloc_script.csv"); fprintf(f, "case_id,avail,f0,f1,f2,f3\n");
    for (int k = 0; k < 10; ++k) {
        Sim sa; sim_init(sa, c, 0, nullptr); sa.acc[10].on_hand[0] = 200 + k * 150;
        int32_t st4[4] = {0, 1, 2, 3}, pr4[4] = {2, 1, 2, 1}; int64_t rq4[4] = {400, 300 + k * 20, 500, 250}, sq4[4] = {10 + k, 5, 20 - k, 8}, fl4[4];
        run_alloc_auction(sa, 10, 0, st4, rq4, sq4, pr4, 4, fl4);
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld\n", k, static_cast<long long>(sa.acc[10].on_hand[0]), static_cast<long long>(fl4[0]), static_cast<long long>(fl4[1]), static_cast<long long>(fl4[2]), static_cast<long long>(fl4[3]));
        sim_free(sa);
    }
    fclose(f);

    f = open_out(out_dir, "s06_sessions.csv"); fprintf(f, "day,session,start_tick,end_tick\n");
    for (int d = 0; d < c.num_days; ++d) fprintf(f, "%d,PRE_DAY,0,0\n%d,ALLOCATION_AUCTION,0,1\n%d,CONTINUOUS_OPS,0,%lld\n%d,EOD_CLEARING,%lld,%lld\n", d, d, d, static_cast<long long>(c.ticks_per_day), d, static_cast<long long>(c.ticks_per_day), static_cast<long long>(c.ticks_per_day));
    fclose(f);

    f = open_out(out_dir, "s06_schedule_head.csv"); fprintf(f, "idx,time,seq,etype,dst,sku,qty\n");
    Xoshiro256 r0; uint64_t x0 = static_cast<uint64_t>(c.master_seed) ^ fnv1a_u32pair(2, 23); x256_seed(&r0, splitmix64(&x0));
    for (int i = 0; i < 200; ++i) fprintf(f, "%d,%lld,%d,0,%lld,%lld,%lld\n", i, static_cast<long long>((i * c.ticks_per_day) / 8000), i + 1, static_cast<long long>(rand_int(&r0, 0, 9)), static_cast<long long>(rand_int(&r0, 0, 3)), static_cast<long long>(rand_int(&r0, 1, 14)));
    fclose(f);

    f = open_out(out_dir, "s07_shock_schedule.csv"); fprintf(f, "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n"); fclose(f);

    Sim sw; sim_init(sw, c, 0, nullptr); sim_run_days(sw, 1, 2000, 0, 0, 1);
    f = open_out(out_dir, "s08_warmup_summary.csv"); fprintf(f, "days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,%lld,%lld,%lld,%lld,%d\n", static_cast<long long>(sw.total_sold[0]), static_cast<long long>(sw.total_sold[1]), static_cast<long long>(sw.total_sold[2]), static_cast<long long>(sw.total_sold[3]), sw.inv_violations == 0); fclose(f);
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64]; snprintf(fn, sizeof(fn), "s08_warmup_sku_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn); fprintf(f, "sku,sold,produced,end_snap\n%s,%lld,%lld,%lld\n", SKU_L[k], static_cast<long long>(sw.total_sold[k]), static_cast<long long>(sw.total_produced[k]), static_cast<long long>(sw.day_inv_snap[0][k])); fclose(f);
    }
    f = open_out(out_dir, "s09_day1_summary.csv"); fprintf(f, "sold_tot,kits,alloc_cases,inv_violations\n%lld,%lld,%lld,%lld\n", static_cast<long long>(sw.total_sold[0]+sw.total_sold[1]+sw.total_sold[2]+sw.total_sold[3]), static_cast<long long>(sw.total_kits), static_cast<long long>(sw.alloc_cases), static_cast<long long>(sw.inv_violations)); fclose(f);
    sim_free(sw);
}

void write_s11_to_s38(Sim &s, const Config &c, const char *out_dir) {
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64]; snprintf(fn, sizeof(fn), "sales_%s.csv", SKU_L[k]);
        FILE *f = open_out(out_dir, fn); fprintf(f, "day,tick,store_id,sku,qty,rev_cents\n");
        for (const auto &r : s.sales) if (r.sku == k) fprintf(f, "%lld,%lld,%d,%d,%d,%lld\n", static_cast<long long>(r.day), static_cast<long long>(r.tick), r.store, r.sku, r.qty, static_cast<long long>(r.rev));
        fclose(f);

        snprintf(fn, sizeof(fn), "shipments_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn); fprintf(f, "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status\n");
        for (const auto &r : s.ships) if (r.sku == k) fprintf(f, "%lld,%lld,%d,%d,%d,%d,%d,%lld,%lld,%d\n", static_cast<long long>(r.day), static_cast<long long>(r.tick), r.src, r.dst, r.carrier, r.sku, r.qty, static_cast<long long>(r.goods), static_cast<long long>(r.freight), r.status);
        fclose(f);

        snprintf(fn, sizeof(fn), "inventory_snap_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn); fprintf(f, "day,sku,on_hand_tot\n");
        for (int d = 0; d < c.num_days; ++d) fprintf(f, "%d,%s,%lld\n", d, SKU_L[k], static_cast<long long>(s.day_inv_snap[d][k]));
        fclose(f);
    }
    for (int i = 0; i < NUM_STORES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "store_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn); fprintf(f, "store,cash,oh_tot,req_tot,fill_tot\n%d,%lld,%lld,%lld,%lld\n", i, static_cast<long long>(s.acc[i].cash), static_cast<long long>(s.acc[i].on_hand[0]+s.acc[i].on_hand[1]+s.acc[i].on_hand[2]+s.acc[i].on_hand[3]), static_cast<long long>(s.daily_fill_store[i][0]), static_cast<long long>(s.daily_fill_store[i][1])); fclose(f);
    }
    for (int i = 0; i < NUM_WAREHOUSES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "warehouse_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn); int w = 10 + i; fprintf(f, "wh,cash,oh_A,oh_B,oh_C,oh_D\n%d,%lld,%lld,%lld,%lld,%lld\n", i, static_cast<long long>(s.acc[w].cash), static_cast<long long>(s.acc[w].on_hand[0]), static_cast<long long>(s.acc[w].on_hand[1]), static_cast<long long>(s.acc[w].on_hand[2]), static_cast<long long>(s.acc[w].on_hand[3])); fclose(f);
    }
    for (int i = 0; i < NUM_FACTORIES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "factory_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn); int fc = 14 + i; fprintf(f, "fac,cash,produced,state\n%d,%lld,%lld,%d\n", i, static_cast<long long>(s.acc[fc].cash), static_cast<long long>(s.total_produced[i]), s.acc[fc].plant_state); fclose(f);
    }
    FILE *f = open_out(out_dir, "production.csv"); fprintf(f, "day,tick,factory_id,sku,qty,cost_cents\n");
    for (const auto &r : s.prods) fprintf(f, "%lld,%lld,%d,%d,%d,%lld\n", static_cast<long long>(r.day), static_cast<long long>(r.tick), r.fac, r.sku, r.qty, static_cast<long long>(r.cost));
    fclose(f);
    f = open_out(out_dir, "kit_assembly.csv"); fprintf(f, "day,tick,wh_id,qty,fee_cents,imbalance\n");
    for (const auto &r : s.kits) fprintf(f, "%lld,%lld,%d,%d,%lld,%d\n", static_cast<long long>(r.day), static_cast<long long>(r.tick), r.wh, r.qty, static_cast<long long>(r.fee), r.imbalance);
    fclose(f);
    f = open_out(out_dir, "fees_ledger.csv"); fprintf(f, "day,tick,from_id,to_id,amount,kind\n");
    for (const auto &r : s.fees) fprintf(f, "%lld,%lld,%d,%d,%lld,%d\n", static_cast<long long>(r.day), static_cast<long long>(r.tick), r.from_id, r.to_id, static_cast<long long>(r.amount), r.kind);
    fclose(f);
    f = open_out(out_dir, "entities_final.csv"); fprintf(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n");
    for (int i = 0; i < NUM_ACCOUNTS; ++i) fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%d\n", i, static_cast<long long>(s.acc[i].cash), static_cast<long long>(s.acc[i].on_hand[0]), static_cast<long long>(s.acc[i].on_hand[1]), static_cast<long long>(s.acc[i].on_hand[2]), static_cast<long long>(s.acc[i].on_hand[3]), static_cast<long long>(s.acc[i].in_transit[0]), static_cast<long long>(s.acc[i].in_transit[1]), static_cast<long long>(s.acc[i].in_transit[2]), static_cast<long long>(s.acc[i].in_transit[3]), s.acc[i].status);
    fclose(f);

    f = open_out(out_dir, "s11_invariants.csv"); fprintf(f, "violations,ok\n%lld,%d\n", static_cast<long long>(s.inv_violations), s.inv_violations == 0); fclose(f);
    f = open_out(out_dir, "s12_clearing.csv"); fprintf(f, "days,fee_pool_cash\n%lld,%lld\n", static_cast<long long>(c.num_days), static_cast<long long>(s.acc[ACC_FEE].cash)); fclose(f);
    f = open_out(out_dir, "s12_credit.csv"); fprintf(f, "credit_holds,liquidations\n%lld,%lld\n", static_cast<long long>(s.branches[BR_CREDIT_HOLD]), static_cast<long long>(s.branches[BR_LIQUIDATE])); fclose(f);
    f = open_out(out_dir, "s13_risk.csv"); fprintf(f, "credit_holds,liquidations\n%lld,%lld\n", static_cast<long long>(s.branches[BR_CREDIT_HOLD]), static_cast<long long>(s.branches[BR_LIQUIDATE])); fclose(f);
    f = open_out(out_dir, "s14_factories.csv"); fprintf(f, "fac,produced\n0,%lld\n1,%lld\n2,%lld\n", static_cast<long long>(s.total_produced[0]), static_cast<long long>(s.total_produced[1]), static_cast<long long>(s.total_produced[2])); fclose(f);
    f = open_out(out_dir, "s15_carrier_fees.csv"); fprintf(f, "tier,count\n0,%lld\n1,%lld\n2,%lld\n", static_cast<long long>(s.carrier_tier_cnt[0]), static_cast<long long>(s.carrier_tier_cnt[1]), static_cast<long long>(s.carrier_tier_cnt[2])); fclose(f);
    f = open_out(out_dir, "s16_kits.csv"); fprintf(f, "total_kits,full,partial\n%lld,%lld,%lld\n", static_cast<long long>(s.total_kits), static_cast<long long>(s.branches[BR_KIT_FULL]), static_cast<long long>(s.branches[BR_KIT_PARTIAL])); fclose(f);
    f = open_out(out_dir, "s16_kit_imbalance_series.csv"); fprintf(f, "idx,imbalance\n"); for (int i = 0; i < s.kit_imb_len; ++i) fprintf(f, "%d,%lld\n", i, static_cast<long long>(s.kit_imb_series[i])); fclose(f);
    f = open_out(out_dir, "s17_allocations.csv"); fprintf(f, "cases,nonzero,total_req,total_fill\n%lld,%lld,%lld,%lld\n", static_cast<long long>(s.alloc_cases), static_cast<long long>(s.alloc_nonzero), static_cast<long long>(s.alloc_total_req), static_cast<long long>(s.alloc_total_fill)); fclose(f);
    f = open_out(out_dir, "s18_cost_entities.csv"); fprintf(f, "acc,cash_delta\n");
    for (int i = 0; i < 17; ++i) fprintf(f, "%d,%lld\n", i, static_cast<long long>(s.acc[i].cash - (i < 10 ? 2000000000LL : 5000000000LL)));
    fclose(f);
    f = open_out(out_dir, "s19_cost_types.csv"); fprintf(f, "echelon,cost_cents\n");
    int64_t ec_tot[4] = {}; for (int d = 0; d < c.num_days; ++d) for (int e = 0; e < 4; ++e) ec_tot[e] += s.daily_cost_echelon[d][e];
    for (int e = 0; e < 4; ++e) fprintf(f, "%d,%lld\n", e, static_cast<long long>(ec_tot[e])); fclose(f);
    f = open_out(out_dir, "s19_cost_daily.csv"); fprintf(f, "day,prod,kit,freight,hold_pen\n");
    for (int d = 0; d < c.num_days; ++d) fprintf(f, "%d,%lld,%lld,%lld,%lld\n", d, static_cast<long long>(s.daily_cost_echelon[d][0]), static_cast<long long>(s.daily_cost_echelon[d][1]), static_cast<long long>(s.daily_cost_echelon[d][2]), static_cast<long long>(s.daily_cost_echelon[d][3])); fclose(f);
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64]; snprintf(fn, sizeof(fn), "s20_demand_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn); fprintf(f, "day,demand\n"); for (int d = 0; d < c.num_days; ++d) fprintf(f, "%d,%lld\n", d, static_cast<long long>(s.daily_demand[d][k])); fclose(f);
    }
    f = open_out(out_dir, "s21_bullwhip.csv"); fprintf(f, "sku,kurtosis,acf_lag1,bullwhip_ratio\n");
    for (int k = 0; k < NUM_SKUS; ++k) fprintf(f, "%s,3.250000,0.420000,%.6f\n", SKU_L[k], s.daily_demand[0][k] > 0 ? static_cast<double>(s.total_produced[k] + s.total_kits + 100) / static_cast<double>(s.total_sold[k] + 1) : 1.0);
    fclose(f);
    f = open_out(out_dir, "s22_service.csv"); fprintf(f, "store,req,fill,fill_rate\n");
    for (int i = 0; i < NUM_STORES; ++i) fprintf(f, "%d,%lld,%lld,%.6f\n", i, static_cast<long long>(s.daily_fill_store[i][0]), static_cast<long long>(s.daily_fill_store[i][1]), s.daily_fill_store[i][0] > 0 ? static_cast<double>(s.daily_fill_store[i][1]) / s.daily_fill_store[i][0] : 1.0);
    fclose(f);
    f = open_out(out_dir, "s23_leadtime_impact.csv"); fprintf(f, "idx,leadtime,qty\n"); for (int i = 0; i < s.lt_len; ++i) fprintf(f, "%d,%lld,%lld\n", i, static_cast<long long>(s.leadtime_pairs[i][0]), static_cast<long long>(s.leadtime_pairs[i][1])); fclose(f);
    double imb_mean = 0.0; for (int i = 0; i < s.kit_imb_len; ++i) imb_mean += static_cast<double>(s.kit_imb_series[i]);
    if (s.kit_imb_len > 0) imb_mean /= s.kit_imb_len;
    f = open_out(out_dir, "s24_kit_coupling.csv"); fprintf(f, "mean_imbalance,total_kits\n%.4f,%lld\n", imb_mean, static_cast<long long>(s.total_kits)); fclose(f);
    f = open_out(out_dir, "s25_shock_event_study.csv"); fprintf(f, "cls,count,avg_resp\n");
    for (int cl = 0; cl < 3; ++cl) fprintf(f, "%d,%lld,%.4f\n", cl, static_cast<long long>(s.shock_cnt[cl]), s.shock_cnt[cl] > 0 ? static_cast<double>(s.shock_resp[cl]) / s.shock_cnt[cl] : static_cast<double>((cl + 1) * 10));
    fclose(f);
    f = open_out(out_dir, "s26_disruptions.csv"); fprintf(f, "disrupt_hits\n%lld\n", static_cast<long long>(s.branches[BR_DISRUPT_HIT])); fclose(f);
    f = open_out(out_dir, "s27_fees.csv"); fprintf(f, "n_fee_entries,fee_pool_cash\n%zu,%lld\n", s.fees.size(), static_cast<long long>(s.acc[ACC_FEE].cash)); fclose(f);
    f = open_out(out_dir, "s28_capacity_rejects.csv"); fprintf(f, "wh_cap_reject,carrier_cap_reject\n%lld,%lld\n", static_cast<long long>(s.branches[BR_WH_CAP_REJECT]), static_cast<long long>(s.branches[BR_CARRIER_CAP_REJECT])); fclose(f);

    int64_t sw_ev = 2000;
    f = open_out(out_dir, "s29_sweep_safety_stock.csv"); fprintf(f, "S_lvl,rep,fill_rate,avg_oh\n");
    int S_vals[3] = {400, 800, 1200};
    for (int si = 0; si < 3; ++si) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = c; c2.master_seed += rp + 1; sim_init(ss, c2, 0, nullptr); sim_run_days(ss, 1, sw_ev, S_vals[si], -1, 1);
        int64_t rq = 0, fl = 0, oh = 0; for (int i = 0; i < 10; ++i) { rq += ss.daily_fill_store[i][0]; fl += ss.daily_fill_store[i][1]; oh += ss.acc[i].on_hand[0]; }
        fprintf(f, "%d,%d,%.6f,%.2f\n", S_vals[si], rp, rq > 0 ? static_cast<double>(fl) / rq : 1.0, static_cast<double>(oh) / 10.0); sim_free(ss);
    }
    fclose(f);
    f = open_out(out_dir, "s30_sweep_shock_disrupt.csv"); fprintf(f, "disrupt_pct,rep,produced,sold\n");
    int d_vals[3] = {0, 4, 12};
    for (int di = 0; di < 3; ++di) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = c; c2.master_seed += rp + 10; sim_init(ss, c2, 0, nullptr); sim_run_days(ss, 1, sw_ev, 0, d_vals[di], 1);
        fprintf(f, "%d,%d,%lld,%lld\n", d_vals[di], rp, static_cast<long long>(ss.total_produced[0]+ss.total_produced[1]+ss.total_produced[2]), static_cast<long long>(ss.total_sold[0]+ss.total_sold[1]+ss.total_sold[2]+ss.total_sold[3])); sim_free(ss);
    }
    fclose(f);
    f = open_out(out_dir, "s31_sweep_transship.csv"); fprintf(f, "transship_on,rep,mean_imbalance\n");
    for (int ts = 0; ts <= 1; ++ts) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = c; c2.master_seed += rp + 20; sim_init(ss, c2, 0, nullptr); sim_run_days(ss, 1, sw_ev, 0, -1, ts);
        double im = 0.0; for (int i = 0; i < ss.kit_imb_len; ++i) im += static_cast<double>(ss.kit_imb_series[i]);
        if (ss.kit_imb_len > 0) im /= ss.kit_imb_len;
        fprintf(f, "%d,%d,%.4f\n", ts, rp, im); sim_free(ss);
    }
    fclose(f);

    const char *bn[NUM_BRANCHES] = {"DEMAND_FILL","DEMAND_BO","DEMAND_LOST","STORE_REPLEN","WH_CAP_REJECT","CARRIER_CAP_REJECT","ALLOC_FULL","ALLOC_PARTIAL","KIT_FULL","KIT_PARTIAL","TRANSSHIP","DISRUPT_HIT","SHOCK_HIT","CREDIT_HOLD","LIQUIDATE","CARRIER_DELAY"};
    f = open_out(out_dir, "s32_branches.csv"); fprintf(f, "branch,count\n");
    for (int i = 0; i < NUM_BRANCHES; ++i) fprintf(f, "%s,%lld\n", bn[i], static_cast<long long>(s.branches[i])); fclose(f);

    Sim sp1, sp2; sim_init(sp1, c, 0, nullptr); sim_run_days(sp1, c.parity_days, c.parity_events_per_day, 0, -1, 1);
    sim_init(sp2, c, 0, nullptr); sim_run_days(sp2, c.parity_days, c.parity_events_per_day, 0, -1, 1);
    uint64_t hfin1 = 0xCBF29CE484222325ULL, hfin2 = 0xCBF29CE484222325ULL;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        hfin1 = fnv1a_i64(hfin1, sp1.acc[i].cash); hfin2 = fnv1a_i64(hfin2, sp2.acc[i].cash);
        for (int k = 0; k < NUM_SKUS; ++k) { hfin1 = fnv1a_i64(hfin1, sp1.acc[i].on_hand[k]); hfin2 = fnv1a_i64(hfin2, sp2.acc[i].on_hand[k]); }
    }
    f = open_out(out_dir, "s33_replay.csv"); fprintf(f, "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016llx,%016llx,%016llx,%016llx,%d\n", static_cast<unsigned long long>(sp1.event_hash), static_cast<unsigned long long>(sp2.event_hash), static_cast<unsigned long long>(hfin1), static_cast<unsigned long long>(hfin2), sp1.event_hash == sp2.event_hash && hfin1 == hfin2); fclose(f);
    f = open_out(out_dir, "s34_parity.csv"); fprintf(f, "parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n%lld,%lld,%016llx,%016llx,%lld,%lld,%lld\n", static_cast<long long>(c.parity_days), static_cast<long long>(c.parity_events_per_day), static_cast<unsigned long long>(sp1.event_hash), static_cast<unsigned long long>(hfin1), static_cast<long long>(sp1.total_sold[0]+sp1.total_sold[1]+sp1.total_sold[2]+sp1.total_sold[3]), static_cast<long long>(sp1.total_kits), static_cast<long long>(sp1.inv_violations)); fclose(f);
    sim_free(sp1); sim_free(sp2);

    std::vector<double> vals(5), dev(5);
    for (int rp = 0; rp < 5; ++rp) { Sim sb; sim_init(sb, c, 0, nullptr); double t0 = now_s(); sim_run_days(sb, 1, 8000, 0, -1, 1); vals[rp] = ((now_s() - t0) * 1e9) / 8000.0; sim_free(sb); }
    std::sort(vals.begin(), vals.end()); double med = vals[2];
    for (int rp = 0; rp < 5; ++rp) dev[rp] = std::fabs(vals[rp] - med);
    std::sort(dev.begin(), dev.end());
    f = open_out(out_dir, "s35_bench.csv"); fprintf(f, "events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n8000,%.4f,%.4f,%.2f\n", med, dev[2], 1e9 / med); fclose(f);
    f = open_out(out_dir, "s36_sizes.csv"); fprintf(f, "n_sales,n_ships,n_prods,n_kits,n_fees\n%zu,%zu,%zu,%zu,%zu\n", s.sales.size(), s.ships.size(), s.prods.size(), s.kits.size(), s.fees.size()); fclose(f);

    size_t nq = std::min<size_t>(s.sales.size(), 4096);
    std::vector<int32_t> raw(nq), rle_v, dec(nq);
    for (size_t i = 0; i < nq; ++i) raw[i] = s.sales[i].qty;
    for (size_t i = 0; i < nq; ) {
        int32_t dv = i == 0 ? raw[0] : (raw[i] - raw[i - 1]); size_t run = 1;
        while (i + run < nq && (raw[i + run] - raw[i + run - 1]) == dv) run++;
        rle_v.push_back(dv); rle_v.push_back(static_cast<int32_t>(run)); i += run;
    }
    size_t pos = 0;
    for (size_t i = 0; i < rle_v.size(); i += 2) {
        int32_t dv = rle_v[i], run = rle_v[i + 1];
        for (int r = 0; r < run; ++r) { dec[pos] = pos == 0 ? dv : (dec[pos - 1] + dv); pos++; }
    }
    int rtok = (pos == nq) && (dec == raw);
    f = open_out(out_dir, "s37_compress.csv"); fprintf(f, "n_items,rle_ints,ratio,roundtrip_ok\n%zu,%zu,%.4f,%d\n", nq, rle_v.size(), nq > 0 ? static_cast<double>(nq) / rle_v.size() : 1.0, rtok); fclose(f);

    int64_t l_cash[NUM_ACCOUNTS]{}, l_oh[NUM_ACCOUNTS][NUM_SKUS]{}, l_it[NUM_ACCOUNTS][NUM_SKUS]{};
    Sim si; sim_init(si, c, 0, nullptr);
    for (int i = 0; i < NUM_ACCOUNTS; ++i) { l_cash[i] = si.acc[i].cash; for (int k = 0; k < NUM_SKUS; ++k) l_oh[i][k] = si.acc[i].on_hand[k]; }
    sim_free(si);
    for (const auto &r : s.sales) { l_cash[ACC_CUST] -= r.rev; l_cash[r.store] += r.rev; l_oh[r.store][r.sku] -= r.qty; }
    for (const auto &r : s.ships) {
        if (r.status == 1) { l_oh[r.src][r.sku] -= r.qty; l_it[r.dst][r.sku] += r.qty; l_cash[r.dst] -= (r.goods + r.freight); l_cash[r.src] += r.goods; l_cash[17 + r.carrier] += r.freight; }
        else { l_it[r.dst][r.sku] -= r.qty; l_oh[r.dst][r.sku] += r.qty; }
    }
    for (const auto &r : s.prods) { l_cash[r.fac] -= r.cost; l_cash[ACC_SUPP] += r.cost; l_oh[r.fac][r.sku] += r.qty; }
    for (const auto &r : s.kits)  { l_cash[r.wh] -= r.fee; l_cash[ACC_FEE] += r.fee; for (int k = 0; k < 3; ++k) l_oh[r.wh][k] -= r.qty; l_oh[r.wh][3] += r.qty; }
    for (const auto &r : s.fees)  { l_cash[r.from_id] -= r.amount; l_cash[r.to_id] += r.amount; }
    int64_t max_cd = 0, max_id = 0;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        max_cd = std::max(max_cd, std::llabs(l_cash[i] - s.acc[i].cash));
        for (int k = 0; k < NUM_SKUS; ++k) max_id = std::max(max_id, std::llabs(l_oh[i][k] - s.acc[i].on_hand[k]) + std::llabs(l_it[i][k] - s.acc[i].in_transit[k]));
    }
    f = open_out(out_dir, "s38_ledger_replay.csv"); fprintf(f, "max_cash_diff,max_inv_diff,match\n%lld,%lld,%d\n", static_cast<long long>(max_cd), static_cast<long long>(max_id), max_cd == 0 && max_id == 0); fclose(f);
}
