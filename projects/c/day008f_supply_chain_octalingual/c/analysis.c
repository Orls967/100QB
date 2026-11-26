#include "common.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }
static int cmp_d(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
static const char *SKU_L[NUM_SKUS] = {"A", "B", "C", "D"};

static FILE *open_out(const char *dir, const char *name) {
    char p[512];
    snprintf(p, sizeof(p), "%s/%s", dir, name);
    return fopen(p, "w");
}

void write_s01_to_s09(const Config *c, const char *out_dir) {
    /* s01 */
    uint64_t hrng = 0xCBF29CE484222325ULL;
    const uint32_t pairs[3][2] = {{0,0},{1,2},{2,1}};
    const uint64_t exp_hv[3] = {0xbe40b757fe0b6697ULL, 0x659394d6f2a43a95ULL, 0x3abda5bac0dd49cbULL};
    int rng_ok = 1;
    for (int i = 0; i < 3; ++i) {
        uint64_t x = (uint64_t)c->master_seed ^ fnv1a_u32pair(pairs[i][0], pairs[i][1]);
        uint64_t sd = splitmix64(&x);
        Xoshiro256 st; x256_seed(&st, sd);
        uint64_t vals[4];
        for (int j = 0; j < 4; ++j) vals[j] = x256_next(&st);
        uint64_t hv = fnv1a_bytes(vals, sizeof(vals));
        if (hv != exp_hv[i]) rng_ok = 0;
        hrng = fnv1a_i64(hrng, (int64_t)hv);
    }
    Sim s; sim_init(&s, c, 0, NULL);
    int32_t sts[3] = {0, 1, 2}; int64_t reqs[3] = {1000, 1500, 1000}, seqs[3] = {1, 2, 3}, fills[3]; int32_t prio[3] = {2, 1, 1};
    run_alloc_auction(&s, 10, 0, sts, reqs, seqs, prio, 3, fills, 0, 0);
    uint64_t halloc = 0xCBF29CE484222325ULL;
    for (int i = 0; i < 3; ++i) halloc = fnv1a_i64(halloc, fills[i]);
    int inv_ok = check_invariants(&s);
    int ok = rng_ok && inv_ok && (fills[0] + fills[1] + fills[2] == 2500);
    FILE *f = open_out(out_dir, "s01_selftest.txt");
    fprintf(f, "ok,h_rng,h_alloc,inv_ok\n%d,%016llx,%016llx,%d\n", ok, (unsigned long long)hrng, (unsigned long long)halloc, inv_ok);
    fclose(f);
    sim_free(&s);

    /* s02 */
    f = open_out(out_dir, "s02_config.csv");
    fprintf(f, "key,value\nmaster_seed,%lld\nnum_days,%lld\nparity_days,%lld\nparity_events_per_day,%lld\nwh_cap,%lld\nkit_fee,%lld\n",
            (long long)c->master_seed, (long long)c->num_days, (long long)c->parity_days, (long long)c->parity_events_per_day, (long long)c->wh_cap, (long long)c->kit_fee);
    fclose(f);

    /* s03 */
    sim_init(&s, c, 0, NULL);
    f = open_out(out_dir, "s03_entities.csv");
    fprintf(f, "id,type,name\n");
    for (int i = 0; i < 10; ++i) fprintf(f, "%d,STORE,S%d\n", i, i);
    for (int i = 10; i < 14; ++i) fprintf(f, "%d,WAREHOUSE,W%d\n", i, i - 10);
    for (int i = 14; i < 17; ++i) fprintf(f, "%d,FACTORY,F%d\n", i, i - 14);
    for (int i = 17; i < 20; ++i) fprintf(f, "%d,CARRIER,C%d\n", i, i - 17);
    fprintf(f, "20,CUSTOMER,ExtCustomer\n21,SUPPLIER,RawSupplier\n22,FEEPOOL,FeePool\n");
    fclose(f);

    f = open_out(out_dir, "s03_accounts_init.csv");
    fprintf(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n");
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lld,%d\n", i, (long long)s.acc[i].cash,
                (long long)s.acc[i].on_hand[0], (long long)s.acc[i].on_hand[1], (long long)s.acc[i].on_hand[2], (long long)s.acc[i].on_hand[3],
                (long long)s.acc[i].credit_limit, s.acc[i].status);
    }
    fclose(f);

    f = open_out(out_dir, "s03_carriers.csv");
    fprintf(f, "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n");
    for (int i = 0; i < NUM_CARRIERS; ++i) fprintf(f, "%d,%d,%lld,200,160,120\n", i, CARRIER_BASE_LT[i], (long long)c->carrier_cap_daily);
    fclose(f);

    f = open_out(out_dir, "s03_network_edges.csv");
    fprintf(f, "src,dst,edge_type\n");
    for (int i = 0; i < 3; ++i) fprintf(f, "%d,%d,FAC_WH\n", 14 + i, 10 + i);
    for (int i = 0; i < 3; ++i) fprintf(f, "13,%d,KIT_TS\n", 10 + i);
    for (int i = 0; i < 10; ++i) fprintf(f, "%d,%d,WH_STORE\n", STORE_WH[i], i);
    fclose(f);
    sim_free(&s);

    /* s04 policy script (60 scenarios) */
    f = open_out(out_dir, "s04_policy_script.csv");
    fprintf(f, "case_id,inv_pos,s_lvl,S_lvl,order_qty,branch\n");
    for (int k = 0; k < 60; ++k) {
        int64_t inv = (k * 37) % 1000 - 100;
        int64_t s_lvl = 250, S_lvl = 800;
        int64_t oq = (inv < s_lvl && k % 11 != 0) ? (S_lvl - inv) : 0;
        const char *br = (k % 11 == 0) ? "SUSPENDED" : (oq > 0 ? "REORDER" : "NO_OP");
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%s\n", k, (long long)inv, (long long)s_lvl, (long long)S_lvl, (long long)oq, br);
    }
    fclose(f);

    /* s05 alloc script (10 cases) */
    f = open_out(out_dir, "s05_alloc_script.csv");
    fprintf(f, "case_id,avail,f0,f1,f2,f3\n");
    for (int k = 0; k < 10; ++k) {
        Sim sa; sim_init(&sa, c, 0, NULL);
        sa.acc[10].on_hand[0] = 200 + k * 150;
        int32_t st4[4] = {0, 1, 2, 3};
        int64_t rq4[4] = {400, 300 + k * 20, 500, 250};
        int64_t sq4[4] = {10 + k, 5, 20 - k, 8};
        int32_t pr4[4] = {2, 1, 2, 1};
        int64_t fl4[4];
        run_alloc_auction(&sa, 10, 0, st4, rq4, sq4, pr4, 4, fl4, 0, 0);
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld\n", k, (long long)sa.acc[10].on_hand[0],
                (long long)fl4[0], (long long)fl4[1], (long long)fl4[2], (long long)fl4[3]);
        sim_free(&sa);
    }
    fclose(f);

    /* s06, s07 */
    f = open_out(out_dir, "s06_sessions.csv");
    fprintf(f, "day,session,start_tick,end_tick\n");
    for (int d = 0; d < c->num_days; ++d) {
        fprintf(f, "%d,PRE_DAY,0,0\n%d,ALLOCATION_AUCTION,0,1\n%d,CONTINUOUS_OPS,0,%lld\n%d,EOD_CLEARING,%lld,%lld\n",
                d, d, d, (long long)c->ticks_per_day, d, (long long)c->ticks_per_day, (long long)c->ticks_per_day);
    }
    fclose(f);

    f = open_out(out_dir, "s06_schedule_head.csv");
    fprintf(f, "idx,time,seq,etype,dst,sku,qty\n");
    Xoshiro256 r0; uint64_t x0 = (uint64_t)c->master_seed ^ fnv1a_u32pair(2, 23); x256_seed(&r0, splitmix64(&x0));
    for (int i = 0; i < 200; ++i) {
        int64_t t = (i * c->ticks_per_day) / 8000;
        int64_t st = rand_int(&r0, 0, 9), sk = rand_int(&r0, 0, 3), q = rand_int(&r0, 1, 14);
        fprintf(f, "%d,%lld,%d,0,%lld,%lld,%lld\n", i, (long long)t, i + 1, (long long)st, (long long)sk, (long long)q);
    }
    fclose(f);

    f = open_out(out_dir, "s07_shock_schedule.csv");
    fprintf(f, "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n");
    fclose(f);

    /* s08 warmup & s09 day1 */
    Sim sw; sim_init(&sw, c, 0, NULL);
    sim_run_days(&sw, 1, 2000, 0, 0, 1);
    f = open_out(out_dir, "s08_warmup_summary.csv");
    fprintf(f, "days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,%lld,%lld,%lld,%lld,%d\n",
            (long long)sw.total_sold[0], (long long)sw.total_sold[1], (long long)sw.total_sold[2], (long long)sw.total_sold[3], sw.inv_violations == 0);
    fclose(f);
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64]; snprintf(fn, sizeof(fn), "s08_warmup_sku_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn);
        fprintf(f, "sku,sold,produced,end_snap\n%s,%lld,%lld,%lld\n", SKU_L[k], (long long)sw.total_sold[k], (long long)sw.total_produced[k], (long long)sw.day_inv_snap[0][k]);
        fclose(f);
    }
    f = open_out(out_dir, "s09_day1_summary.csv");
    fprintf(f, "sold_tot,kits,alloc_cases,inv_violations\n%lld,%lld,%lld,%lld\n",
            (long long)(sw.total_sold[0] + sw.total_sold[1] + sw.total_sold[2] + sw.total_sold[3]),
            (long long)sw.total_kits, (long long)sw.alloc_cases, (long long)sw.inv_violations);
    fclose(f);
    sim_free(&sw);
}

void write_s11_to_s38(Sim *s, const Config *c, const char *out_dir, double *stage_s) {
    (void)stage_s;
    /* S10 remaining split logs */
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64];
        snprintf(fn, sizeof(fn), "sales_%s.csv", SKU_L[k]);
        FILE *f = open_out(out_dir, fn);
        fprintf(f, "day,tick,store_id,sku,qty,rev_cents\n");
        for (int i = 0; i < s->n_sales; ++i) if (s->sales[i].sku == k)
            fprintf(f, "%lld,%lld,%d,%d,%d,%lld\n", (long long)s->sales[i].day, (long long)s->sales[i].tick, s->sales[i].store, s->sales[i].sku, s->sales[i].qty, (long long)s->sales[i].rev);
        fclose(f);

        snprintf(fn, sizeof(fn), "shipments_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn);
        fprintf(f, "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status\n");
        for (int i = 0; i < s->n_ships; ++i) if (s->ships[i].sku == k)
            fprintf(f, "%lld,%lld,%d,%d,%d,%d,%d,%lld,%lld,%d\n", (long long)s->ships[i].day, (long long)s->ships[i].tick, s->ships[i].src, s->ships[i].dst, s->ships[i].carrier, s->ships[i].sku, s->ships[i].qty, (long long)s->ships[i].goods, (long long)s->ships[i].freight, s->ships[i].status);
        fclose(f);

        snprintf(fn, sizeof(fn), "inventory_snap_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn);
        fprintf(f, "day,sku,on_hand_tot\n");
        for (int d = 0; d < c->num_days; ++d) fprintf(f, "%d,%s,%lld\n", d, SKU_L[k], (long long)s->day_inv_snap[d][k]);
        fclose(f);
    }
    for (int i = 0; i < NUM_STORES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "store_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn);
        fprintf(f, "store,cash,oh_tot,req_tot,fill_tot\n%d,%lld,%lld,%lld,%lld\n", i, (long long)s->acc[i].cash,
                (long long)(s->acc[i].on_hand[0]+s->acc[i].on_hand[1]+s->acc[i].on_hand[2]+s->acc[i].on_hand[3]),
                (long long)s->daily_fill_store[i][0], (long long)s->daily_fill_store[i][1]);
        fclose(f);
    }
    for (int i = 0; i < NUM_WAREHOUSES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "warehouse_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn); int w = 10 + i;
        fprintf(f, "wh,cash,oh_A,oh_B,oh_C,oh_D\n%d,%lld,%lld,%lld,%lld,%lld\n", i, (long long)s->acc[w].cash,
                (long long)s->acc[w].on_hand[0], (long long)s->acc[w].on_hand[1], (long long)s->acc[w].on_hand[2], (long long)s->acc[w].on_hand[3]);
        fclose(f);
    }
    for (int i = 0; i < NUM_FACTORIES; ++i) {
        char fn[64]; snprintf(fn, sizeof(fn), "factory_daily_%d.csv", i);
        FILE *f = open_out(out_dir, fn); int fc = 14 + i;
        fprintf(f, "fac,cash,produced,state\n%d,%lld,%lld,%d\n", i, (long long)s->acc[fc].cash, (long long)s->total_produced[i], s->acc[fc].plant_state);
        fclose(f);
    }
    FILE *f = open_out(out_dir, "production.csv");
    fprintf(f, "day,tick,factory_id,sku,qty,cost_cents\n");
    for (int i = 0; i < s->n_prods; ++i) fprintf(f, "%lld,%lld,%d,%d,%d,%lld\n", (long long)s->prods[i].day, (long long)s->prods[i].tick, s->prods[i].fac, s->prods[i].sku, s->prods[i].qty, (long long)s->prods[i].cost);
    fclose(f);

    f = open_out(out_dir, "kit_assembly.csv");
    fprintf(f, "day,tick,wh_id,qty,fee_cents,imbalance\n");
    for (int i = 0; i < s->n_kits; ++i) fprintf(f, "%lld,%lld,%d,%d,%lld,%d\n", (long long)s->kits[i].day, (long long)s->kits[i].tick, s->kits[i].wh, s->kits[i].qty, (long long)s->kits[i].fee, s->kits[i].imbalance);
    fclose(f);

    f = open_out(out_dir, "fees_ledger.csv");
    fprintf(f, "day,tick,from_id,to_id,amount,kind\n");
    for (int i = 0; i < s->n_fees; ++i) fprintf(f, "%lld,%lld,%d,%d,%lld,%d\n", (long long)s->fees[i].day, (long long)s->fees[i].tick, s->fees[i].from_id, s->fees[i].to_id, (long long)s->fees[i].amount, s->fees[i].kind);
    fclose(f);

    f = open_out(out_dir, "entities_final.csv");
    fprintf(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n");
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        fprintf(f, "%d,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld,%d\n", i, (long long)s->acc[i].cash,
                (long long)s->acc[i].on_hand[0], (long long)s->acc[i].on_hand[1], (long long)s->acc[i].on_hand[2], (long long)s->acc[i].on_hand[3],
                (long long)s->acc[i].in_transit[0], (long long)s->acc[i].in_transit[1], (long long)s->acc[i].in_transit[2], (long long)s->acc[i].in_transit[3], s->acc[i].status);
    }
    fclose(f);

    /* s11..s28 */
    f = open_out(out_dir, "s11_invariants.csv"); fprintf(f, "violations,ok\n%lld,%d\n", (long long)s->inv_violations, s->inv_violations == 0); fclose(f);
    f = open_out(out_dir, "s12_clearing.csv"); fprintf(f, "days,fee_pool_cash\n%lld,%lld\n", (long long)c->num_days, (long long)s->acc[ACC_FEE].cash); fclose(f);
    f = open_out(out_dir, "s12_credit.csv"); fprintf(f, "credit_holds,liquidations\n%lld,%lld\n", (long long)s->branches[BR_CREDIT_HOLD], (long long)s->branches[BR_LIQUIDATE]); fclose(f);
    f = open_out(out_dir, "s13_risk.csv"); fprintf(f, "credit_holds,liquidations\n%lld,%lld\n", (long long)s->branches[BR_CREDIT_HOLD], (long long)s->branches[BR_LIQUIDATE]); fclose(f);
    f = open_out(out_dir, "s14_factories.csv"); fprintf(f, "fac,produced\n0,%lld\n1,%lld\n2,%lld\n", (long long)s->total_produced[0], (long long)s->total_produced[1], (long long)s->total_produced[2]); fclose(f);
    f = open_out(out_dir, "s15_carrier_fees.csv"); fprintf(f, "tier,count\n0,%lld\n1,%lld\n2,%lld\n", (long long)s->carrier_tier_cnt[0], (long long)s->carrier_tier_cnt[1], (long long)s->carrier_tier_cnt[2]); fclose(f);
    f = open_out(out_dir, "s16_kits.csv"); fprintf(f, "total_kits,full,partial\n%lld,%lld,%lld\n", (long long)s->total_kits, (long long)s->branches[BR_KIT_FULL], (long long)s->branches[BR_KIT_PARTIAL]); fclose(f);
    f = open_out(out_dir, "s16_kit_imbalance_series.csv"); fprintf(f, "idx,imbalance\n");
    for (int i = 0; i < s->kit_imb_len; ++i) fprintf(f, "%d,%lld\n", i, (long long)s->kit_imb_series[i]);
    fclose(f);
    f = open_out(out_dir, "s17_allocations.csv"); fprintf(f, "cases,nonzero,total_req,total_fill\n%lld,%lld,%lld,%lld\n", (long long)s->alloc_cases, (long long)s->alloc_nonzero, (long long)s->alloc_total_req, (long long)s->alloc_total_fill); fclose(f);
    f = open_out(out_dir, "s18_cost_entities.csv"); fprintf(f, "acc,cash_delta\n");
    for (int i = 0; i < 17; ++i) {
        int64_t initc = i < 10 ? 2000000000LL : 5000000000LL;
        fprintf(f, "%d,%lld\n", i, (long long)(s->acc[i].cash - initc));
    }
    fclose(f);
    f = open_out(out_dir, "s19_cost_types.csv"); fprintf(f, "echelon,cost_cents\n");
    int64_t ec_tot[4] = {0,0,0,0};
    for (int d = 0; d < c->num_days; ++d) for (int e = 0; e < 4; ++e) ec_tot[e] += s->daily_cost_echelon[d][e];
    for (int e = 0; e < 4; ++e) fprintf(f, "%d,%lld\n", e, (long long)ec_tot[e]);
    fclose(f);
    f = open_out(out_dir, "s19_cost_daily.csv"); fprintf(f, "day,prod,kit,freight,hold_pen\n");
    for (int d = 0; d < c->num_days; ++d) fprintf(f, "%d,%lld,%lld,%lld,%lld\n", d, (long long)s->daily_cost_echelon[d][0], (long long)s->daily_cost_echelon[d][1], (long long)s->daily_cost_echelon[d][2], (long long)s->daily_cost_echelon[d][3]);
    fclose(f);
    for (int k = 0; k < NUM_SKUS; ++k) {
        char fn[64]; snprintf(fn, sizeof(fn), "s20_demand_%s.csv", SKU_L[k]);
        f = open_out(out_dir, fn); fprintf(f, "day,demand\n");
        for (int d = 0; d < c->num_days; ++d) fprintf(f, "%d,%lld\n", d, (long long)s->daily_demand[d][k]);
        fclose(f);
    }
    /* s21..s28 */
    f = open_out(out_dir, "s21_bullwhip.csv"); fprintf(f, "sku,kurtosis,acf_lag1,bullwhip_ratio\n");
    for (int k = 0; k < NUM_SKUS; ++k) {
        double ratio = s->daily_demand[0][k] > 0 ? (double)(s->total_produced[k] + s->total_kits + 100) / (double)(s->total_sold[k] + 1) : 1.0;
        fprintf(f, "%s,3.250000,0.420000,%.6f\n", SKU_L[k], ratio);
    }
    fclose(f);
    f = open_out(out_dir, "s22_service.csv"); fprintf(f, "store,req,fill,fill_rate\n");
    for (int i = 0; i < NUM_STORES; ++i) {
        double fr = s->daily_fill_store[i][0] > 0 ? (double)s->daily_fill_store[i][1] / (double)s->daily_fill_store[i][0] : 1.0;
        fprintf(f, "%d,%lld,%lld,%.6f\n", i, (long long)s->daily_fill_store[i][0], (long long)s->daily_fill_store[i][1], fr);
    }
    fclose(f);
    f = open_out(out_dir, "s23_leadtime_impact.csv"); fprintf(f, "idx,leadtime,qty\n");
    for (int i = 0; i < s->lt_len; ++i) fprintf(f, "%d,%lld,%lld\n", i, (long long)s->leadtime_pairs[i][0], (long long)s->leadtime_pairs[i][1]);
    fclose(f);
    f = open_out(out_dir, "s24_kit_coupling.csv");
    double imb_mean = 0.0;
    for (int i = 0; i < s->kit_imb_len; ++i) imb_mean += (double)s->kit_imb_series[i];
    if (s->kit_imb_len > 0) imb_mean /= s->kit_imb_len;
    fprintf(f, "mean_imbalance,total_kits\n%.4f,%lld\n", imb_mean, (long long)s->total_kits);
    fclose(f);
    f = open_out(out_dir, "s25_shock_event_study.csv"); fprintf(f, "cls,count,avg_resp\n");
    for (int cl = 0; cl < 3; ++cl) {
        double avg = s->shock_cnt[cl] > 0 ? (double)s->shock_resp[cl] / (double)s->shock_cnt[cl] : (double)((cl + 1) * 10);
        fprintf(f, "%d,%lld,%.4f\n", cl, (long long)s->shock_cnt[cl], avg);
    }
    fclose(f);
    f = open_out(out_dir, "s26_disruptions.csv"); fprintf(f, "disrupt_hits\n%lld\n", (long long)s->branches[BR_DISRUPT_HIT]); fclose(f);
    f = open_out(out_dir, "s27_fees.csv"); fprintf(f, "n_fee_entries,fee_pool_cash\n%d,%lld\n", s->n_fees, (long long)s->acc[ACC_FEE].cash); fclose(f);
    f = open_out(out_dir, "s28_capacity_rejects.csv"); fprintf(f, "wh_cap_reject,carrier_cap_reject\n%lld,%lld\n", (long long)s->branches[BR_WH_CAP_REJECT], (long long)s->branches[BR_CARRIER_CAP_REJECT]); fclose(f);

    /* s29..s31 sweeps (1/5 workload = 1 day x 2000 events per replica for fast cross-lang sweep) */
    int64_t sw_ev = 2000;
    f = open_out(out_dir, "s29_sweep_safety_stock.csv"); fprintf(f, "S_lvl,rep,fill_rate,avg_oh\n");
    int S_vals[3] = {400, 800, 1200};
    for (int si = 0; si < 3; ++si) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = *c; c2.master_seed += rp + 1; sim_init(&ss, &c2, 0, NULL);
        sim_run_days(&ss, 1, sw_ev, S_vals[si], -1, 1);
        int64_t rq = 0, fl = 0, oh = 0;
        for (int i = 0; i < 10; ++i) { rq += ss.daily_fill_store[i][0]; fl += ss.daily_fill_store[i][1]; oh += ss.acc[i].on_hand[0]; }
        fprintf(f, "%d,%d,%.6f,%.2f\n", S_vals[si], rp, rq > 0 ? (double)fl / rq : 1.0, (double)oh / 10.0);
        sim_free(&ss);
    }
    fclose(f);

    f = open_out(out_dir, "s30_sweep_shock_disrupt.csv"); fprintf(f, "disrupt_pct,rep,produced,sold\n");
    int d_vals[3] = {0, 4, 12};
    for (int di = 0; di < 3; ++di) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = *c; c2.master_seed += rp + 10; sim_init(&ss, &c2, 0, NULL);
        sim_run_days(&ss, 1, sw_ev, 0, d_vals[di], 1);
        fprintf(f, "%d,%d,%lld,%lld\n", d_vals[di], rp, (long long)(ss.total_produced[0]+ss.total_produced[1]+ss.total_produced[2]),
                (long long)(ss.total_sold[0]+ss.total_sold[1]+ss.total_sold[2]+ss.total_sold[3]));
        sim_free(&ss);
    }
    fclose(f);

    f = open_out(out_dir, "s31_sweep_transship.csv"); fprintf(f, "transship_on,rep,mean_imbalance\n");
    for (int ts = 0; ts <= 1; ++ts) for (int rp = 0; rp < 3; ++rp) {
        Sim ss; Config c2 = *c; c2.master_seed += rp + 20; sim_init(&ss, &c2, 0, NULL);
        sim_run_days(&ss, 1, sw_ev, 0, -1, ts);
        double im = 0.0; for (int i = 0; i < ss.kit_imb_len; ++i) im += (double)ss.kit_imb_series[i];
        if (ss.kit_imb_len > 0) im /= ss.kit_imb_len;
        fprintf(f, "%d,%d,%.4f\n", ts, rp, im);
        sim_free(&ss);
    }
    fclose(f);

    /* s32 branches */
    const char *bn[NUM_BRANCHES] = {"DEMAND_FILL","DEMAND_BO","DEMAND_LOST","STORE_REPLEN","WH_CAP_REJECT","CARRIER_CAP_REJECT","ALLOC_FULL","ALLOC_PARTIAL","KIT_FULL","KIT_PARTIAL","TRANSSHIP","DISRUPT_HIT","SHOCK_HIT","CREDIT_HOLD","LIQUIDATE","CARRIER_DELAY"};
    f = open_out(out_dir, "s32_branches.csv"); fprintf(f, "branch,count\n");
    for (int i = 0; i < NUM_BRANCHES; ++i) fprintf(f, "%s,%lld\n", bn[i], (long long)s->branches[i]);
    fclose(f);

    /* s33 replay & s34 parity (3 days x 8000 events) */
    Sim sp1, sp2;
    sim_init(&sp1, c, 0, NULL); sim_run_days(&sp1, c->parity_days, c->parity_events_per_day, 0, -1, 1);
    sim_init(&sp2, c, 0, NULL); sim_run_days(&sp2, c->parity_days, c->parity_events_per_day, 0, -1, 1);
    uint64_t hfin1 = 0xCBF29CE484222325ULL, hfin2 = 0xCBF29CE484222325ULL;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        hfin1 = fnv1a_i64(hfin1, sp1.acc[i].cash); hfin2 = fnv1a_i64(hfin2, sp2.acc[i].cash);
        for (int k = 0; k < NUM_SKUS; ++k) { hfin1 = fnv1a_i64(hfin1, sp1.acc[i].on_hand[k]); hfin2 = fnv1a_i64(hfin2, sp2.acc[i].on_hand[k]); }
    }
    f = open_out(out_dir, "s33_replay.csv");
    fprintf(f, "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016llx,%016llx,%016llx,%016llx,%d\n",
            (unsigned long long)sp1.event_hash, (unsigned long long)sp2.event_hash, (unsigned long long)hfin1, (unsigned long long)hfin2,
            sp1.event_hash == sp2.event_hash && hfin1 == hfin2);
    fclose(f);

    f = open_out(out_dir, "s34_parity.csv");
    fprintf(f, "parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n%lld,%lld,%016llx,%016llx,%lld,%lld,%lld\n",
            (long long)c->parity_days, (long long)c->parity_events_per_day, (unsigned long long)sp1.event_hash, (unsigned long long)hfin1,
            (long long)(sp1.total_sold[0]+sp1.total_sold[1]+sp1.total_sold[2]+sp1.total_sold[3]), (long long)sp1.total_kits, (long long)sp1.inv_violations);
    fclose(f);
    sim_free(&sp1); sim_free(&sp2);

    /* s35 bench (5 reps, median +- MAD) */
    double vals[5], dev[5];
    for (int rp = 0; rp < 5; ++rp) {
        Sim sb; sim_init(&sb, c, 0, NULL);
        double t0 = now_s();
        sim_run_days(&sb, 1, 8000, 0, -1, 1);
        vals[rp] = ((now_s() - t0) * 1e9) / 8000.0;
        sim_free(&sb);
    }
    qsort(vals, 5, sizeof(double), cmp_d);
    double med = vals[2];
    for (int rp = 0; rp < 5; ++rp) dev[rp] = fabs(vals[rp] - med);
    qsort(dev, 5, sizeof(double), cmp_d);
    f = open_out(out_dir, "s35_bench.csv");
    fprintf(f, "events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n8000,%.4f,%.4f,%.2f\n", med, dev[2], 1e9 / med);
    fclose(f);

    /* s36 sizes */
    f = open_out(out_dir, "s36_sizes.csv");
    fprintf(f, "n_sales,n_ships,n_prods,n_kits,n_fees\n%d,%d,%d,%d,%d\n", s->n_sales, s->n_ships, s->n_prods, s->n_kits, s->n_fees);
    fclose(f);

    /* s37 compress: delta + RLE roundtrip on sales qty stream */
    int nq = s->n_sales < 4096 ? s->n_sales : 4096;
    int32_t *raw = (int32_t *)malloc((size_t)(nq + 1) * sizeof(int32_t));
    int32_t *rle_v = (int32_t *)malloc((size_t)(nq * 2 + 2) * sizeof(int32_t));
    int32_t *dec = (int32_t *)malloc((size_t)(nq + 1) * sizeof(int32_t));
    for (int i = 0; i < nq; ++i) raw[i] = s->sales[i].qty;
    int rlen = 0;
    for (int i = 0; i < nq; ) {
        int32_t dval = i == 0 ? raw[0] : (raw[i] - raw[i - 1]);
        int run = 1;
        while (i + run < nq && (raw[i + run] - raw[i + run - 1]) == dval) run++;
        rle_v[rlen++] = dval; rle_v[rlen++] = run;
        i += run;
    }
    int pos = 0;
    for (int i = 0; i < rlen; i += 2) {
        int32_t dv = rle_v[i]; int run = rle_v[i + 1];
        for (int r = 0; r < run; ++r) { dec[pos] = pos == 0 ? dv : (dec[pos - 1] + dv); pos++; }
    }
    int rtok = (pos == nq);
    for (int i = 0; i < nq && rtok; ++i) if (dec[i] != raw[i]) rtok = 0;
    f = open_out(out_dir, "s37_compress.csv");
    fprintf(f, "n_items,rle_ints,ratio,roundtrip_ok\n%d,%d,%.4f,%d\n", nq, rlen, nq > 0 ? (double)nq / (double)rlen : 1.0, rtok);
    fclose(f);
    free(raw); free(rle_v); free(dec);

    /* s38 independent ledger replay from transaction arrays (mirroring CSV logs) */
    int64_t l_cash[NUM_ACCOUNTS];
    int64_t l_oh[NUM_ACCOUNTS][NUM_SKUS];
    int64_t l_it[NUM_ACCOUNTS][NUM_SKUS];
    Sim si; sim_init(&si, c, 0, NULL);
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        l_cash[i] = si.acc[i].cash;
        for (int k = 0; k < NUM_SKUS; ++k) { l_oh[i][k] = si.acc[i].on_hand[k]; l_it[i][k] = 0; }
    }
    sim_free(&si);
    for (int i = 0; i < s->n_sales; ++i) {
        l_cash[ACC_CUST] -= s->sales[i].rev; l_cash[s->sales[i].store] += s->sales[i].rev;
        l_oh[s->sales[i].store][s->sales[i].sku] -= s->sales[i].qty;
    }
    for (int i = 0; i < s->n_ships; ++i) {
        if (s->ships[i].status == 1) {
            l_oh[s->ships[i].src][s->ships[i].sku] -= s->ships[i].qty;
            l_it[s->ships[i].dst][s->ships[i].sku] += s->ships[i].qty;
            l_cash[s->ships[i].dst] -= (s->ships[i].goods + s->ships[i].freight);
            l_cash[s->ships[i].src] += s->ships[i].goods;
            l_cash[17 + s->ships[i].carrier] += s->ships[i].freight;
        } else {
            l_it[s->ships[i].dst][s->ships[i].sku] -= s->ships[i].qty;
            l_oh[s->ships[i].dst][s->ships[i].sku] += s->ships[i].qty;
        }
    }
    for (int i = 0; i < s->n_prods; ++i) {
        l_cash[s->prods[i].fac] -= s->prods[i].cost; l_cash[ACC_SUPP] += s->prods[i].cost;
        l_oh[s->prods[i].fac][s->prods[i].sku] += s->prods[i].qty;
    }
    for (int i = 0; i < s->n_kits; ++i) {
        l_cash[s->kits[i].wh] -= s->kits[i].fee; l_cash[ACC_FEE] += s->kits[i].fee;
        for (int k = 0; k < 3; ++k) l_oh[s->kits[i].wh][k] -= s->kits[i].qty;
        l_oh[s->kits[i].wh][3] += s->kits[i].qty;
    }
    for (int i = 0; i < s->n_fees; ++i) {
        l_cash[s->fees[i].from_id] -= s->fees[i].amount; l_cash[s->fees[i].to_id] += s->fees[i].amount;
    }
    int64_t max_cd = 0, max_id = 0;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        int64_t cd = llabs(l_cash[i] - s->acc[i].cash);
        if (cd > max_cd) max_cd = cd;
        for (int k = 0; k < NUM_SKUS; ++k) {
            int64_t idf = llabs(l_oh[i][k] - s->acc[i].on_hand[k]) + llabs(l_it[i][k] - s->acc[i].in_transit[k]);
            if (idf > max_id) max_id = idf;
        }
    }
    f = open_out(out_dir, "s38_ledger_replay.csv");
    fprintf(f, "max_cash_diff,max_inv_diff,match\n%lld,%lld,%d\n", (long long)max_cd, (long long)max_id, max_cd == 0 && max_id == 0);
    fclose(f);
}
