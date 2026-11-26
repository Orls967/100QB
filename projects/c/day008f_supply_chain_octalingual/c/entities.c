#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int load_config(const char *path, Config *c) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    memset(c, 0, sizeof(*c));
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        long long v = atoll(eq + 1);
        if (!strcmp(line, "master_seed")) c->master_seed = v;
        else if (!strcmp(line, "num_days")) c->num_days = v;
        else if (!strcmp(line, "parity_days")) c->parity_days = v;
        else if (!strcmp(line, "parity_events_per_day")) c->parity_events_per_day = v;
        else if (!strcmp(line, "ticks_per_day")) c->ticks_per_day = v;
        else if (!strcmp(line, "num_skus")) c->num_skus = v;
        else if (!strcmp(line, "num_factories")) c->num_factories = v;
        else if (!strcmp(line, "num_warehouses")) c->num_warehouses = v;
        else if (!strcmp(line, "num_stores")) c->num_stores = v;
        else if (!strcmp(line, "num_carriers")) c->num_carriers = v;
        else if (!strcmp(line, "ewma_alpha_milli")) c->ewma_alpha_milli = v;
        else if (!strcmp(line, "store_s_base")) c->store_s_base = v;
        else if (!strcmp(line, "store_S_base")) c->store_S_base = v;
        else if (!strcmp(line, "wh_s_base")) c->wh_s_base = v;
        else if (!strcmp(line, "wh_S_base")) c->wh_S_base = v;
        else if (!strcmp(line, "wh_cap")) c->wh_cap = v;
        else if (!strcmp(line, "factory_cap_norm")) c->factory_cap_norm = v;
        else if (!strcmp(line, "factory_cap_ot")) c->factory_cap_ot = v;
        else if (!strcmp(line, "carrier_cap_daily")) c->carrier_cap_daily = v;
        else if (!strcmp(line, "kit_fee")) c->kit_fee = v;
        else if (!strcmp(line, "holding_cost_per_day")) c->holding_cost_per_day = v;
        else if (!strcmp(line, "backorder_penalty")) c->backorder_penalty = v;
        else if (!strcmp(line, "disrupt_prob_pct")) c->disrupt_prob_pct = v;
        else if (!strcmp(line, "shock_prob_pct")) c->shock_prob_pct = v;
        else if (!strcmp(line, "credit_margin_bps")) c->credit_margin_bps = v;
        else if (!strcmp(line, "k_chk")) c->k_chk = v;
    }
    fclose(f);
    return 1;
}

void sim_init(Sim *s, const Config *c, int32_t log_en, const char *out_dir) {
    memset(s, 0, sizeof(*s));
    s->cfg = *c;
    s->log_enabled = log_en;
    s->next_seq = 1;
    s->event_hash = 0xCBF29CE484222325ULL;
    heap_init(&s->heap);
    for (int i = 0; i < NUM_ACCOUNTS + 2; ++i) {
        uint64_t x = (uint64_t)c->master_seed ^ fnv1a_u32pair((uint32_t)(i / 10), (uint32_t)i);
        x256_seed(&s->rng[i], splitmix64(&x));
    }
    for (int i = 0; i < NUM_STORES; ++i) {
        s->acc[i].cash = 2000000000LL;
        s->acc[i].credit_limit = 500000000LL;
        for (int k = 0; k < 3; ++k) s->acc[i].on_hand[k] = 400;
        s->acc[i].on_hand[3] = 120;
        for (int k = 0; k < NUM_SKUS; ++k) s->acc[i].ewma_milli[k] = 12000;
        s->acc[i].safety_mul = 100;
    }
    for (int i = 10; i < 14; ++i) {
        s->acc[i].cash = 5000000000LL;
        s->acc[i].credit_limit = 1000000000LL;
        for (int k = 0; k < 3; ++k) s->acc[i].on_hand[k] = 2500;
        s->acc[i].on_hand[3] = 800;
    }
    for (int i = 14; i < 17; ++i) {
        s->acc[i].cash = 5000000000LL;
        s->acc[i].credit_limit = 2000000000LL;
        s->acc[i].on_hand[i - 14] = 4000;
    }
    for (int i = 17; i < 20; ++i) {
        s->acc[i].cash = 500000000LL;
        s->acc[i].credit_limit = 200000000LL;
    }
    s->acc[ACC_CUST].cash = 50000000000LL;
    s->acc[ACC_SUPP].cash = 50000000000LL;
    s->acc[ACC_FEE].cash  = 1000000000LL;
    for (int i = 0; i < NUM_ACCOUNTS; ++i) {
        s->init_cash += s->acc[i].cash;
        for (int k = 0; k < NUM_SKUS; ++k) s->init_stock[k] += s->acc[i].on_hand[k];
    }
    if (log_en && out_dir) {
        char p[512];
        snprintf(p, sizeof(p), "%s/orders_log.csv", out_dir); s->f_orders = fopen(p, "w");
        fprintf(s->f_orders, "day,tick,seq,etype,src,dst,sku,qty\n");
        snprintf(p, sizeof(p), "%s/allocations.csv", out_dir); s->f_alloc = fopen(p, "w");
        fprintf(s->f_alloc, "day,wh,sku,store,req,fill\n");
        snprintf(p, sizeof(p), "%s/disruptions.csv", out_dir); s->f_disrupt = fopen(p, "w");
        fprintf(s->f_disrupt, "day,tick,fac,state,dur\n");
        snprintf(p, sizeof(p), "%s/credit_calls.csv", out_dir); s->f_credit = fopen(p, "w");
        fprintf(s->f_credit, "day,acc,cash,req_margin,status\n");
        snprintf(p, sizeof(p), "%s/risk_actions.csv", out_dir); s->f_risk = fopen(p, "w");
        fprintf(s->f_risk, "day,tick,acc,action,amount\n");
        snprintf(p, sizeof(p), "%s/transshipments.csv", out_dir); s->f_trans = fopen(p, "w");
        fprintf(s->f_trans, "day,tick,src,dst,sku,qty\n");
        snprintf(p, sizeof(p), "%s/shocks_log.csv", out_dir); s->f_shocks = fopen(p, "w");
        fprintf(s->f_shocks, "day,tick,cls,store,mag\n");
        snprintf(p, sizeof(p), "%s/accounts_series.csv", out_dir); s->f_acc_ser = fopen(p, "w");
        fprintf(s->f_acc_ser, "day,acc,cash,oh_tot,it_tot,status\n");
    }
}

void sim_free(Sim *s) {
    heap_free(&s->heap);
    free(s->sales); free(s->ships); free(s->prods); free(s->kits); free(s->fees);
    if (s->f_orders) fclose(s->f_orders);
    if (s->f_alloc) fclose(s->f_alloc);
    if (s->f_disrupt) fclose(s->f_disrupt);
    if (s->f_credit) fclose(s->f_credit);
    if (s->f_risk) fclose(s->f_risk);
    if (s->f_trans) fclose(s->f_trans);
    if (s->f_shocks) fclose(s->f_shocks);
    if (s->f_acc_ser) fclose(s->f_acc_ser);
}
