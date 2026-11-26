#ifndef SC_COMMON_H
#define SC_COMMON_H

#include <stdint.h>
#include <stdio.h>

#define NUM_SKUS 4
#define NUM_FACTORIES 3
#define NUM_WAREHOUSES 4
#define NUM_STORES 10
#define NUM_CARRIERS 3
#define NUM_ACCOUNTS 23
#define ACC_CUST 20
#define ACC_SUPP 21
#define ACC_FEE 22

#define ST_ACTIVE 0
#define ST_CREDIT_HOLD 1
#define ST_LIQUIDATION 2
#define ST_SUSPENDED 3

#define PLANT_OPEN 0
#define PLANT_DISRUPTED 1
#define PLANT_RECOVERY 2

#define EV_DEMAND 0
#define EV_REPLENISH 1
#define EV_DELIVER 2
#define EV_PRODUCE 3
#define EV_KIT 4
#define EV_TRANSSHIP 5
#define EV_SHOCK 6
#define EV_DISRUPT 7
#define EV_RISK_CHK 8

#define BR_DEMAND_FILL 0
#define BR_DEMAND_BO 1
#define BR_DEMAND_LOST 2
#define BR_STORE_REPLEN 3
#define BR_WH_CAP_REJECT 4
#define BR_CARRIER_CAP_REJECT 5
#define BR_ALLOC_FULL 6
#define BR_ALLOC_PARTIAL 7
#define BR_KIT_FULL 8
#define BR_KIT_PARTIAL 9
#define BR_TRANSSHIP 10
#define BR_DISRUPT_HIT 11
#define BR_SHOCK_HIT 12
#define BR_CREDIT_HOLD 13
#define BR_LIQUIDATE 14
#define BR_CARRIER_DELAY 15
#define NUM_BRANCHES 16

typedef struct {
    uint64_t s[4];
} Xoshiro256;

typedef struct {
    int64_t master_seed, num_days, parity_days, parity_events_per_day;
    int64_t ticks_per_day, num_skus, num_factories, num_warehouses;
    int64_t num_stores, num_carriers, ewma_alpha_milli;
    int64_t store_s_base, store_S_base, wh_s_base, wh_S_base, wh_cap;
    int64_t factory_cap_norm, factory_cap_ot, carrier_cap_daily;
    int64_t kit_fee, holding_cost_per_day, backorder_penalty;
    int64_t disrupt_prob_pct, shock_prob_pct, credit_margin_bps, k_chk;
} Config;

typedef struct {
    int64_t cash;
    int64_t on_hand[NUM_SKUS];
    int64_t in_transit[NUM_SKUS];
    int64_t backorder[NUM_SKUS];
    int64_t credit_limit;
    int32_t status;
    int32_t ewma_milli[NUM_SKUS];
    int32_t safety_mul;
    int32_t daily_vol;
    int32_t plant_state;
    int32_t disrupt_rem;
} Account;

typedef struct {
    int64_t time;
    uint64_t seq;
    int32_t etype, src, dst, sku, qty, aux;
} SimEvent;

typedef struct {
    SimEvent *data;
    int cap, len;
} MinHeap;

typedef struct {
    int64_t day, tick;
    int32_t store, sku, qty;
    int64_t rev;
} SaleRec;

typedef struct {
    int64_t day, tick;
    int32_t src, dst, carrier, sku, qty, status;
    int64_t goods, freight;
} ShipRec;

typedef struct {
    int64_t day, tick;
    int32_t fac, sku, qty;
    int64_t cost;
} ProdRec;

typedef struct {
    int64_t day, tick;
    int32_t wh, qty, imbalance;
    int64_t fee;
} KitRec;

typedef struct {
    int64_t day, tick;
    int32_t from_id, to_id, kind;
    int64_t amount;
} FeeRec;

typedef struct {
    Config cfg;
    Account acc[NUM_ACCOUNTS];
    Xoshiro256 rng[NUM_ACCOUNTS + 2];
    MinHeap heap;
    uint64_t next_seq;
    int64_t init_cash;
    int64_t init_stock[NUM_SKUS];
    int64_t total_produced[NUM_SKUS];
    int64_t total_sold[NUM_SKUS];
    int64_t total_kits;
    int64_t branches[NUM_BRANCHES];
    int64_t carrier_day_vol[NUM_CARRIERS];
    int64_t carrier_tier_cnt[3];
    int64_t alloc_total_req, alloc_total_fill, alloc_cases, alloc_nonzero;
    int64_t inv_violations;
    int32_t log_enabled;
    int64_t daily_demand[6][NUM_SKUS];
    int64_t daily_cost_echelon[6][4];
    int64_t daily_fill_store[10][2];
    int64_t shock_window[25][2];
    int64_t shock_cnt[3];
    int64_t shock_resp[3];
    int64_t day_inv_snap[6][NUM_SKUS];
    int64_t kit_imb_series[64];
    int kit_imb_len;
    int64_t leadtime_pairs[256][2];
    int lt_len;
    SaleRec *sales; int n_sales, cap_sales;
    ShipRec *ships; int n_ships, cap_ships;
    ProdRec *prods; int n_prods, cap_prods;
    KitRec *kits;   int n_kits, cap_kits;
    FeeRec *fees;   int n_fees, cap_fees;
    FILE *f_orders, *f_alloc, *f_disrupt, *f_credit, *f_risk, *f_trans, *f_shocks, *f_acc_ser;
    uint64_t event_hash;
} Sim;

extern const int64_t SKU_RETAIL[NUM_SKUS];
extern const int64_t SKU_WHOLESALE[NUM_SKUS];
extern const int64_t SKU_PROD[NUM_SKUS];
extern const int32_t STORE_WH[NUM_STORES];
extern const int32_t CARRIER_BASE_LT[NUM_CARRIERS];

uint64_t splitmix64(uint64_t *x);
uint64_t fnv1a_bytes(const void *buf, size_t n);
uint64_t fnv1a_u32pair(uint32_t a, uint32_t b);
uint64_t fnv1a_i64(uint64_t h, int64_t v);
void x256_seed(Xoshiro256 *st, uint64_t seed);
uint64_t x256_next(Xoshiro256 *st);
int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi);

void heap_init(MinHeap *h);
void heap_free(MinHeap *h);
void heap_push(MinHeap *h, SimEvent ev);
SimEvent heap_pop(MinHeap *h);

int load_config(const char *path, Config *c);
void sim_init(Sim *s, const Config *c, int32_t log_en, const char *out_dir);
void sim_free(Sim *s);
int check_invariants(Sim *s);
void transfer_cash(Sim *s, int32_t from_id, int32_t to_id, int64_t amt, int32_t kind, int64_t day, int64_t tick);
int64_t carrier_rate(int64_t day_vol);
void run_alloc_auction(Sim *s, int32_t wh, int32_t sku, const int32_t *stores, const int64_t *reqs, const int64_t *seqs, const int32_t *prio, int n, int64_t *fills, int64_t day, int64_t tick);
void sim_run_days(Sim *s, int64_t days, int64_t ev_per_day, int safety_override, int disrupt_override, int transship_on);
void sim_eod_clearing(Sim *s, int64_t day);
void sim_risk_check(Sim *s, int64_t day, int64_t tick);

void write_s01_to_s09(const Config *c, const char *out_dir);
void write_s11_to_s38(Sim *s, const Config *c, const char *out_dir, double *stage_s);
void write_s39_svgs(Sim *s, const char *out_dir);
void write_s40_manifest(const char *out_dir, const char *toolchain, const double *stage_s, int64_t ev_total);

#endif
