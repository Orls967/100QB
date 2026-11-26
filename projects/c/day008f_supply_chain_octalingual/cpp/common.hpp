#pragma once
#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>

constexpr int NUM_SKUS = 4;
constexpr int NUM_FACTORIES = 3;
constexpr int NUM_WAREHOUSES = 4;
constexpr int NUM_STORES = 10;
constexpr int NUM_CARRIERS = 3;
constexpr int NUM_ACCOUNTS = 23;
constexpr int ACC_CUST = 20;
constexpr int ACC_SUPP = 21;
constexpr int ACC_FEE  = 22;

enum { ST_ACTIVE = 0, ST_CREDIT_HOLD = 1, ST_LIQUIDATION = 2, ST_SUSPENDED = 3 };
enum { PLANT_OPEN = 0, PLANT_DISRUPTED = 1, PLANT_RECOVERY = 2 };
enum { EV_DEMAND = 0, EV_REPLENISH = 1, EV_DELIVER = 2, EV_PRODUCE = 3, EV_KIT = 4, EV_TRANSSHIP = 5, EV_SHOCK = 6, EV_DISRUPT = 7, EV_RISK_CHK = 8 };
enum {
    BR_DEMAND_FILL = 0, BR_DEMAND_BO, BR_DEMAND_LOST, BR_STORE_REPLEN,
    BR_WH_CAP_REJECT, BR_CARRIER_CAP_REJECT, BR_ALLOC_FULL, BR_ALLOC_PARTIAL,
    BR_KIT_FULL, BR_KIT_PARTIAL, BR_TRANSSHIP, BR_DISRUPT_HIT,
    BR_SHOCK_HIT, BR_CREDIT_HOLD, BR_LIQUIDATE, BR_CARRIER_DELAY, NUM_BRANCHES
};

struct Xoshiro256 { uint64_t s[4]; };

struct Config {
    int64_t master_seed, num_days, parity_days, parity_events_per_day;
    int64_t ticks_per_day, num_skus, num_factories, num_warehouses;
    int64_t num_stores, num_carriers, ewma_alpha_milli;
    int64_t store_s_base, store_S_base, wh_s_base, wh_S_base, wh_cap;
    int64_t factory_cap_norm, factory_cap_ot, carrier_cap_daily;
    int64_t kit_fee, holding_cost_per_day, backorder_penalty;
    int64_t disrupt_prob_pct, shock_prob_pct, credit_margin_bps, k_chk;
};

struct Account {
    int64_t cash = 0;
    int64_t on_hand[NUM_SKUS]{};
    int64_t in_transit[NUM_SKUS]{};
    int64_t backorder[NUM_SKUS]{};
    int64_t credit_limit = 0;
    int32_t status = 0;
    int32_t ewma_milli[NUM_SKUS]{};
    int32_t safety_mul = 0;
    int32_t daily_vol = 0;
    int32_t plant_state = 0;
    int32_t disrupt_rem = 0;
};

struct SimEvent {
    int64_t time;
    uint64_t seq;
    int32_t etype, src, dst, sku, qty, aux;
};

struct SaleRec { int64_t day, tick; int32_t store, sku, qty; int64_t rev; };
struct ShipRec { int64_t day, tick; int32_t src, dst, carrier, sku, qty, status; int64_t goods, freight; };
struct ProdRec { int64_t day, tick; int32_t fac, sku, qty; int64_t cost; };
struct KitRec  { int64_t day, tick; int32_t wh, qty, imbalance; int64_t fee; };
struct FeeRec  { int64_t day, tick; int32_t from_id, to_id, kind; int64_t amount; };

struct Sim {
    Config cfg{};
    Account acc[NUM_ACCOUNTS]{};
    Xoshiro256 rng[NUM_ACCOUNTS + 2]{};
    std::vector<SimEvent> heap;
    uint64_t next_seq = 1;
    int64_t init_cash = 0;
    int64_t init_stock[NUM_SKUS]{};
    int64_t total_produced[NUM_SKUS]{};
    int64_t total_sold[NUM_SKUS]{};
    int64_t total_kits = 0;
    int64_t branches[NUM_BRANCHES]{};
    int64_t carrier_day_vol[NUM_CARRIERS]{};
    int64_t carrier_tier_cnt[3]{};
    int64_t alloc_total_req = 0, alloc_total_fill = 0, alloc_cases = 0, alloc_nonzero = 0;
    int64_t inv_violations = 0;
    int32_t log_enabled = 0;
    int64_t daily_demand[6][NUM_SKUS]{};
    int64_t daily_cost_echelon[6][4]{};
    int64_t daily_fill_store[10][2]{};
    int64_t shock_cnt[3]{};
    int64_t shock_resp[3]{};
    int64_t day_inv_snap[6][NUM_SKUS]{};
    int64_t kit_imb_series[64]{};
    int kit_imb_len = 0;
    int64_t leadtime_pairs[256][2]{};
    int lt_len = 0;
    int64_t bin_demand[NUM_SKUS][100]{};
    int64_t bin_orders[NUM_SKUS][100]{};
    int64_t bin_store_dem[NUM_STORES][100]{};
    int32_t shock_ev_bin[128]{};
    int32_t shock_ev_cls[128]{};
    int32_t shock_ev_st[128]{};
    int n_shock_evs = 0;
    int64_t n_transships = 0;
    int64_t dem_events_tot = 0, dem_bo_events = 0;
    int64_t dem_units_tot = 0, dem_units_filled = 0;
    int64_t carrier_cap_tot = 0, carrier_vol_tot = 0;
    int64_t cap_scale = 1;
    std::vector<SaleRec> sales;
    std::vector<ShipRec> ships;
    std::vector<ProdRec> prods;
    std::vector<KitRec> kits;
    std::vector<FeeRec> fees;
    FILE *f_orders = nullptr, *f_alloc = nullptr, *f_disrupt = nullptr, *f_credit = nullptr;
    FILE *f_risk = nullptr, *f_trans = nullptr, *f_shocks = nullptr, *f_acc_ser = nullptr;
    uint64_t event_hash = 0xCBF29CE484222325ULL;
};

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
void heap_push(std::vector<SimEvent> &h, SimEvent ev);
SimEvent heap_pop(std::vector<SimEvent> &h);

bool load_config(const char *path, Config &c);
void sim_init(Sim &s, const Config &c, int32_t log_en, const char *out_dir);
void sim_init_scaled(Sim &s, const Config &c, int32_t log_en, const char *out_dir, int64_t ev_per_day);
void sim_free(Sim &s);
int check_invariants(Sim &s);
void transfer_cash(Sim &s, int32_t from_id, int32_t to_id, int64_t amt, int32_t kind, int64_t day, int64_t tick);
int64_t carrier_rate(int64_t day_vol);
void run_alloc_auction(Sim &s, int32_t wh, int32_t sku, const int32_t *stores, const int64_t *reqs, const int64_t *seqs, const int32_t *prio, int n, int64_t *fills);
void sim_run_days(Sim &s, int64_t days, int64_t ev_per_day, int safety_override, int disrupt_override, int transship_on);
void sim_eod_clearing(Sim &s, int64_t day);
void sim_risk_check(Sim &s, int64_t day, int64_t tick);

void write_s01_to_s09(const Config &c, const char *out_dir);
void write_s11_to_s38(Sim &s, const Config &c, const char *out_dir);
void write_s39_svgs(Sim &s, const char *out_dir);
void write_s40_manifest(const char *out_dir, const char *toolchain, const double *stage_s, int64_t ev_total);
