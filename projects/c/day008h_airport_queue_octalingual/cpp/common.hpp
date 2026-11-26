#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

constexpr int64_t BP = 10000;
constexpr int64_t SLOTS_PER_DAY = 8640;
constexpr int64_t SLOTS_PER_HOUR = 360;
constexpr int N_STN = 6;
constexpr int MAX_SRV = 16;
constexpr int MAX_Q = 4096;
constexpr int MAX_POOL = 8192;
constexpr int AUDIT_DAYS = 3;

enum { STN_CHECKIN=0, STN_SECURITY=1, STN_IMIGRASI=2, STN_GATE=3, STN_BELT=4, STN_CUSTOMS=5 };
enum { CLS_ECON=0, CLS_BIZ=1, CLS_TRANSIT=2 };
enum { DEST_DOM=0, DEST_INTL=1 };
enum { SRV_IDLE=0, SRV_BUSY=1, SRV_DOWN=2 };
enum {
    BR_ARR_ACC=0, BR_ARR_REJ=1, BR_SVC_DONE=2, BR_SEC_EXTRA=3, BR_SRV_BREAK=4,
    BR_SRV_REPAIR=5, BR_ABANDON=6, BR_MISSED_FLT=7, BR_REBOOK=8, BR_BAG_MISS=9,
    BR_BELT_BLOCK=10, BR_ROSTER_DROP=11, BR_ROUTE_INTL=12, BR_ROUTE_TRANSIT=13, BR_BIZ_PRIO=14
};

struct Xoshiro256 { uint64_t s[4]{}; };

struct Config {
    int64_t master_seed{}, slots_per_day{}, slots_per_hour{}, parity_days{}, arr_scale_bp{};
    int64_t srv_cnt[N_STN]{}, mu_bp[N_STN]{}, security_extra_prob_bp{}, security_extra_mu_bp{};
    int64_t belt_cap{}, break_prob_bp{}, repair_prob_bp{}, patience_slots{}, gate_close_lead{};
    int64_t bag_miss_bp{}, intl_prob_bp{}, biz_prob_bp{}, transit_prob_bp{}, rebook_prob_bp{};
};

struct RunParams {
    int64_t days{}, arr_scale_bp{}, break_prob_bp{}, repair_prob_bp{}, belt_cap{}, gate_close_lead{};
    int64_t srv_ovr[N_STN]{-1,-1,-1,-1,-1,-1}, mu_ovr[N_STN]{-1,-1,-1,-1,-1,-1};
    int32_t priority_on{1}, breakdown_on{1}, roster_on{1}, customs_on{0}, permute_srv{0}, timescale_div2{0}, audit_on{0};
};

struct Agent {
    int32_t active{}, cls{}, dest{}, flight_hr{}, n_bags{}, cur_stn{}, extra_sec{}, is_bag{};
    int64_t arr_slot{}, q_enter_slot{}, dep_slot{};
};

struct Server { int32_t status{}, enabled{}, agent_idx{-1}, extra_sec{}; };

struct Station {
    int32_t q[MAX_Q]{}, head{}, tail{}, len{}, cap{}, has_prio{};
    int32_t c_max{}, c_active{}, busy_cnt{};
    Server srv[MAX_SRV]{};
    int64_t arr_cnt{}, dep_cnt{}, block_cnt{}, aban_cnt{}, sum_L{}, sum_W{}, wait_cnt{}, busy_slots{};
};

struct StnEvent { int32_t slot, stn, delta; };

struct Sim {
    Config cfg{}; RunParams rp{};
    Xoshiro256 rng[8]{};
    Agent pool[MAX_POOL]{};
    int32_t free_stk[MAX_POOL]{}, free_top{};
    Station stn[N_STN]{};
    int32_t arr_prof[24]{};
    int32_t roster[24][N_STN]{};
    int64_t branches[15]{};
    int64_t total_arr{}, total_done{}, total_lost{}, total_missed{}, total_rebook{}, total_aban{};
    int64_t total_bags{}, total_bags_miss{}, total_belt_block{}, total_breakdowns{};
    int64_t sys_sum_L{}, sys_sum_sojourn{}, inv_violations{};
    uint64_t event_hash{0xcbf29ce484222325ULL};
    int64_t cls_wait_sum[3]{}, cls_wait_cnt[3]{}, cls_sojourn_sum[3]{}, cls_sojourn_cnt[3]{};
    int64_t stn_wait_sum[N_STN]{}, stn_wait_cnt[N_STN]{};
    int64_t hr_arr[24]{}, hr_util_sum[24][N_STN]{};
    std::vector<int32_t> wait_samples, sojourn_samples, q_series, srv_series;
    std::vector<StnEvent> stn_ev[N_STN];
    int32_t audit_slots{};
};

uint64_t splitmix64(uint64_t *x);
uint64_t fnv1a_bytes(const void *buf, size_t n);
uint64_t fnv1a_u32pair(uint32_t a, uint32_t b);
uint64_t fnv1a_i64(uint64_t h, int64_t v);
Xoshiro256 xoshiro_seed(uint64_t sd);
uint64_t xoshiro_next(Xoshiro256 *st);
int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi);
int32_t bernoulli(Xoshiro256 *st, int64_t p_bp);

Config load_config(const char *path);
RunParams default_rp(const Config &c, int64_t days);

struct DtmcRes { double p_block{}, L{}, L_q{}, U{}, W{}, W_q{}, lam_eff{}; };
DtmcRes solve_dtmc_single(int c, int K, double lam, double mu);
void solve_dtmc_tandem(int K1, int K2, double lam, double mu1, double mu2, DtmcRes &r1, DtmcRes &r2, double &p_blk_tnd);
DtmcRes solve_geo_closed_form(double lam, double mu);

void stn_init(Station &st, int c_max, int cap, int has_prio);
void stn_push_tail(Station &st, int32_t idx);
void stn_push_head(Station &st, int32_t idx);
int32_t stn_pop_next(Station &st, const Agent *pool, int32_t priority_on, int64_t *prio_br);
int32_t stn_remove_at(Station &st, int pos);
int32_t pool_alloc(Sim &s);
void pool_free(Sim &s, int32_t idx);
std::unique_ptr<Sim> sim_create(const Config &c, const RunParams &rp);
void sim_run(Sim &s);
void sim_step_slot(Sim &s, int64_t slot);

void write_s01_to_s13(const Config &c, const char *out_dir);
void write_s14_to_s45(const Config &c, const char *out_dir, Sim &s_main);
void write_s46_svgs(const Config &c, const char *out_dir, Sim &s_main);
