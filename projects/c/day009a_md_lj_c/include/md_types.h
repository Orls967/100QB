#ifndef MD_TYPES_H
#define MD_TYPES_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define MD_PI 3.14159265358979323846264338327950288
#define MD_TWO_PI 6.28318530717958647692528676655900577
#define MD_MASTER_SEED 20261127ULL
#define MD_WCA_RC 1.12246204830937298143353304967917952 /* 2^(1/6) */

#define MD_MAX_STAGES 64
#define MD_MAX_CRITERIA 32
#define MD_MAX_FILES 512
#define MD_MAX_FLAGS 128

typedef struct {
    double x;
    double y;
    double z;
} md_vec3;

typedef struct {
    double lx;
    double ly;
    double lz;
    double inv_lx;
    double inv_ly;
    double inv_lz;
    int dim; /* 2 or 3 */
} md_box;

typedef enum {
    POT_LJ_SHIFT = 0,
    POT_LJ_TAIL = 1,
    POT_WCA = 2,
    POT_MORSE = 3,
    POT_YUKAWA = 4,
    POT_COUNT = 5
} md_pot_kind;

typedef struct {
    md_pot_kind kind;
    double rc;
    double rc2;
    double eps;
    double sigma;
    /* Morse parameters: D_e, alpha, r0 */
    double morse_de;
    double morse_alpha;
    double morse_r0;
    /* Yukawa parameters: A, kappa */
    double yukawa_a;
    double yukawa_kappa;
    /* Shift value u(rc) */
    double u_shift;
    /* Tail correction prefactors (per rho and per rho^2) */
    bool use_tail;
    /* Optional lookup table */
    bool use_table;
    bool use_cubic_table;
    int tab_n;
    double tab_rmin2;
    double tab_rmax2;
    double tab_dr2;
    double tab_inv_dr2;
    double *tab_u;
    double *tab_du;   /* dU/d(r^2) */
    double *tab_fdivr; /* -(1/r) dU/dr = -2 dU/d(r^2) */
    double *tab_dfdivr;
} md_potential;

typedef enum {
    INT_VVERLET = 0,
    INT_LEAPFROG = 1,
    INT_EULER = 2,
    INT_RESPA = 3
} md_integrator_kind;

typedef enum {
    THERM_NONE = 0,       /* NVE */
    THERM_LANGEVIN = 1,   /* BAOAB zero-COM */
    THERM_NOSEHOOVER = 2, /* Chain-1 Nose-Hoover */
    THERM_BERENDSEN = 3,  /* Non-canonical rescaling */
    THERM_ANDERSEN = 4    /* Pairwise momentum-conserving */
} md_thermostat_kind;

typedef struct {
    md_thermostat_kind kind;
    double target_t;
    double gamma;       /* Langevin friction */
    double tau_t;       /* Berendsen relaxation time */
    double nu_coll;     /* Andersen collision frequency per particle */
    double nh_q;        /* Nose-Hoover thermal inertia */
    double nh_xi;       /* Nose-Hoover friction variable */
    double nh_eta;      /* Nose-Hoover integral of xi (for conserved H_NH) */
    /* Optional Berendsen barostat */
    bool use_barostat;
    double target_p;
    double tau_p;
    double beta_t;      /* Isothermal compressibility estimate */
} md_thermostat;

typedef enum {
    NBR_BRUTE = 0,
    NBR_CELL = 1,
    NBR_VERLET = 2,
    NBR_SOA_VERLET = 3,
    NBR_THREADED = 4,
    NBR_CELL_SORTED_SOA = 5
} md_nbr_kind;

typedef struct {
    int n;
    md_box box;
    md_vec3 *pos;       /* Wrapped positions in [0, L) */
    md_vec3 *unwrapped; /* Unwrapped continuous positions for MSD/Lindemann */
    md_vec3 *vel;
    md_vec3 *force;
    /* Instantaneous observables from last force call */
    double pe;
    double pe_tail;
    double virial;      /* sum_{i<j} r_ij . f_ij */
    double p_tail;
    /* Stress virial components (before dividing by V): sum_{i<j} r_a * f_b */
    double w_xx;
    double w_yy;
    double w_zz;
    double w_xy;
    double w_xz;
    double w_yz;
} md_system;

typedef struct {
    int n;
    double *x;
    double *y;
    double *z;
    double *ux;
    double *uy;
    double *uz;
    double *vx;
    double *vy;
    double *vz;
    double *fx;
    double *fy;
    double *fz;
    int *perm; /* Mapping from sorted slot to original particle id */
} md_soa;

typedef struct {
    uint64_t s[4];
    bool has_spare_normal;
    double spare_normal;
} md_rng;

typedef struct {
    double mean;
    double sd;
    double se;
    double min_val;
    double max_val;
    int n_blocks;
    int block_size;
} md_block_stats;

typedef struct {
    double slope;
    double intercept;
    double slope_se;
    double r2;
    int n_pts;
} md_ols_result;

typedef struct {
    char id[16];
    char name[64];
    double elapsed_sec;
    int files_written;
    bool completed;
} md_stage_record;

typedef struct {
    char path[256];
    char kind[16]; /* "csv", "svg", "json" */
    uint64_t fnv1a;
    size_t bytes;
    int rows;
    int cols;
} md_file_record;

typedef struct {
    char id[16];
    char status[32];    /* tercapai / sebagian / belum tercapai / hampa / belum diketahui */
    char measured[256];
    char threshold[192];
} md_criterion_record;

typedef struct {
    char stage[16];
    char file[128];
    char reason[256];
} md_vacuous_flag;

typedef struct {
    bool small_scale;
    int n_threads;
    int c_lines;
    char cc_ver[128];
    char arch[64];
    int cores;

    md_stage_record stages[MD_MAX_STAGES];
    int n_stages;

    md_file_record files[MD_MAX_FILES];
    int n_files;

    md_criterion_record criteria[MD_MAX_CRITERIA];
    int n_criteria;

    md_vacuous_flag flags[MD_MAX_FLAGS];
    int n_flags;

    /* Global diagnostics aggregated across stages */
    uint64_t total_nbr_overflows;
    uint64_t total_late_rebuilds;
    uint64_t combined_csv_hash;
    double total_elapsed_sec;

    /* Key measured values shared across stages for criteria evaluation */
    double c1_max_rel_err;
    double c2_max_rel_err;
    double c3_p_tot_1e5;
    double c4_rev_max_pos_err;
    double c5_vv_slope;
    double c5_vv_slope_se;
    double c5_euler_ratio;
    bool c6_all_within_3se;
    double c6_max_z;
    double c7_min_ks_p;
    bool c7_all_pass;
    bool c8_canon_within_10pct;
    bool c8_noncanon_smaller_3se;
    double c8_berendsen_z;
    double c8_nve_z;
    double c8_canon_max_rel;
    bool c9_all_within_3se;
    double c9_max_z;
    double c10_rel_err;
    double c10_z;
    double c11_diff_rel_err;
    bool c11_mono_3se;
    double c11_min_step_z;
    double c12_tail_dev;
    double c12_core_max;
    double c13_eta;
    double c13_eta_se;
    double c14_q6_solid;
    double c14_q6_liquid;
    double c14_q6_drop_z;
    double c14_cold_q6;
    bool c15_t15_mono_3se;
    double c15_t15_min_z;
    bool c15_t20_mono;
    bool c16_threads_bit_exact;
    bool c16_runs_bit_exact;
    uint64_t c16_hash_1t;
    uint64_t c16_hash_2t;
    uint64_t c16_hash_4t;
    uint64_t c16_hash_8t;
    double c17_brute_slope;
    double c17_brute_se;
    double c17_cell_slope;
    double c17_cell_se;
    double c17_verlet_slope;
    double c17_verlet_se;
    double c18_medians_us[6];
    double c18_mads_us[6];
    bool c18_monotone_faster;
    bool c20_selftest_ok;
} md_run_ctx;

#endif /* MD_TYPES_H */
