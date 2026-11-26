#ifndef MD_NEIGHBOR_H
#define MD_NEIGHBOR_H

#include "md_types.h"

/* 3D / 2D linked-cell list */
typedef struct {
    int nx, ny, nz;
    int n_cells;
    int cap_n;
    int *head; /* [n_cells] */
    int *next; /* [n] */
    int *cell_of; /* [n] */
} md_cell_list;

/* Verlet neighbor list (supports half-list for serial Newton-3rd and full sorted list for deterministic threads) */
typedef struct {
    int n;
    double r_cut;
    double r_skin;
    double r_list;
    double r_list2;

    /* Half neighbor list: only j > i */
    int *half_offsets; /* [n + 1] */
    int *half_nbrs;
    int half_cap;
    int half_count;

    /* Full neighbor list: all j != i in strictly ascending order of j */
    int *full_offsets; /* [n + 1] */
    int *full_nbrs;
    int full_cap;
    int full_count;

    /* Reference positions at last rebuild for displacement trigger */
    md_vec3 *ref_pos;
    bool built;
    bool need_full;

    /* Diagnostics & invariants (c23) */
    uint64_t rebuild_count;
    uint64_t check_count;
    uint64_t overflow_count;
    uint64_t late_rebuild_count;
    double max_disp_before_rebuild;
    double last_max_disp;

    /* Internal cell list used to accelerate Verlet rebuild in O(N) */
    md_cell_list cells;
} md_verlet_list;

void md_cell_init(md_cell_list *cl, int n);
void md_cell_free(md_cell_list *cl);
void md_cell_build(md_cell_list *cl, const md_system *sys, double r_cut);

void md_verlet_init(md_verlet_list *vl, int n, double r_cut, double r_skin, bool need_full);
void md_verlet_free(md_verlet_list *vl);
void md_verlet_build(md_verlet_list *vl, const md_system *sys);

/* Returns true if rebuild was triggered; audits displacement invariant <= r_skin/2 */
bool md_verlet_check_and_rebuild(md_verlet_list *vl, const md_system *sys);

/* Reorder SoA arrays by cell index for V5 cache locality */
void md_soa_sort_by_cell(md_soa *soa, const md_system *sys, double r_cut);

#endif /* MD_NEIGHBOR_H */
