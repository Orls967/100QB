#ifndef MD_FORCES_H
#define MD_FORCES_H

#include "md_neighbor.h"
#include "md_potential.h"
#include "md_types.h"

/* V0: Brute-force O(N^2) AoS with Newton's 3rd law */
void md_forces_brute(md_system *sys, const md_potential *pot);

/* V1: Linked-cell list O(N) AoS with Newton's 3rd law */
void md_forces_cell(md_system *sys, const md_potential *pot, md_cell_list *cl);

/* V2: Verlet half-list AoS with Newton's 3rd law */
void md_forces_verlet(md_system *sys, const md_potential *pot, md_verlet_list *vl);

/* V3: Verlet half-list on SoA layout */
void md_forces_soa_verlet(md_system *sys, md_soa *soa, const md_potential *pot, md_verlet_list *vl);

/* V4: Deterministic pthread parallel force evaluation using sorted full neighbor list.
 * Guaranteed bitwise identical across 1, 2, 4, 8 threads! */
void md_forces_threaded(md_system *sys, const md_potential *pot, md_verlet_list *vl, int n_threads);

/* V5: Cell-sorted SoA layout with Verlet list (supports 1 or multi-threaded) */
void md_forces_cell_sorted_soa(md_system *sys, md_soa *soa, const md_potential *pot,
                               md_verlet_list *vl, int n_threads);

/* Unified dispatcher by neighbor method */
void md_forces_compute(md_system *sys, const md_potential *pot, md_nbr_kind method,
                       md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads);

/* r-RESPA split forces: fast short-range WCA and slow attractive LJ tail */
void md_forces_respa_fast(md_system *sys, const md_potential *pot, md_verlet_list *vl);
void md_forces_respa_slow(md_system *sys, const md_potential *pot, md_verlet_list *vl);

#endif /* MD_FORCES_H */
