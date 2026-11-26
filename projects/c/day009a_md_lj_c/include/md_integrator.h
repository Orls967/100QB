#ifndef MD_INTEGRATOR_H
#define MD_INTEGRATOR_H

#include "md_forces.h"
#include "md_types.h"

/* Initialize thermostat presets */
void md_thermostat_init(md_thermostat *th, md_thermostat_kind kind, double target_t,
                        int n_particles, int dim);
void md_thermostat_enable_barostat(md_thermostat *th, double target_p, double tau_p, double beta_t);

/* Single MD integration step supporting NVE and all thermostats/barostats */
void md_step_vverlet(md_system *sys, const md_potential *pot, md_nbr_kind method,
                     md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                     double dt, md_thermostat *th, md_rng *rng);

void md_step_leapfrog(md_system *sys, const md_potential *pot, md_nbr_kind method,
                      md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                      double dt);

void md_step_euler(md_system *sys, const md_potential *pot, md_nbr_kind method,
                   md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                   double dt);

/* Two-timescale r-RESPA step: n_inner fast WCA steps per outer attractive step */
void md_step_respa(md_system *sys, const md_potential *pot, md_verlet_list *vl,
                   double dt_outer, int n_inner);

/* Extended Nose-Hoover conserved pseudo-Hamiltonian H_NH = K + U + 0.5*Q*xi^2 + N_f*T*eta */
double md_nosehoover_conserved(const md_system *sys, const md_thermostat *th);

#endif /* MD_INTEGRATOR_H */
