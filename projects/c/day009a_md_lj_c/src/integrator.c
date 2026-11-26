#include "md_integrator.h"
#include "md_init.h"
#include "md_observables.h"
#include "md_rng.h"
#include "md_vec.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void md_thermostat_init(md_thermostat *th, md_thermostat_kind kind, double target_t,
                        int n_particles, int dim) {
    memset(th, 0, sizeof(*th));
    th->kind = kind;
    th->target_t = target_t;
    th->gamma = 1.0;
    th->tau_t = 0.2;
    th->nu_coll = 2.0;
    int dof = dim * (n_particles > 1 ? (n_particles - 1) : 1);
    /* Nose-Hoover inertial mass Q = N_f * T * tau^2 with tau = 0.5 */
    th->nh_q = (double)dof * target_t * 0.25;
    if (th->nh_q < 1.0) th->nh_q = 1.0;
    th->nh_xi = 0.0;
    th->nh_eta = 0.0;
    th->use_barostat = false;
    th->target_p = 1.0;
    th->tau_p = 1.0;
    th->beta_t = 0.1;
}

void md_thermostat_enable_barostat(md_thermostat *th, double target_p, double tau_p, double beta_t) {
    th->use_barostat = true;
    th->target_p = target_p;
    th->tau_p = tau_p;
    th->beta_t = beta_t;
}

static void apply_langevin_o_step(md_system *sys, const md_thermostat *th, double dt, md_rng *rng) {
    /* Zero-COM Ornstein-Uhlenbeck step:
     * v_i <- c1 * v_i + c2 * (z_i - mean(z))
     * Orthogonal projection onto sum_i v_i = 0 preserves zero COM momentum strictly
     * and yields exact canonical stationary moments <K> = d(N-1)T/2, Var(K) = d(N-1)T^2/2! */
    int n = sys->n;
    double c1 = exp(-th->gamma * dt);
    double c2 = sqrt(th->target_t * (1.0 - c1 * c1));

    double *zx = (double *)malloc((size_t)(3 * n) * sizeof(double));
    double *zy = zx + n;
    double *zz = zx + 2 * n;
    double mx = 0.0, my = 0.0, mz = 0.0;

    for (int i = 0; i < n; ++i) {
        zx[i] = md_rng_normal(rng);
        zy[i] = md_rng_normal(rng);
        zz[i] = (sys->box.dim == 3) ? md_rng_normal(rng) : 0.0;
        mx += zx[i];
        my += zy[i];
        mz += zz[i];
    }
    double inv_n = 1.0 / (double)n;
    mx *= inv_n;
    my *= inv_n;
    mz *= inv_n;

    for (int i = 0; i < n; ++i) {
        sys->vel[i].x = c1 * sys->vel[i].x + c2 * (zx[i] - mx);
        sys->vel[i].y = c1 * sys->vel[i].y + c2 * (zy[i] - my);
        if (sys->box.dim == 3) {
            sys->vel[i].z = c1 * sys->vel[i].z + c2 * (zz[i] - mz);
        } else {
            sys->vel[i].z = 0.0;
        }
    }
    free(zx);
}

static void apply_nosehoover_half(md_system *sys, md_thermostat *th, double dt_half) {
    /* Martyna-Tuckerman-Tobias-Klein (MTTK) symmetric half-step for 1-chain Nose-Hoover:
     * 1) xi += (dt/4) * (2*K - N_f*T) / Q
     * 2) v  *= exp(-xi * dt/2); eta += xi * dt/2
     * 3) xi += (dt/4) * (2*K_new - N_f*T) / Q */
    int dof = md_dof(sys);
    double ke1 = md_kinetic_energy(sys);
    double g1 = (2.0 * ke1 - (double)dof * th->target_t) / th->nh_q;
    th->nh_xi += 0.5 * dt_half * g1;

    double sc = exp(-th->nh_xi * dt_half);
    for (int i = 0; i < sys->n; ++i) {
        sys->vel[i].x *= sc;
        sys->vel[i].y *= sc;
        sys->vel[i].z *= sc;
    }
    th->nh_eta += th->nh_xi * dt_half;

    double ke2 = md_kinetic_energy(sys);
    double g2 = (2.0 * ke2 - (double)dof * th->target_t) / th->nh_q;
    th->nh_xi += 0.5 * dt_half * g2;
}

static void apply_berendsen_step(md_system *sys, const md_thermostat *th, double dt) {
    double t_inst = md_temperature(sys);
    if (t_inst <= 1e-12) return;
    double ratio = 1.0 + (dt / th->tau_t) * (th->target_t / t_inst - 1.0);
    if (ratio < 0.25) ratio = 0.25;
    if (ratio > 4.0) ratio = 4.0;
    double lam = sqrt(ratio);
    for (int i = 0; i < sys->n; ++i) {
        sys->vel[i].x *= lam;
        sys->vel[i].y *= lam;
        sys->vel[i].z *= lam;
    }
}

static void apply_andersen_step(md_system *sys, const md_thermostat *th, double dt, md_rng *rng) {
    /* Pairwise momentum-conserving Andersen collisions:
     * Expected pair collisions per step = 0.5 * N * nu * dt.
     * For each chosen pair (i, j), preserve v_+ = (v_i + v_j)/2 and draw
     * relative velocity v_- ~ N(0, 2*T*I), which satisfies exact canonical detailed
     * balance on the zero-momentum manifold sum_i v_i = 0! */
    int n = sys->n;
    if (n < 2) return;
    double p_coll = th->nu_coll * dt;
    double s2t = sqrt(2.0 * th->target_t);
    int n_trials = n / 2;
    for (int k = 0; k < n_trials; ++k) {
        if (md_rng_uniform(rng) < p_coll) {
            int i = (int)(md_rng_u64(rng) % (uint64_t)n);
            int j = (int)(md_rng_u64(rng) % (uint64_t)(n - 1));
            if (j >= i) j++;

            double vpx = 0.5 * (sys->vel[i].x + sys->vel[j].x);
            double vpy = 0.5 * (sys->vel[i].y + sys->vel[j].y);
            double vpz = 0.5 * (sys->vel[i].z + sys->vel[j].z);

            double vmx = s2t * md_rng_normal(rng);
            double vmy = s2t * md_rng_normal(rng);
            double vmz = (sys->box.dim == 3) ? (s2t * md_rng_normal(rng)) : 0.0;

            sys->vel[i].x = vpx + 0.5 * vmx;
            sys->vel[i].y = vpy + 0.5 * vmy;
            sys->vel[i].z = (sys->box.dim == 3) ? (vpz + 0.5 * vmz) : 0.0;

            sys->vel[j].x = vpx - 0.5 * vmx;
            sys->vel[j].y = vpy - 0.5 * vmy;
            sys->vel[j].z = (sys->box.dim == 3) ? (vpz - 0.5 * vmz) : 0.0;
        }
    }
}

static void apply_barostat_berendsen(md_system *sys, const md_thermostat *th,
                                     md_verlet_list *vl, double dt) {
    double p_inst = md_pressure(sys);
    double factor = 1.0 - (th->beta_t * dt / th->tau_p) * (th->target_p - p_inst);
    if (factor < 0.95) factor = 0.95;
    if (factor > 1.05) factor = 1.05;
    double mu = (sys->box.dim == 3) ? cbrt(factor) : sqrt(factor);
    md_system_rescale_box(sys, mu);
    if (vl) {
        /* Rescaling alters all interparticle distances; invalidate Verlet list cleanly */
        vl->built = false;
    }
}

void md_step_vverlet(md_system *sys, const md_potential *pot, md_nbr_kind method,
                     md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                     double dt, md_thermostat *th, md_rng *rng) {
    int n = sys->n;
    double dt_half = 0.5 * dt;
    md_thermostat_kind tk = th ? th->kind : THERM_NONE;

    if (tk == THERM_NOSEHOOVER) {
        apply_nosehoover_half(sys, th, dt_half);
    }

    /* B: half velocity kick */
    for (int i = 0; i < n; ++i) {
        sys->vel[i].x += dt_half * sys->force[i].x;
        sys->vel[i].y += dt_half * sys->force[i].y;
        sys->vel[i].z += dt_half * sys->force[i].z;
    }

    if (tk == THERM_LANGEVIN) {
        /* BAOAB splitting: A(dt/2) -> O(dt) -> A(dt/2) */
        for (int i = 0; i < n; ++i) {
            sys->unwrapped[i].x += dt_half * sys->vel[i].x;
            sys->unwrapped[i].y += dt_half * sys->vel[i].y;
            sys->unwrapped[i].z += dt_half * sys->vel[i].z;
        }
        apply_langevin_o_step(sys, th, dt, rng);
        for (int i = 0; i < n; ++i) {
            sys->unwrapped[i].x += dt_half * sys->vel[i].x;
            sys->unwrapped[i].y += dt_half * sys->vel[i].y;
            sys->unwrapped[i].z += dt_half * sys->vel[i].z;
            sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);
        }
    } else {
        /* Standard A(dt) drift */
        for (int i = 0; i < n; ++i) {
            sys->unwrapped[i].x += dt * sys->vel[i].x;
            sys->unwrapped[i].y += dt * sys->vel[i].y;
            sys->unwrapped[i].z += dt * sys->vel[i].z;
            sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);
        }
    }

    /* Recompute forces at new positions */
    md_forces_compute(sys, pot, method, cl, vl, soa, n_threads);

    /* Second B: half velocity kick */
    for (int i = 0; i < n; ++i) {
        sys->vel[i].x += dt_half * sys->force[i].x;
        sys->vel[i].y += dt_half * sys->force[i].y;
        sys->vel[i].z += dt_half * sys->force[i].z;
    }

    if (tk == THERM_NOSEHOOVER) {
        apply_nosehoover_half(sys, th, dt_half);
    } else if (tk == THERM_BERENDSEN) {
        apply_berendsen_step(sys, th, dt);
    } else if (tk == THERM_ANDERSEN) {
        apply_andersen_step(sys, th, dt, rng);
    }

    if (th && th->use_barostat) {
        apply_barostat_berendsen(sys, th, vl, dt);
        md_forces_compute(sys, pot, method, cl, vl, soa, n_threads);
    }
}

void md_step_leapfrog(md_system *sys, const md_potential *pot, md_nbr_kind method,
                      md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                      double dt) {
    /* Synchronized Leapfrog: v(t + dt/2) = v(t) + (dt/2)a(t), x(t+dt) = x(t) + dt*v(t+dt/2),
     * compute a(t+dt), reconstruct synchronized v(t+dt) = v(t+dt/2) + (dt/2)a(t+dt).
     * For standalone leapfrog energy diagnostics we store half-step kinetic via explicit staggered
     * update using v_{n+1/2} = v_{n-1/2} + dt * a_n when caller kicks -dt/2 initially,
     * OR here we perform canonical drift-kick-drift leapfrog (DKD) which is mathematically
     * dual to Velocity Verlet (KDK) and gives distinct phase-space shadowing! */
    int n = sys->n;
    double dt_half = 0.5 * dt;

    /* D(dt/2) */
    for (int i = 0; i < n; ++i) {
        sys->unwrapped[i].x += dt_half * sys->vel[i].x;
        sys->unwrapped[i].y += dt_half * sys->vel[i].y;
        sys->unwrapped[i].z += dt_half * sys->vel[i].z;
        sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);
    }
    /* Force at midpoint x(t + dt/2) */
    md_forces_compute(sys, pot, method, cl, vl, soa, n_threads);
    /* K(dt) */
    for (int i = 0; i < n; ++i) {
        sys->vel[i].x += dt * sys->force[i].x;
        sys->vel[i].y += dt * sys->force[i].y;
        sys->vel[i].z += dt * sys->force[i].z;
    }
    /* D(dt/2) */
    for (int i = 0; i < n; ++i) {
        sys->unwrapped[i].x += dt_half * sys->vel[i].x;
        sys->unwrapped[i].y += dt_half * sys->vel[i].y;
        sys->unwrapped[i].z += dt_half * sys->vel[i].z;
        sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);
    }
    /* Update forces at full step for accurate potential energy readout */
    md_forces_compute(sys, pot, method, cl, vl, soa, n_threads);
}

void md_step_euler(md_system *sys, const md_potential *pot, md_nbr_kind method,
                   md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads,
                   double dt) {
    /* Explicit forward Euler (non-symplectic 1st-order negative control) */
    int n = sys->n;
    for (int i = 0; i < n; ++i) {
        sys->unwrapped[i].x += dt * sys->vel[i].x;
        sys->unwrapped[i].y += dt * sys->vel[i].y;
        sys->unwrapped[i].z += dt * sys->vel[i].z;
        sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);

        sys->vel[i].x += dt * sys->force[i].x;
        sys->vel[i].y += dt * sys->force[i].y;
        sys->vel[i].z += dt * sys->force[i].z;
    }
    md_forces_compute(sys, pot, method, cl, vl, soa, n_threads);
}

void md_step_respa(md_system *sys, const md_potential *pot, md_verlet_list *vl,
                   double dt_outer, int n_inner) {
    int n = sys->n;
    double dt_in = dt_outer / (double)n_inner;
    double half_out = 0.5 * dt_outer;
    double half_in = 0.5 * dt_in;

    /* Outer slow force half kick */
    md_forces_respa_slow(sys, pot, vl);
    for (int i = 0; i < n; ++i) {
        sys->vel[i].x += half_out * sys->force[i].x;
        sys->vel[i].y += half_out * sys->force[i].y;
        sys->vel[i].z += half_out * sys->force[i].z;
    }

    /* Inner fast WCA velocity Verlet loop */
    md_forces_respa_fast(sys, pot, vl);
    for (int s = 0; s < n_inner; ++s) {
        for (int i = 0; i < n; ++i) {
            sys->vel[i].x += half_in * sys->force[i].x;
            sys->vel[i].y += half_in * sys->force[i].y;
            sys->vel[i].z += half_in * sys->force[i].z;

            sys->unwrapped[i].x += dt_in * sys->vel[i].x;
            sys->unwrapped[i].y += dt_in * sys->vel[i].y;
            sys->unwrapped[i].z += dt_in * sys->vel[i].z;
            sys->pos[i] = md_box_wrap(&sys->box, sys->unwrapped[i]);
        }
        md_forces_respa_fast(sys, pot, vl);
        for (int i = 0; i < n; ++i) {
            sys->vel[i].x += half_in * sys->force[i].x;
            sys->vel[i].y += half_in * sys->force[i].y;
            sys->vel[i].z += half_in * sys->force[i].z;
        }
    }

    /* Second outer slow force half kick */
    double pe_fast = sys->pe;
    double vir_fast = sys->virial;
    md_forces_respa_slow(sys, pot, vl);
    for (int i = 0; i < n; ++i) {
        sys->vel[i].x += half_out * sys->force[i].x;
        sys->vel[i].y += half_out * sys->force[i].y;
        sys->vel[i].z += half_out * sys->force[i].z;
    }
    /* Combined PE and virial */
    sys->pe += pe_fast;
    sys->virial += vir_fast;
}

double md_nosehoover_conserved(const md_system *sys, const md_thermostat *th) {
    int dof = md_dof(sys);
    double ke = md_kinetic_energy(sys);
    return ke + sys->pe + 0.5 * th->nh_q * th->nh_xi * th->nh_xi +
           (double)dof * th->target_t * th->nh_eta;
}
