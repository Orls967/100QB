#include "md_forces.h"
#include "md_init.h"
#include "md_neighbor.h"
#include "md_observables.h"
#include "md_potential.h"
#include "md_rng.h"
#include "md_stages.h"
#include "md_vec.h"
#include <math.h>
#include <stdio.h>

int md_run_selftest(void) {
    int failures = 0;
    printf("[selftest] Starting built-in self-test suite...\n");

    /* 1. Two-particle analytical force & energy at exact distances */
    {
        md_potential pot;
        md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
        double test_r[4] = {1.0, MD_WCA_RC, 1.5, 2.0};
        for (int k = 0; k < 4; ++k) {
            double r = test_r[k];
            md_system sys;
            md_init_fixture_dimer(&sys, r, 10.0);
            md_forces_brute(&sys, &pot);

            double inv_r2 = 1.0 / (r * r);
            double inv_r6 = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;
            double rc_inv6 = pow(1.0 / 2.5, 6.0);
            double u_shift = 4.0 * (rc_inv6 * rc_inv6 - rc_inv6);
            double u_exact = 4.0 * (inv_r12 - inv_r6) - u_shift;
            /* Particle 0 is at x = L/2 - r/2, particle 1 at x = L/2 + r/2, so dx = -r.
             * Repulsive force pushes particle 0 in -x direction: fx0 = -24*(2*r^-12 - r^-6)/r */
            double f_mag = 24.0 * (2.0 * inv_r12 - inv_r6) / r;
            double fx0_exact = -f_mag;

            double err_u = fabs(sys.pe - u_exact);
            double err_f = fabs(sys.force[0].x - fx0_exact) + fabs(sys.force[1].x + fx0_exact);
            if (err_u > 1e-12 || err_f > 1e-12) {
                printf("  [FAIL] dimer r=%.6f err_u=%.3e err_f=%.3e\n", r, err_u, err_f);
                failures++;
            }
            md_system_free(&sys);
        }
        printf("  [OK] 1. Analytical dimer force & energy (4 distances, LJ-shift)\n");
    }

    /* 2. Small FCC lattice (N=108, nc=3, rho=1.0) static energy vs exact shell formula */
    {
        md_system sys;
        md_init_fcc(&sys, 3, 1.0);
        md_potential pot;
        md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
        md_forces_brute(&sys, &pot);

        /* For FCC with lattice parameter a = (4/rho)^(1/3):
         * Shell distances: sqrt(1/2)*a (12), a (6), sqrt(3/2)*a (24), sqrt(2)*a (12), sqrt(5/2)*a (24) */
        double a = cbrt(4.0 / 1.0);
        double shell_r[5] = {
            a * sqrt(0.5), a * 1.0, a * sqrt(1.5), a * sqrt(2.0), a * sqrt(2.5)
        };
        int shell_cnt[5] = {12, 6, 24, 12, 24};
        double u_formula_per_n = 0.0;
        for (int s = 0; s < 5; ++s) {
            if (shell_r[s] < pot.rc) {
                double fdr = 0.0;
                double u_pair = md_pot_eval_exact(&pot, shell_r[s] * shell_r[s], &fdr);
                u_formula_per_n += 0.5 * (double)shell_cnt[s] * u_pair;
            }
        }
        double u_meas_per_n = sys.pe / (double)sys.n;
        double err_fcc = fabs(u_meas_per_n - u_formula_per_n);
        double max_f = 0.0;
        for (int i = 0; i < sys.n; ++i) {
            double fn = v3_norm(sys.force[i]);
            if (fn > max_f) max_f = fn;
        }
        double q6 = md_compute_q6(&sys, 1.35, NULL);
        if (err_fcc > 1e-11 || max_f > 1e-11 || fabs(q6 - 0.57452426) > 1e-5) {
            printf("  [FAIL] FCC static check: err_u=%.3e max_f=%.3e q6=%.8f\n",
                   err_fcc, max_f, q6);
            failures++;
        } else {
            printf("  [OK] 2. FCC N=108 lattice energy (u=%.6f, err=%.2e, Q6=%.6f)\n",
                   u_meas_per_n, err_fcc, q6);
        }
        md_system_free(&sys);
    }

    /* 3. Zero COM momentum & exact kinetic temperature initialization */
    {
        md_system sys;
        md_init_fcc(&sys, 4, 0.8);
        md_rng rng;
        md_rng_seed(&rng, MD_MASTER_SEED);
        md_init_velocities(&sys, 1.5, true, &rng);
        md_vec3 p_tot = md_total_momentum(&sys);
        double p_norm = v3_norm(p_tot);
        double t_meas = md_temperature(&sys);
        if (p_norm > 1e-12 || fabs(t_meas - 1.5) > 1e-12) {
            printf("  [FAIL] zero momentum check: |P|=%.3e T=%.12f\n", p_norm, t_meas);
            failures++;
        } else {
            printf("  [OK] 3. Zero COM momentum (|P|=%.2e, T=%.6f)\n", p_norm, t_meas);
        }
        md_system_free(&sys);
    }

    /* 4. Potential lookup table (cubic Hermite) vs analytical across all 5 potentials */
    {
        double worst_u = 0.0, worst_f = 0.0;
        for (int k = 0; k < POT_COUNT; ++k) {
            md_potential pot;
            md_potential_init(&pot, (md_pot_kind)k, 2.5);
            md_pot_build_table(&pot, 4096, 0.65, true);
            for (int i = 0; i < 100; ++i) {
                double r = 0.85 + (pot.rc - 0.90) * ((double)i + 0.37) / 100.0;
                double r2 = r * r;
                double fe = 0.0, ft = 0.0;
                double ue = md_pot_eval_exact(&pot, r2, &fe);
                double ut = md_pot_eval(&pot, r2, &ft);
                double du = fabs(ut - ue) / (fabs(ue) + 1.0);
                double df = fabs(ft - fe) / (fabs(fe) + 1.0);
                if (du > worst_u) worst_u = du;
                if (df > worst_f) worst_f = df;
            }
            md_potential_free(&pot);
        }
        if (worst_u > 1e-7 || worst_f > 1e-6) {
            printf("  [FAIL] potential table error: worst_u=%.3e worst_f=%.3e\n", worst_u, worst_f);
            failures++;
        } else {
            printf("  [OK] 4. Potential lookup tables (5 potentials, max_du=%.2e, max_df=%.2e)\n",
                   worst_u, worst_f);
        }
    }

    printf("[selftest] Finished with %d failure(s).\n", failures);
    return failures == 0 ? 0 : 1;
}
