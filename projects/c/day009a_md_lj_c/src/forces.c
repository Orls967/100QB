#include "md_forces.h"
#include "md_vec.h"
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double pe, vir, wxx, wyy, wzz, wxy, wxz, wyz;
} force_acc;

static void zero_system_forces(md_system *sys) {
    memset(sys->force, 0, (size_t)sys->n * sizeof(md_vec3));
    sys->pe = sys->pe_tail = sys->virial = sys->p_tail = 0.0;
    sys->w_xx = sys->w_yy = sys->w_zz = sys->w_xy = sys->w_xz = sys->w_yz = 0.0;
}

static void finish_system_forces(md_system *sys, const md_potential *pot, const force_acc *a) {
    sys->pe = a->pe;
    sys->virial = a->vir;
    sys->w_xx = a->wxx; sys->w_yy = a->wyy; sys->w_zz = a->wzz;
    sys->w_xy = a->wxy; sys->w_xz = a->wxz; sys->w_yz = a->wyz;
    sys->pe_tail = 0.0;
    sys->p_tail = 0.0;
    if (pot->use_tail && sys->box.dim == 3) {
        double ut = 0.0, pt = 0.0;
        md_pot_tail_corrections(pot, sys->n, md_box_volume(&sys->box), &ut, &pt);
        sys->pe_tail = ut;
        sys->p_tail = pt;
        sys->pe += ut;
    }
}

static inline void accum_pair_n3(md_system *sys, int i, int j, md_vec3 d, double r2,
                                 const md_potential *pot, force_acc *a) {
    double fdr = 0.0;
    double u = md_pot_eval(pot, r2, &fdr);
    double fx = fdr * d.x, fy = fdr * d.y, fz = fdr * d.z;
    sys->force[i].x += fx; sys->force[i].y += fy; sys->force[i].z += fz;
    sys->force[j].x -= fx; sys->force[j].y -= fy; sys->force[j].z -= fz;
    a->pe += u;
    a->vir += d.x * fx + d.y * fy + d.z * fz;
    a->wxx += d.x * fx; a->wyy += d.y * fy; a->wzz += d.z * fz;
    a->wxy += d.x * fy; a->wxz += d.x * fz; a->wyz += d.y * fz;
}

void md_forces_brute(md_system *sys, const md_potential *pot) {
    zero_system_forces(sys);
    force_acc a = {0};
    double rc2 = pot->rc2;
    for (int i = 0; i < sys->n; ++i) {
        md_vec3 pi = sys->pos[i];
        for (int j = i + 1; j < sys->n; ++j) {
            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
            double r2 = v3_norm2(d);
            if (r2 < rc2) accum_pair_n3(sys, i, j, d, r2, pot, &a);
        }
    }
    finish_system_forces(sys, pot, &a);
}

void md_forces_cell(md_system *sys, const md_potential *pot, md_cell_list *cl) {
    md_cell_build(cl, sys, pot->rc);
    int nx = cl->nx, ny = cl->ny, nz = cl->nz;
    if (nx < 3 || ny < 3 || (sys->box.dim == 3 && nz < 3)) {
        md_forces_brute(sys, pot);
        return;
    }
    zero_system_forces(sys);
    force_acc a = {0};
    double rc2 = pot->rc2;
    static const int off3[13][3] = {
        {1,0,0},{1,1,0},{0,1,0},{-1,1,0},{0,0,1},{1,0,1},{1,1,1},
        {0,1,1},{-1,1,1},{-1,0,1},{-1,-1,1},{0,-1,1},{1,-1,1}
    };
    static const int off2[4][3] = {{1,0,0},{1,1,0},{0,1,0},{-1,1,0}};
    int n_off = (sys->box.dim == 3) ? 13 : 4;
    const int (*offs)[3] = (sys->box.dim == 3) ? off3 : off2;

    for (int cz = 0; cz < nz; ++cz) {
        for (int cy = 0; cy < ny; ++cy) {
            for (int cx = 0; cx < nx; ++cx) {
                int cid = (cz * ny + cy) * nx + cx;
                for (int i = cl->head[cid]; i != -1; i = cl->next[i]) {
                    md_vec3 pi = sys->pos[i];
                    for (int j = cl->next[i]; j != -1; j = cl->next[j]) {
                        md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
                        double r2 = v3_norm2(d);
                        if (r2 < rc2) accum_pair_n3(sys, i, j, d, r2, pot, &a);
                    }
                }
                for (int o = 0; o < n_off; ++o) {
                    int ncx = (cx + offs[o][0] + nx) % nx;
                    int ncy = (cy + offs[o][1] + ny) % ny;
                    int ncz = (cz + offs[o][2] + nz) % nz;
                    int ncid = (ncz * ny + ncy) * nx + ncx;
                    for (int i = cl->head[cid]; i != -1; i = cl->next[i]) {
                        md_vec3 pi = sys->pos[i];
                        for (int j = cl->head[ncid]; j != -1; j = cl->next[j]) {
                            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
                            double r2 = v3_norm2(d);
                            if (r2 < rc2) accum_pair_n3(sys, i, j, d, r2, pot, &a);
                        }
                    }
                }
            }
        }
    }
    finish_system_forces(sys, pot, &a);
}

void md_forces_verlet(md_system *sys, const md_potential *pot, md_verlet_list *vl) {
    (void)md_verlet_check_and_rebuild(vl, sys);
    zero_system_forces(sys);
    force_acc a = {0};
    double rc2 = pot->rc2;
    for (int i = 0; i < sys->n; ++i) {
        md_vec3 pi = sys->pos[i];
        int start = vl->half_offsets[i], end = vl->half_offsets[i + 1];
        for (int k = start; k < end; ++k) {
            int j = vl->half_nbrs[k];
            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
            double r2 = v3_norm2(d);
            if (r2 < rc2) accum_pair_n3(sys, i, j, d, r2, pot, &a);
        }
    }
    finish_system_forces(sys, pot, &a);
}

void md_forces_soa_verlet(md_system *sys, md_soa *soa, const md_potential *pot, md_verlet_list *vl) {
    (void)md_verlet_check_and_rebuild(vl, sys);
    md_soa_from_system(soa, sys);
    force_acc a = {0};
    double rc2 = pot->rc2;
    double lx = sys->box.lx, ly = sys->box.ly, lz = sys->box.lz;
    double ilx = sys->box.inv_lx, ily = sys->box.inv_ly, ilz = sys->box.inv_lz;
    int dim = sys->box.dim;
    const double *restrict px = soa->x, *restrict py = soa->y, *restrict pz = soa->z;
    double *restrict fx_a = soa->fx, *restrict fy_a = soa->fy, *restrict fz_a = soa->fz;

    for (int i = 0; i < sys->n; ++i) {
        double xi = px[i], yi = py[i], zi = pz[i];
        double fxi = fx_a[i], fyi = fy_a[i], fzi = fz_a[i];
        int start = vl->half_offsets[i], end = vl->half_offsets[i + 1];
        for (int k = start; k < end; ++k) {
            int j = vl->half_nbrs[k];
            double dx = md_mic_1d(xi - px[j], lx, ilx);
            double dy = md_mic_1d(yi - py[j], ly, ily);
            double dz = (dim == 3) ? md_mic_1d(zi - pz[j], lz, ilz) : 0.0;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                double fdr = 0.0;
                double u = md_pot_eval(pot, r2, &fdr);
                double fx = fdr * dx, fy = fdr * dy, fz = fdr * dz;
                fxi += fx; fyi += fy; fzi += fz;
                fx_a[j] -= fx; fy_a[j] -= fy; fz_a[j] -= fz;
                a.pe += u;
                a.vir += dx * fx + dy * fy + dz * fz;
                a.wxx += dx * fx; a.wyy += dy * fy; a.wzz += dz * fz;
                a.wxy += dx * fy; a.wxz += dx * fz; a.wyz += dy * fz;
            }
        }
        fx_a[i] = fxi; fy_a[i] = fyi; fz_a[i] = fzi;
    }
    md_soa_forces_to_system(sys, soa);
    finish_system_forces(sys, pot, &a);
}

typedef struct {
    int i_start, i_end;
    const md_system *sys;
    const md_potential *pot;
    const md_verlet_list *vl;
    md_vec3 *force_out;
    force_acc *per_i;
} md_thr_work;

static void *forces_worker(void *arg) {
    md_thr_work *w = (md_thr_work *)arg;
    const md_system *sys = w->sys;
    const md_potential *pot = w->pot;
    const md_verlet_list *vl = w->vl;
    double rc2 = pot->rc2;

    for (int i = w->i_start; i < w->i_end; ++i) {
        md_vec3 pi = sys->pos[i];
        double fxi = 0.0, fyi = 0.0, fzi = 0.0;
        force_acc a = {0};
        int start = vl->full_offsets[i], end = vl->full_offsets[i + 1];
        for (int k = start; k < end; ++k) {
            int j = vl->full_nbrs[k];
            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
            double r2 = v3_norm2(d);
            if (r2 < rc2) {
                double fdr = 0.0;
                double u = md_pot_eval(pot, r2, &fdr);
                double fx = fdr * d.x, fy = fdr * d.y, fz = fdr * d.z;
                fxi += fx; fyi += fy; fzi += fz;
                a.pe += 0.5 * u;
                a.vir += 0.5 * (d.x * fx + d.y * fy + d.z * fz);
                a.wxx += 0.5 * d.x * fx; a.wyy += 0.5 * d.y * fy; a.wzz += 0.5 * d.z * fz;
                a.wxy += 0.5 * d.x * fy; a.wxz += 0.5 * d.x * fz; a.wyz += 0.5 * d.y * fz;
            }
        }
        w->force_out[i] = v3_make(fxi, fyi, fzi);
        w->per_i[i] = a;
    }
    return NULL;
}

void md_forces_threaded(md_system *sys, const md_potential *pot, md_verlet_list *vl, int n_threads) {
    if (!vl->need_full) {
        vl->need_full = true;
        vl->full_offsets = (int *)calloc((size_t)(sys->n + 1), sizeof(int));
        vl->full_cap = sys->n * 192 + 512;
        vl->full_nbrs = (int *)malloc((size_t)vl->full_cap * sizeof(int));
        vl->built = false;
    }
    (void)md_verlet_check_and_rebuild(vl, sys);
    if (n_threads < 1) n_threads = 1;
    if (n_threads > 16) n_threads = 16;
    if (n_threads > sys->n) n_threads = sys->n;

    int n = sys->n;
    force_acc *per_i = (force_acc *)calloc((size_t)n, sizeof(force_acc));
    md_thr_work works[16];
    pthread_t tids[16];
    for (int t = 0; t < n_threads; ++t) {
        works[t] = (md_thr_work){
            (t * n) / n_threads, ((t + 1) * n) / n_threads, sys, pot, vl, sys->force, per_i
        };
    }
    if (n_threads == 1) {
        forces_worker(&works[0]);
    } else {
        for (int t = 1; t < n_threads; ++t) pthread_create(&tids[t], NULL, forces_worker, &works[t]);
        forces_worker(&works[0]);
        for (int t = 1; t < n_threads; ++t) pthread_join(tids[t], NULL);
    }

    force_acc tot = {0};
    for (int i = 0; i < n; ++i) {
        tot.pe += per_i[i].pe; tot.vir += per_i[i].vir;
        tot.wxx += per_i[i].wxx; tot.wyy += per_i[i].wyy; tot.wzz += per_i[i].wzz;
        tot.wxy += per_i[i].wxy; tot.wxz += per_i[i].wxz; tot.wyz += per_i[i].wyz;
    }
    free(per_i);
    finish_system_forces(sys, pot, &tot);
}

typedef struct {
    int i_start, i_end;
    const md_soa *soa;
    md_box box;
    const md_potential *pot;
    const md_verlet_list *vl;
    const int *inv_perm;
    double *fx_out, *fy_out, *fz_out, *pe_i, *vir_i;
} md_v5_work;

static void *forces_v5_worker(void *arg) {
    md_v5_work *w = (md_v5_work *)arg;
    double rc2 = w->pot->rc2;
    double lx = w->box.lx, ly = w->box.ly, lz = w->box.lz;
    double ilx = w->box.inv_lx, ily = w->box.inv_ly, ilz = w->box.inv_lz;
    int dim = w->box.dim;
    const double *restrict px = w->soa->x, *restrict py = w->soa->y, *restrict pz = w->soa->z;

    for (int s = w->i_start; s < w->i_end; ++s) {
        int orig_i = w->soa->perm[s];
        double xi = px[s], yi = py[s], zi = pz[s];
        double fxi = 0.0, fyi = 0.0, fzi = 0.0, pe = 0.0, vir = 0.0;
        int start = w->vl->full_offsets[orig_i], end = w->vl->full_offsets[orig_i + 1];
        for (int k = start; k < end; ++k) {
            int sj = w->inv_perm[w->vl->full_nbrs[k]];
            double dx = md_mic_1d(xi - px[sj], lx, ilx);
            double dy = md_mic_1d(yi - py[sj], ly, ily);
            double dz = (dim == 3) ? md_mic_1d(zi - pz[sj], lz, ilz) : 0.0;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                double fdr = 0.0;
                double u = md_pot_eval(w->pot, r2, &fdr);
                double fx = fdr * dx, fy = fdr * dy, fz = fdr * dz;
                fxi += fx; fyi += fy; fzi += fz;
                pe += 0.5 * u;
                vir += 0.5 * (dx * fx + dy * fy + dz * fz);
            }
        }
        w->fx_out[s] = fxi; w->fy_out[s] = fyi; w->fz_out[s] = fzi;
        w->pe_i[orig_i] = pe; w->vir_i[orig_i] = vir;
    }
    return NULL;
}

void md_forces_cell_sorted_soa(md_system *sys, md_soa *soa, const md_potential *pot,
                               md_verlet_list *vl, int n_threads) {
    if (!vl->need_full) {
        vl->need_full = true;
        vl->full_offsets = (int *)calloc((size_t)(sys->n + 1), sizeof(int));
        vl->full_cap = sys->n * 192 + 512;
        vl->full_nbrs = (int *)malloc((size_t)vl->full_cap * sizeof(int));
        vl->built = false;
    }
    bool rebuilt = md_verlet_check_and_rebuild(vl, sys);
    if (rebuilt || soa->perm[0] == 0) {
        md_soa_sort_by_cell(soa, sys, vl->r_list);
    } else {
        for (int s = 0; s < sys->n; ++s) {
            int o = soa->perm[s];
            soa->x[s] = sys->pos[o].x; soa->y[s] = sys->pos[o].y; soa->z[s] = sys->pos[o].z;
        }
    }
    int n = sys->n;
    int *inv_perm = (int *)malloc((size_t)n * sizeof(int));
    for (int s = 0; s < n; ++s) inv_perm[soa->perm[s]] = s;
    if (n_threads < 1) n_threads = 1;
    if (n_threads > 16) n_threads = 16;
    if (n_threads > n) n_threads = n;

    double *pe_i = (double *)calloc((size_t)n, sizeof(double));
    double *vir_i = (double *)calloc((size_t)n, sizeof(double));
    md_v5_work works[16];
    pthread_t tids[16];
    for (int t = 0; t < n_threads; ++t) {
        works[t] = (md_v5_work){
            (t * n) / n_threads, ((t + 1) * n) / n_threads, soa, sys->box, pot, vl, inv_perm,
            soa->fx, soa->fy, soa->fz, pe_i, vir_i
        };
    }
    if (n_threads == 1) {
        forces_v5_worker(&works[0]);
    } else {
        for (int t = 1; t < n_threads; ++t) pthread_create(&tids[t], NULL, forces_v5_worker, &works[t]);
        forces_v5_worker(&works[0]);
        for (int t = 1; t < n_threads; ++t) pthread_join(tids[t], NULL);
    }
    md_soa_forces_to_system(sys, soa);
    force_acc tot = {0};
    for (int i = 0; i < n; ++i) { tot.pe += pe_i[i]; tot.vir += vir_i[i]; }
    free(inv_perm); free(pe_i); free(vir_i);
    finish_system_forces(sys, pot, &tot);
}

void md_forces_compute(md_system *sys, const md_potential *pot, md_nbr_kind method,
                       md_cell_list *cl, md_verlet_list *vl, md_soa *soa, int n_threads) {
    switch (method) {
        case NBR_BRUTE:           md_forces_brute(sys, pot); break;
        case NBR_CELL:            md_forces_cell(sys, pot, cl); break;
        case NBR_VERLET:          md_forces_verlet(sys, pot, vl); break;
        case NBR_SOA_VERLET:      md_forces_soa_verlet(sys, soa, pot, vl); break;
        case NBR_THREADED:        md_forces_threaded(sys, pot, vl, n_threads); break;
        case NBR_CELL_SORTED_SOA: md_forces_cell_sorted_soa(sys, soa, pot, vl, n_threads); break;
    }
}

static void forces_respa_common(md_system *sys, const md_potential *pot, md_verlet_list *vl, bool fast) {
    (void)md_verlet_check_and_rebuild(vl, sys);
    zero_system_forces(sys);
    md_potential wca;
    md_potential_init(&wca, POT_WCA, MD_WCA_RC);
    wca.eps = pot->eps; wca.sigma = pot->sigma;
    double rlim2 = fast ? wca.rc2 : pot->rc2;
    double pe = 0.0, vir = 0.0;

    for (int i = 0; i < sys->n; ++i) {
        md_vec3 pi = sys->pos[i];
        int start = vl->half_offsets[i], end = vl->half_offsets[i + 1];
        for (int k = start; k < end; ++k) {
            int j = vl->half_nbrs[k];
            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
            double r2 = v3_norm2(d);
            if (r2 < rlim2) {
                double fdr = 0.0;
                double u = fast ? md_pot_eval_exact(&wca, r2, &fdr) : md_pot_lj_attr_eval(pot, r2, &fdr);
                double fx = fdr * d.x, fy = fdr * d.y, fz = fdr * d.z;
                sys->force[i].x += fx; sys->force[i].y += fy; sys->force[i].z += fz;
                sys->force[j].x -= fx; sys->force[j].y -= fy; sys->force[j].z -= fz;
                pe += u; vir += d.x * fx + d.y * fy + d.z * fz;
            }
        }
    }
    sys->pe = pe;
    sys->virial = vir;
}

void md_forces_respa_fast(md_system *sys, const md_potential *pot, md_verlet_list *vl) {
    forces_respa_common(sys, pot, vl, true);
}

void md_forces_respa_slow(md_system *sys, const md_potential *pot, md_verlet_list *vl) {
    forces_respa_common(sys, pot, vl, false);
}
