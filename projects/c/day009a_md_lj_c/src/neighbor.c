#include "md_neighbor.h"
#include "md_vec.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void md_cell_init(md_cell_list *cl, int n) {
    memset(cl, 0, sizeof(*cl));
    cl->cap_n = n;
    cl->next = (int *)calloc((size_t)n, sizeof(int));
    cl->cell_of = (int *)calloc((size_t)n, sizeof(int));
}

void md_cell_free(md_cell_list *cl) {
    free(cl->head);
    free(cl->next);
    free(cl->cell_of);
    memset(cl, 0, sizeof(*cl));
}

void md_cell_build(md_cell_list *cl, const md_system *sys, double r_cut) {
    if (sys->n > cl->cap_n) {
        free(cl->next);
        free(cl->cell_of);
        cl->cap_n = sys->n;
        cl->next = (int *)calloc((size_t)sys->n, sizeof(int));
        cl->cell_of = (int *)calloc((size_t)sys->n, sizeof(int));
    }
    int nx = (int)floor(sys->box.lx / r_cut);
    int ny = (int)floor(sys->box.ly / r_cut);
    int nz = (sys->box.dim == 3) ? (int)floor(sys->box.lz / r_cut) : 1;
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    if (nz < 1) nz = 1;
    /* Cap excessive cells for very dilute boxes to keep memory compact */
    if (nx > 64) nx = 64;
    if (ny > 64) ny = 64;
    if (nz > 64) nz = 64;

    int total = nx * ny * nz;
    if (total != cl->n_cells || cl->head == NULL) {
        free(cl->head);
        cl->head = (int *)malloc((size_t)total * sizeof(int));
        cl->n_cells = total;
    }
    cl->nx = nx;
    cl->ny = ny;
    cl->nz = nz;

    for (int c = 0; c < total; ++c) {
        cl->head[c] = -1;
    }

    double sx = (double)nx * sys->box.inv_lx;
    double sy = (double)ny * sys->box.inv_ly;
    double sz = (sys->box.dim == 3) ? ((double)nz * sys->box.inv_lz) : 0.0;

    for (int i = 0; i < sys->n; ++i) {
        md_vec3 p = md_box_wrap(&sys->box, sys->pos[i]);
        int cx = (int)floor(p.x * sx);
        int cy = (int)floor(p.y * sy);
        int cz = (sys->box.dim == 3) ? (int)floor(p.z * sz) : 0;
        if (cx < 0) cx = 0; else if (cx >= nx) cx = nx - 1;
        if (cy < 0) cy = 0; else if (cy >= ny) cy = ny - 1;
        if (cz < 0) cz = 0; else if (cz >= nz) cz = nz - 1;
        int cid = (cz * ny + cy) * nx + cx;
        cl->cell_of[i] = cid;
        cl->next[i] = cl->head[cid];
        cl->head[cid] = i;
    }
}

static void verlet_ensure_half_cap(md_verlet_list *vl, int needed) {
    if (needed <= vl->half_cap) return;
    int new_cap = (vl->half_cap > 0) ? vl->half_cap * 2 : 1024;
    while (new_cap < needed) new_cap *= 2;
    int *p = (int *)realloc(vl->half_nbrs, (size_t)new_cap * sizeof(int));
    if (!p) {
        vl->overflow_count++;
        return;
    }
    vl->half_nbrs = p;
    vl->half_cap = new_cap;
}

static void verlet_ensure_full_cap(md_verlet_list *vl, int needed) {
    if (needed <= vl->full_cap) return;
    int new_cap = (vl->full_cap > 0) ? vl->full_cap * 2 : 2048;
    while (new_cap < needed) new_cap *= 2;
    int *p = (int *)realloc(vl->full_nbrs, (size_t)new_cap * sizeof(int));
    if (!p) {
        vl->overflow_count++;
        return;
    }
    vl->full_nbrs = p;
    vl->full_cap = new_cap;
}

static int cmp_int_asc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}

void md_verlet_init(md_verlet_list *vl, int n, double r_cut, double r_skin, bool need_full) {
    memset(vl, 0, sizeof(*vl));
    vl->n = n;
    vl->r_cut = r_cut;
    vl->r_skin = r_skin;
    vl->r_list = r_cut + r_skin;
    vl->r_list2 = vl->r_list * vl->r_list;
    vl->need_full = need_full;

    vl->half_offsets = (int *)calloc((size_t)(n + 1), sizeof(int));
    vl->half_cap = n * 96 + 256;
    vl->half_nbrs = (int *)malloc((size_t)vl->half_cap * sizeof(int));

    if (need_full) {
        vl->full_offsets = (int *)calloc((size_t)(n + 1), sizeof(int));
        vl->full_cap = n * 192 + 512;
        vl->full_nbrs = (int *)malloc((size_t)vl->full_cap * sizeof(int));
    }

    vl->ref_pos = (md_vec3 *)calloc((size_t)n, sizeof(md_vec3));
    md_cell_init(&vl->cells, n);
}

void md_verlet_free(md_verlet_list *vl) {
    free(vl->half_offsets);
    free(vl->half_nbrs);
    free(vl->full_offsets);
    free(vl->full_nbrs);
    free(vl->ref_pos);
    md_cell_free(&vl->cells);
    memset(vl, 0, sizeof(*vl));
}

void md_verlet_build(md_verlet_list *vl, const md_system *sys) {
    int n = sys->n;
    double rl2 = vl->r_list2;

    md_cell_build(&vl->cells, sys, vl->r_list);
    int nx = vl->cells.nx;
    int ny = vl->cells.ny;
    int nz = vl->cells.nz;
    bool use_cells = (nx >= 3 && ny >= 3 && (sys->box.dim == 2 || nz >= 3));

    int half_pos = 0;
    int full_pos = 0;
    vl->half_offsets[0] = 0;
    if (vl->need_full) {
        vl->full_offsets[0] = 0;
    }

    for (int i = 0; i < n; ++i) {
        md_vec3 pi = sys->pos[i];
        int half_start = half_pos;
        int full_start = full_pos;

        if (!use_cells) {
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
                if (v3_norm2(d) < rl2) {
                    if (j > i) {
                        verlet_ensure_half_cap(vl, half_pos + 1);
                        if (half_pos < vl->half_cap) vl->half_nbrs[half_pos++] = j;
                    }
                    if (vl->need_full) {
                        verlet_ensure_full_cap(vl, full_pos + 1);
                        if (full_pos < vl->full_cap) vl->full_nbrs[full_pos++] = j;
                    }
                }
            }
        } else {
            int cid = vl->cells.cell_of[i];
            int cx = cid % nx;
            int cy = (cid / nx) % ny;
            int cz = cid / (nx * ny);
            int dz_lo = (sys->box.dim == 3) ? -1 : 0;
            int dz_hi = (sys->box.dim == 3) ? 1 : 0;

            for (int dz = dz_lo; dz <= dz_hi; ++dz) {
                int ncz = (cz + dz + nz) % nz;
                for (int dy = -1; dy <= 1; ++dy) {
                    int ncy = (cy + dy + ny) % ny;
                    for (int dx = -1; dx <= 1; ++dx) {
                        int ncx = (cx + dx + nx) % nx;
                        int ncid = (ncz * ny + ncy) * nx + ncx;
                        for (int j = vl->cells.head[ncid]; j != -1; j = vl->cells.next[j]) {
                            if (j == i) continue;
                            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
                            if (v3_norm2(d) < rl2) {
                                if (j > i) {
                                    verlet_ensure_half_cap(vl, half_pos + 1);
                                    if (half_pos < vl->half_cap) vl->half_nbrs[half_pos++] = j;
                                }
                                if (vl->need_full) {
                                    verlet_ensure_full_cap(vl, full_pos + 1);
                                    if (full_pos < vl->full_cap) vl->full_nbrs[full_pos++] = j;
                                }
                            }
                        }
                    }
                }
            }
            /* Sort per-particle lists in strictly ascending index order for locked summation */
            if (half_pos - half_start > 1) {
                qsort(vl->half_nbrs + half_start, (size_t)(half_pos - half_start),
                      sizeof(int), cmp_int_asc);
            }
            if (vl->need_full && (full_pos - full_start > 1)) {
                qsort(vl->full_nbrs + full_start, (size_t)(full_pos - full_start),
                      sizeof(int), cmp_int_asc);
            }
        }
        vl->half_offsets[i + 1] = half_pos;
        if (vl->need_full) {
            vl->full_offsets[i + 1] = full_pos;
        }
        vl->ref_pos[i] = sys->unwrapped[i];
    }

    vl->half_count = half_pos;
    vl->full_count = full_pos;
    vl->built = true;
    vl->rebuild_count++;
    vl->last_max_disp = 0.0;
}

bool md_verlet_check_and_rebuild(md_verlet_list *vl, const md_system *sys) {
    vl->check_count++;
    if (!vl->built) {
        md_verlet_build(vl, sys);
        return true;
    }
    double max_d2 = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        md_vec3 d = v3_sub(sys->unwrapped[i], vl->ref_pos[i]);
        double d2 = v3_norm2(d);
        if (d2 > max_d2) max_d2 = d2;
    }
    double max_d = sqrt(max_d2);
    vl->last_max_disp = max_d;
    if (max_d > vl->max_disp_before_rebuild) {
        vl->max_disp_before_rebuild = max_d;
    }
    /* Strict mathematical safety limit: max_disp <= r_skin / 2.
     * Trigger rebuild proactively at 0.35 * r_skin so even fast thermal tails never cross r_skin/2. */
    if (max_d > 0.5 * vl->r_skin) {
        vl->late_rebuild_count++;
    }
    if (max_d >= 0.35 * vl->r_skin) {
        md_verlet_build(vl, sys);
        return true;
    }
    return false;
}

void md_soa_sort_by_cell(md_soa *soa, const md_system *sys, double r_cut) {
    md_cell_list cl;
    md_cell_init(&cl, sys->n);
    md_cell_build(&cl, sys, r_cut);

    int slot = 0;
    for (int c = 0; c < cl.n_cells; ++c) {
        for (int i = cl.head[c]; i != -1; i = cl.next[i]) {
            soa->x[slot] = sys->pos[i].x;
            soa->y[slot] = sys->pos[i].y;
            soa->z[slot] = sys->pos[i].z;
            soa->ux[slot] = sys->unwrapped[i].x;
            soa->uy[slot] = sys->unwrapped[i].y;
            soa->uz[slot] = sys->unwrapped[i].z;
            soa->vx[slot] = sys->vel[i].x;
            soa->vy[slot] = sys->vel[i].y;
            soa->vz[slot] = sys->vel[i].z;
            soa->fx[slot] = 0.0;
            soa->fy[slot] = 0.0;
            soa->fz[slot] = 0.0;
            soa->perm[slot] = i;
            slot++;
        }
    }
    md_cell_free(&cl);
}
