#include "md_vec.h"
#include <stdlib.h>
#include <string.h>

md_box md_box_cubic(double l) {
    md_box b;
    b.lx = l;
    b.ly = l;
    b.lz = l;
    b.inv_lx = 1.0 / l;
    b.inv_ly = 1.0 / l;
    b.inv_lz = 1.0 / l;
    b.dim = 3;
    return b;
}

md_box md_box_ortho(double lx, double ly, double lz, int dim) {
    md_box b;
    b.lx = lx;
    b.ly = ly;
    b.lz = (dim == 3) ? lz : 1.0;
    b.inv_lx = 1.0 / lx;
    b.inv_ly = 1.0 / ly;
    b.inv_lz = (dim == 3) ? (1.0 / lz) : 1.0;
    b.dim = dim;
    return b;
}

double md_box_volume(const md_box *b) {
    if (b->dim == 2) {
        return b->lx * b->ly;
    }
    return b->lx * b->ly * b->lz;
}

void md_system_alloc(md_system *sys, int n, md_box box) {
    memset(sys, 0, sizeof(*sys));
    sys->n = n;
    sys->box = box;
    sys->pos = (md_vec3 *)calloc((size_t)n, sizeof(md_vec3));
    sys->unwrapped = (md_vec3 *)calloc((size_t)n, sizeof(md_vec3));
    sys->vel = (md_vec3 *)calloc((size_t)n, sizeof(md_vec3));
    sys->force = (md_vec3 *)calloc((size_t)n, sizeof(md_vec3));
}

void md_system_free(md_system *sys) {
    free(sys->pos);
    free(sys->unwrapped);
    free(sys->vel);
    free(sys->force);
    memset(sys, 0, sizeof(*sys));
}

void md_system_copy(md_system *dst, const md_system *src) {
    if (dst->n != src->n) {
        md_system_free(dst);
        md_system_alloc(dst, src->n, src->box);
    }
    dst->box = src->box;
    size_t sz = (size_t)src->n * sizeof(md_vec3);
    memcpy(dst->pos, src->pos, sz);
    memcpy(dst->unwrapped, src->unwrapped, sz);
    memcpy(dst->vel, src->vel, sz);
    memcpy(dst->force, src->force, sz);
    dst->pe = src->pe;
    dst->pe_tail = src->pe_tail;
    dst->virial = src->virial;
    dst->p_tail = src->p_tail;
    dst->w_xx = src->w_xx;
    dst->w_yy = src->w_yy;
    dst->w_zz = src->w_zz;
    dst->w_xy = src->w_xy;
    dst->w_xz = src->w_xz;
    dst->w_yz = src->w_yz;
}

void md_system_wrap_all(md_system *sys) {
    for (int i = 0; i < sys->n; ++i) {
        sys->pos[i] = md_box_wrap(&sys->box, sys->pos[i]);
    }
}

void md_system_rescale_box(md_system *sys, double scale) {
    sys->box.lx *= scale;
    sys->box.ly *= scale;
    sys->box.inv_lx = 1.0 / sys->box.lx;
    sys->box.inv_ly = 1.0 / sys->box.ly;
    if (sys->box.dim == 3) {
        sys->box.lz *= scale;
        sys->box.inv_lz = 1.0 / sys->box.lz;
    }
    for (int i = 0; i < sys->n; ++i) {
        sys->pos[i].x *= scale;
        sys->pos[i].y *= scale;
        sys->unwrapped[i].x *= scale;
        sys->unwrapped[i].y *= scale;
        if (sys->box.dim == 3) {
            sys->pos[i].z *= scale;
            sys->unwrapped[i].z *= scale;
        }
    }
    md_system_wrap_all(sys);
}

void md_soa_alloc(md_soa *soa, int n) {
    memset(soa, 0, sizeof(*soa));
    soa->n = n;
    soa->x = (double *)calloc((size_t)n, sizeof(double));
    soa->y = (double *)calloc((size_t)n, sizeof(double));
    soa->z = (double *)calloc((size_t)n, sizeof(double));
    soa->ux = (double *)calloc((size_t)n, sizeof(double));
    soa->uy = (double *)calloc((size_t)n, sizeof(double));
    soa->uz = (double *)calloc((size_t)n, sizeof(double));
    soa->vx = (double *)calloc((size_t)n, sizeof(double));
    soa->vy = (double *)calloc((size_t)n, sizeof(double));
    soa->vz = (double *)calloc((size_t)n, sizeof(double));
    soa->fx = (double *)calloc((size_t)n, sizeof(double));
    soa->fy = (double *)calloc((size_t)n, sizeof(double));
    soa->fz = (double *)calloc((size_t)n, sizeof(double));
    soa->perm = (int *)calloc((size_t)n, sizeof(int));
    for (int i = 0; i < n; ++i) {
        soa->perm[i] = i;
    }
}

void md_soa_free(md_soa *soa) {
    free(soa->x);
    free(soa->y);
    free(soa->z);
    free(soa->ux);
    free(soa->uy);
    free(soa->uz);
    free(soa->vx);
    free(soa->vy);
    free(soa->vz);
    free(soa->fx);
    free(soa->fy);
    free(soa->fz);
    free(soa->perm);
    memset(soa, 0, sizeof(*soa));
}

void md_soa_from_system(md_soa *soa, const md_system *sys) {
    for (int i = 0; i < sys->n; ++i) {
        soa->x[i] = sys->pos[i].x;
        soa->y[i] = sys->pos[i].y;
        soa->z[i] = sys->pos[i].z;
        soa->ux[i] = sys->unwrapped[i].x;
        soa->uy[i] = sys->unwrapped[i].y;
        soa->uz[i] = sys->unwrapped[i].z;
        soa->vx[i] = sys->vel[i].x;
        soa->vy[i] = sys->vel[i].y;
        soa->vz[i] = sys->vel[i].z;
        soa->fx[i] = 0.0;
        soa->fy[i] = 0.0;
        soa->fz[i] = 0.0;
        soa->perm[i] = i;
    }
}

void md_soa_forces_to_system(md_system *sys, const md_soa *soa) {
    for (int i = 0; i < sys->n; ++i) {
        int orig = soa->perm[i];
        sys->force[orig].x = soa->fx[i];
        sys->force[orig].y = soa->fy[i];
        sys->force[orig].z = soa->fz[i];
    }
}
