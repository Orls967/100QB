#ifndef MD_VEC_H
#define MD_VEC_H

#include "md_types.h"
#include <math.h>

static inline md_vec3 v3_make(double x, double y, double z) {
    md_vec3 v = {x, y, z};
    return v;
}

static inline md_vec3 v3_zero(void) {
    md_vec3 v = {0.0, 0.0, 0.0};
    return v;
}

static inline md_vec3 v3_add(md_vec3 a, md_vec3 b) {
    return v3_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline md_vec3 v3_sub(md_vec3 a, md_vec3 b) {
    return v3_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline md_vec3 v3_scale(md_vec3 a, double s) {
    return v3_make(a.x * s, a.y * s, a.z * s);
}

static inline double v3_dot(md_vec3 a, md_vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline double v3_norm2(md_vec3 a) {
    return a.x * a.x + a.y * a.y + a.z * a.z;
}

static inline double v3_norm(md_vec3 a) {
    return sqrt(v3_norm2(a));
}

/* Box initialization and volume */
md_box md_box_cubic(double l);
md_box md_box_ortho(double lx, double ly, double lz, int dim);
double md_box_volume(const md_box *b);

/* Wrap coordinate into [0, L) */
static inline double md_wrap_1d(double x, double l, double inv_l) {
    double w = x - floor(x * inv_l) * l;
    if (w >= l) w -= l;
    if (w < 0.0) w += l;
    return w;
}

static inline md_vec3 md_box_wrap(const md_box *b, md_vec3 p) {
    md_vec3 w;
    w.x = md_wrap_1d(p.x, b->lx, b->inv_lx);
    w.y = md_wrap_1d(p.y, b->ly, b->inv_ly);
    w.z = (b->dim == 3) ? md_wrap_1d(p.z, b->lz, b->inv_lz) : 0.0;
    return w;
}

/* Minimum image convention delta = r_i - r_j in [-L/2, L/2) */
static inline double md_mic_1d(double d, double l, double inv_l) {
    return d - round(d * inv_l) * l;
}

static inline md_vec3 md_box_mic(const md_box *b, md_vec3 a, md_vec3 c) {
    md_vec3 d;
    d.x = md_mic_1d(a.x - c.x, b->lx, b->inv_lx);
    d.y = md_mic_1d(a.y - c.y, b->ly, b->inv_ly);
    d.z = (b->dim == 3) ? md_mic_1d(a.z - c.z, b->lz, b->inv_lz) : 0.0;
    return d;
}

/* System lifecycle and helpers */
void md_system_alloc(md_system *sys, int n, md_box box);
void md_system_free(md_system *sys);
void md_system_copy(md_system *dst, const md_system *src);
void md_system_wrap_all(md_system *sys);
void md_system_rescale_box(md_system *sys, double scale);

/* SoA lifecycle and synchronization */
void md_soa_alloc(md_soa *soa, int n);
void md_soa_free(md_soa *soa);
void md_soa_from_system(md_soa *soa, const md_system *sys);
void md_soa_forces_to_system(md_system *sys, const md_soa *soa);

#endif /* MD_VEC_H */
