#include "sandpile.h"
#include <stdlib.h>
#include <string.h>

void fmt_i128(i128 v, char *buf) {
    if (v == 0) { strcpy(buf, "0"); return; }
    int neg = v < 0, p = 0;
    __extension__ unsigned __int128 u = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
    char tmp[48];
    while (u > 0) { tmp[p++] = (char)('0' + (u % 10)); u /= 10; }
    if (neg) tmp[p++] = '-';
    for (int i = 0; i < p; i++) buf[i] = tmp[p - 1 - i];
    buf[p] = '\0';
}

static i128 bareiss_laplacian_det(int n, int *ovf_out) {
    int m = n * n;
    i128 a[49][49];
    memset(a, 0, sizeof(a));
    for (int i = 0; i < m; i++) {
        a[i][i] = 4;
        int r = i / n, c = i % n;
        if (r > 0) a[i][i - n] = -1;
        if (r + 1 < n) a[i][i + n] = -1;
        if (c > 0) a[i][i - 1] = -1;
        if (c + 1 < n) a[i][i + 1] = -1;
    }
    *ovf_out = 0;
    i128 prev = 1;
    for (int k = 0; k < m - 1; k++) {
        i128 piv = a[k][k];
        for (int i = k + 1; i < m; i++) {
            for (int j = k + 1; j < m; j++) {
                i128 p1 = 0, p2 = 0, df = 0;
                if (__builtin_mul_overflow(piv, a[i][j], &p1) ||
                    __builtin_mul_overflow(a[i][k], a[k][j], &p2) ||
                    __builtin_sub_overflow(p1, p2, &df)) { *ovf_out = 1; return -1; }
                a[i][j] = df / prev;
            }
        }
        prev = piv;
    }
    return a[m - 1][m - 1];
}

void compute_identity(int n, int32_t *e_out) {
    int m = n * n;
    Grid g; grid_init(&g, n);
    for (int i = 0; i < m; i++) g.z[i] = 6;
    stabilize_det(&g);
    for (int i = 0; i < m; i++) g.z[i] = 6 - g.z[i];
    stabilize_det(&g);
    memcpy(e_out, g.z, (size_t)m * sizeof(int32_t));
    grid_free(&g);
}

static int check_d4_symmetry(int n, const int32_t *e) {
    for (int r = 0; r < n; r++)
        for (int c = 0; c < n; c++) {
            int32_t v = e[r * n + c];
            if (e[r * n + (n - 1 - c)] != v || e[(n - 1 - r) * n + c] != v ||
                e[(n - 1 - r) * n + (n - 1 - c)] != v || e[c * n + r] != v ||
                e[c * n + (n - 1 - r)] != v || e[(n - 1 - c) * n + r] != v ||
                e[(n - 1 - c) * n + (n - 1 - r)] != v) return 0;
        }
    return 1;
}

int cmd_abelian(Metrics *m, const char *out_dir) {
    int ns[3] = {8, 16, 24}, tot_cfg = 0, mis = 0;
    uint64_t digest = fnv1a_init();
    Rng rng_init, rng_ord;
    rng_seed(&rng_init, MASTER_SEED ^ 0xA11CEULL);
    rng_seed(&rng_ord, MASTER_SEED ^ 0xB0B55ULL);
    char path[128]; snprintf(path, sizeof(path), "%s/abelian_check.csv", out_dir);
    FILE *fp = fopen(path, "w");
    if (fp) fprintf(fp, "n,configs,orders_per_cfg,mismatches,digest_fnv1a\n");
    for (int idx = 0; idx < 3; idx++) {
        int n = ns[idx], sz = n * n, n_mis = 0;
        int32_t *init_z = (int32_t *)malloc((size_t)sz * sizeof(int32_t));
        Grid gd, gr; grid_init(&gd, n); grid_init(&gr, n);
        uint64_t sub_dig = fnv1a_init();
        for (int c = 0; c < 200; c++) {
            for (int i = 0; i < sz; i++) init_z[i] = rng_int(&rng_init, 0, 9);
            grid_clear(&gd); memcpy(gd.z, init_z, (size_t)sz * sizeof(int32_t));
            stabilize_det(&gd);
            uint64_t h_ref = grid_state_hash(&gd);
            grid_clear(&gr); memcpy(gr.z, init_z, (size_t)sz * sizeof(int32_t));
            stabilize_waves(&gr);
            if (grid_state_hash(&gr) != h_ref) n_mis++;
            for (int o = 0; o < 100; o++) {
                grid_clear(&gr); memcpy(gr.z, init_z, (size_t)sz * sizeof(int32_t));
                stabilize_rand(&gr, &rng_ord);
                if (grid_state_hash(&gr) != h_ref) n_mis++;
            }
            sub_dig = fnv1a_bytes(sub_dig, &h_ref, sizeof(h_ref));
            digest = fnv1a_bytes(digest, &h_ref, sizeof(h_ref));
            tot_cfg++;
        }
        mis += n_mis;
        if (fp) fprintf(fp, "%d,200,100,%d,%016llx\n", n, n_mis, (unsigned long long)sub_dig);
        grid_free(&gd); grid_free(&gr); free(init_z);
    }
    if (fp) fclose(fp);
    m->abelian_configs = tot_cfg; m->abelian_orders = 100;
    m->abelian_mismatches = mis; m->abelian_digest = digest;
    return mis == 0 ? 0 : 1;
}

int cmd_group(Metrics *m, const char *out_dir) {
    m->safe_det_n = 0;
    for (int n = 2; n <= 7; n++) {
        int ov = 0; i128 dv = bareiss_laplacian_det(n, &ov);
        if (n <= 6) { m->det_val[n] = dv; m->det_ovf[n] = ov; }
        if (!ov && dv > 0) m->safe_det_n = n;
    }
    for (int n = 2; n <= 3; n++) {
        int sz = n * n, tot = 1 << (2 * sz), stk[16]; uint8_t deg[16];
        Grid g; grid_init(&g, n);
        int64_t cnt = 0;
        for (int mask = 0; mask < tot; mask++) {
            for (int i = 0; i < sz; i++) g.z[i] = (mask >> (2 * i)) & 3;
            if (is_recurrent(&g, stk, deg)) cnt++;
        }
        if (n == 2) m->bf_n2 = cnt; else m->bf_n3 = cnt;
        grid_free(&g);
    }
    char path[128], dbuf[48]; snprintf(path, sizeof(path), "%s/group_order.csv", out_dir);
    FILE *fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "n,det_bareiss,i128_overflow,brute_force,match_bf\n");
        for (int n = 2; n <= 6; n++) {
            fmt_i128(m->det_val[n], dbuf);
            int64_t bf = (n == 2) ? m->bf_n2 : (n == 3 ? m->bf_n3 : -1);
            int ok = (bf < 0) ? (!m->det_ovf[n] && m->det_val[n] > 0) : (m->det_val[n] == bf);
            fprintf(fp, "%d,%s,%d,%lld,%d\n", n, dbuf, m->det_ovf[n], (long long)bf, ok);
        }
        fclose(fp);
    }
    int n_id = 24, sz = n_id * n_id;
    int32_t *e = (int32_t *)malloc((size_t)sz * sizeof(int32_t));
    int32_t *stk = (int32_t *)malloc((size_t)sz * sizeof(int32_t));
    uint8_t *deg = (uint8_t *)malloc((size_t)sz * sizeof(uint8_t));
    compute_identity(n_id, e);
    Grid g, gr; grid_init(&g, n_id); grid_init(&gr, n_id);
    memcpy(g.z, e, (size_t)sz * sizeof(int32_t));
    m->id_recurrent = is_recurrent(&g, stk, deg);
    m->id_d4_ok = check_d4_symmetry(n_id, e);
    for (int i = 0; i < sz; i++) g.z[i] = e[i] + e[i];
    stabilize_det(&g);
    m->id_ee_ok = (memcmp(g.z, e, (size_t)sz * sizeof(int32_t)) == 0);
    Rng rng; rng_seed(&rng, MASTER_SEED ^ 0x1DE17ULL);
    m->id_er_tested = 200; m->id_er_ok = 0;
    for (int k = 0; k < 200; k++) {
        do {
            for (int i = 0; i < sz; i++) gr.z[i] = rng_int(&rng, 0, 9);
            stabilize_det(&gr);
        } while (!is_recurrent(&gr, stk, deg));
        for (int i = 0; i < sz; i++) g.z[i] = e[i] + gr.z[i];
        stabilize_det(&g);
        if (memcmp(g.z, gr.z, (size_t)sz * sizeof(int32_t)) == 0) m->id_er_ok++;
    }
    grid_free(&g); grid_free(&gr); free(e); free(stk); free(deg);

    int32_t *e128 = (int32_t *)malloc(128 * 128 * sizeof(int32_t));
    compute_identity(128, e128);
    snprintf(path, sizeof(path), "%s/identity_e_128.csv", out_dir);
    fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "row,col_csv\n");
        for (int r = 0; r < 128; r++) {
            fprintf(fp, "%d,", r);
            for (int c = 0; c < 128; c++) fprintf(fp, "%d%s", e128[r * 128 + c], c == 127 ? "\n" : ",");
        }
        fclose(fp);
    }
    free(e128);

    m->ord_tested = 0; m->ord_div_ok = 0;
    for (int n = 2; n <= 4; n++) {
        int s = n * n, en[16], stkn[16]; uint8_t degn[16];
        int64_t grp_ord = (int64_t)m->det_val[n];
        compute_identity(n, en);
        Grid ga, gc; grid_init(&ga, n); grid_init(&gc, n);
        for (int t = 0; t < 5; t++) {
            do {
                for (int i = 0; i < s; i++) ga.z[i] = rng_int(&rng, 0, 9);
                stabilize_det(&ga);
            } while (!is_recurrent(&ga, stkn, degn));
            memcpy(gc.z, ga.z, (size_t)s * sizeof(int32_t));
            int64_t ord = 1;
            while (memcmp(gc.z, en, (size_t)s * sizeof(int32_t)) != 0 && ord <= grp_ord) {
                for (int i = 0; i < s; i++) gc.z[i] += ga.z[i];
                stabilize_det(&gc); ord++;
            }
            m->ord_tested++;
            if (ord <= grp_ord && (grp_ord % ord) == 0) m->ord_div_ok++;
        }
        grid_free(&ga); grid_free(&gc);
    }
    return 0;
}

/* Hand n=2: L = 4I_4 - A(C_4) = [[4,-1,-1,0],[-1,4,0,-1],[-1,0,4,-1],[0,-1,-1,4]].
 * eig(A(C_4)) = {2, 0, 0, -2} => eig(L) = {2, 4, 4, 6} => det(L) = 2*4*4*6 = 192. */
int cmd_selftest(Metrics *m) {
    int ovf2 = 0, ovf3 = 0;
    i128 d2 = bareiss_laplacian_det(2, &ovf2), d3 = bareiss_laplacian_det(3, &ovf3);
    if (ovf2 || d2 != 192 || ovf3) return 1;
    int32_t e3[9], stk[9]; uint8_t deg[9];
    compute_identity(3, e3);
    Grid g, gr; grid_init(&g, 3); grid_init(&gr, 3);
    int64_t rec_cnt = 0;
    for (int mask = 0; mask < (1 << 18); mask++) {
        for (int i = 0; i < 9; i++) gr.z[i] = (mask >> (2 * i)) & 3;
        if (is_recurrent(&gr, stk, deg)) {
            rec_cnt++;
            for (int i = 0; i < 9; i++) g.z[i] = e3[i] + gr.z[i];
            stabilize_det(&g);
            if (memcmp(g.z, gr.z, 9 * sizeof(int32_t)) != 0) return 2;
        }
    }
    if (rec_cnt != (int64_t)d3) return 3;
    grid_clear(&g); g.z[4] = 16; stabilize_det(&g);
    const int32_t exp_z[9] = {2, 1, 2, 1, 0, 1, 2, 1, 2};
    const int64_t exp_o[9] = {0, 1, 0, 1, 5, 1, 0, 1, 0};
    if (memcmp(g.z, exp_z, sizeof(exp_z)) != 0 ||
        memcmp(g.odom, exp_o, sizeof(exp_o)) != 0 || g.sink != 4) return 4;
    grid_free(&g); grid_free(&gr); m->selftest_ok = 1;
    printf("selftest OK: det(n=2)=192, det(n=3)=%lld, e+r=r for all %lld recurrent configs, hand 3x3 OK\n",
           (long long)d3, (long long)rec_cnt);
    return 0;
}
