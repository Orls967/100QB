#include "common.h"

static void binom_probs(int m, double mu, double *out) {
    out[0] = 1.0;
    for (int i = 1; i <= m; ++i) out[i] = 0.0;
    for (int s = 0; s < m; ++s) {
        for (int k = m; k >= 1; --k) out[k] = out[k] * (1.0 - mu) + out[k - 1] * mu;
        out[0] *= (1.0 - mu);
    }
}

DtmcRes solve_dtmc_single(int c, int K, double lam, double mu) {
    double P[64][64]; memset(P, 0, sizeof(P));
    double bp[8];
    for (int n = 0; n <= K; ++n) {
        int m = n < c ? n : c;
        binom_probs(m, mu, bp);
        for (int d = 0; d <= m; ++d) {
            if (n < K) {
                P[n][n - d] += (1.0 - lam) * bp[d];
                P[n][n + 1 - d] += lam * bp[d];
            } else {
                P[n][n - d] += bp[d]; /* arrival blocked when n == K */
            }
        }
    }
    double pi[64];
    for (int i = 0; i <= K; ++i) pi[i] = 1.0 / (double)(K + 1);
    for (int it = 0; it < 25000; ++it) {
        double maxd = 0.0, sum = 0.0;
        for (int j = 0; j <= K; ++j) {
            double nv = 0.0;
            for (int i = 0; i <= K; ++i) nv += pi[i] * P[i][j];
            double df = fabs(nv - pi[j]);
            if (df > maxd) maxd = df;
            pi[j] = nv; sum += nv;
        }
        for (int j = 0; j <= K; ++j) pi[j] /= sum;
        if (maxd < 1e-14) break;
    }
    DtmcRes r; memset(&r, 0, sizeof(r));
    r.p_block = pi[K];
    r.lam_eff = lam * (1.0 - pi[K]);
    for (int n = 0; n <= K; ++n) {
        int m = n < c ? n : c;
        r.L += (double)n * pi[n];
        r.L_q += (double)(n - m) * pi[n];
        r.U += ((double)m / (double)c) * pi[n];
    }
    r.W = r.lam_eff > 0.0 ? r.L / r.lam_eff : 0.0;
    r.W_q = r.lam_eff > 0.0 ? r.L_q / r.lam_eff : 0.0;
    return r;
}

DtmcRes solve_geo_closed_form(double lam, double mu) {
    DtmcRes r; memset(&r, 0, sizeof(r));
    r.p_block = 0.0;
    r.lam_eff = lam;
    r.U = lam / mu;
    r.L = (lam * (1.0 - lam)) / (mu - lam);
    r.W = (1.0 - lam) / (mu - lam);
    r.W_q = r.W - 1.0 / mu;
    r.L_q = lam * r.W_q;
    return r;
}

void solve_dtmc_tandem(int K1, int K2, double lam, double mu1, double mu2, DtmcRes *r1, DtmcRes *r2, double *p_blk_tnd) {
    int S = (K1 + 1) * (K2 + 1);
    double P[64][64]; memset(P, 0, sizeof(P));
    for (int n1 = 0; n1 <= K1; ++n1) {
        for (int n2 = 0; n2 <= K2; ++n2) {
            int s_idx = n1 * (K2 + 1) + n2;
            for (int a = 0; a <= 1; ++a) {
                double pa = a ? lam : (1.0 - lam);
                for (int d1 = 0; d1 <= (n1 > 0 ? 1 : 0); ++d1) {
                    double pd1 = (n1 == 0) ? 1.0 : (d1 ? mu1 : (1.0 - mu1));
                    for (int d2 = 0; d2 <= (n2 > 0 ? 1 : 0); ++d2) {
                        double pd2 = (n2 == 0) ? 1.0 : (d2 ? mu2 : (1.0 - mu2));
                        int nn1 = (n1 < K1) ? (n1 + a - d1) : (n1 - d1);
                        int n2_mid = n2 - d2;
                        int nn2 = (n2_mid < K2) ? (n2_mid + d1) : n2_mid;
                        P[s_idx][nn1 * (K2 + 1) + nn2] += pa * pd1 * pd2;
                    }
                }
            }
        }
    }
    double pi[64];
    for (int i = 0; i < S; ++i) pi[i] = 1.0 / (double)S;
    for (int it = 0; it < 25000; ++it) {
        double maxd = 0.0, sum = 0.0;
        for (int j = 0; j < S; ++j) {
            double nv = 0.0;
            for (int i = 0; i < S; ++i) nv += pi[i] * P[i][j];
            double df = fabs(nv - pi[j]);
            if (df > maxd) maxd = df;
            pi[j] = nv; sum += nv;
        }
        for (int j = 0; j < S; ++j) pi[j] /= sum;
        if (maxd < 1e-14) break;
    }
    memset(r1, 0, sizeof(*r1)); memset(r2, 0, sizeof(*r2));
    double blk1 = 0.0, blk2 = 0.0;
    for (int n1 = 0; n1 <= K1; ++n1) {
        for (int n2 = 0; n2 <= K2; ++n2) {
            double p = pi[n1 * (K2 + 1) + n2];
            r1->L += (double)n1 * p; r2->L += (double)n2 * p;
            if (n1 > 0) { r1->U += p; r1->L_q += (double)(n1 - 1) * p; }
            if (n2 > 0) { r2->U += p; r2->L_q += (double)(n2 - 1) * p; }
            if (n1 == K1) blk1 += p * lam;
            if (n1 > 0 && n2 == K2) blk2 += p * mu1 * (1.0 - mu2);
        }
    }
    r1->p_block = blk1; r1->lam_eff = lam - blk1;
    r1->W = r1->lam_eff > 0.0 ? r1->L / r1->lam_eff : 0.0;
    r1->W_q = r1->lam_eff > 0.0 ? r1->L_q / r1->lam_eff : 0.0;
    r2->p_block = blk2; r2->lam_eff = r1->lam_eff - blk2;
    r2->W = r2->lam_eff > 0.0 ? r2->L / r2->lam_eff : 0.0;
    r2->W_q = r2->lam_eff > 0.0 ? r2->L_q / r2->lam_eff : 0.0;
    *p_blk_tnd = blk2 / (r1->lam_eff > 0.0 ? r1->lam_eff : 1.0);
}
