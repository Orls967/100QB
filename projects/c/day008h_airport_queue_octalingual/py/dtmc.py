from dataclasses import dataclass


@dataclass
class DtmcRes:
    p_block: float = 0.0
    l: float = 0.0
    l_q: float = 0.0
    u: float = 0.0
    w: float = 0.0
    w_q: float = 0.0
    lam_eff: float = 0.0


def binom_probs(m: int, mu: float) -> list[float]:
    out = [0.0] * 8
    out[0] = 1.0
    for _ in range(m):
        for k in range(m, 0, -1):
            out[k] = out[k] * (1.0 - mu) + out[k - 1] * mu
        out[0] *= 1.0 - mu
    return out


def solve_dtmc_single(c: int, k_cap: int, lam: float, mu: float) -> DtmcRes:
    p = [[0.0] * 64 for _ in range(64)]
    for n in range(k_cap + 1):
        m = min(n, c)
        bp = binom_probs(m, mu)
        for d in range(m + 1):
            if n < k_cap:
                p[n][n - d] += (1.0 - lam) * bp[d]
                p[n][n + 1 - d] += lam * bp[d]
            else:
                p[n][n - d] += bp[d]
    pi = [0.0] * 64
    for i in range(k_cap + 1):
        pi[i] = 1.0 / (k_cap + 1)
    for _ in range(25000):
        maxd = 0.0
        total = 0.0
        for j in range(k_cap + 1):
            nv = 0.0
            for i in range(k_cap + 1):
                nv += pi[i] * p[i][j]
            df = abs(nv - pi[j])
            if df > maxd:
                maxd = df
            pi[j] = nv
            total += nv
        for j in range(k_cap + 1):
            pi[j] /= total
        if maxd < 1e-14:
            break
    r = DtmcRes(p_block=pi[k_cap], lam_eff=lam * (1.0 - pi[k_cap]))
    for n in range(k_cap + 1):
        m = min(n, c)
        r.l += n * pi[n]
        r.l_q += (n - m) * pi[n]
        r.u += (m / float(c)) * pi[n]
    r.w = r.l / r.lam_eff if r.lam_eff > 0.0 else 0.0
    r.w_q = r.l_q / r.lam_eff if r.lam_eff > 0.0 else 0.0
    return r


def solve_geo_closed_form(lam: float, mu: float) -> DtmcRes:
    l_val = (lam * (1.0 - lam)) / (mu - lam)
    w_val = (1.0 - lam) / (mu - lam)
    wq = w_val - 1.0 / mu
    return DtmcRes(
        p_block=0.0,
        l=l_val,
        l_q=lam * wq,
        u=lam / mu,
        w=w_val,
        w_q=wq,
        lam_eff=lam,
    )


def solve_dtmc_tandem(k1: int, k2: int, lam: float, mu1: float, mu2: float) -> tuple[DtmcRes, DtmcRes, float]:
    s_tot = (k1 + 1) * (k2 + 1)
    p = [[0.0] * 64 for _ in range(64)]
    for n1 in range(k1 + 1):
        for n2 in range(k2 + 1):
            s_idx = n1 * (k2 + 1) + n2
            for a in (0, 1):
                pa = lam if a != 0 else (1.0 - lam)
                for d1 in range(2 if n1 > 0 else 1):
                    pd1 = 1.0 if n1 == 0 else (mu1 if d1 != 0 else (1.0 - mu1))
                    for d2 in range(2 if n2 > 0 else 1):
                        pd2 = 1.0 if n2 == 0 else (mu2 if d2 != 0 else (1.0 - mu2))
                        nn1 = (n1 + a - d1) if n1 < k1 else (n1 - d1)
                        n2_mid = n2 - d2
                        nn2 = (n2_mid + d1) if n2_mid < k2 else n2_mid
                        p[s_idx][nn1 * (k2 + 1) + nn2] += pa * pd1 * pd2
    pi = [0.0] * 64
    for i in range(s_tot):
        pi[i] = 1.0 / s_tot
    for _ in range(25000):
        maxd = 0.0
        total = 0.0
        for j in range(s_tot):
            nv = 0.0
            for i in range(s_tot):
                nv += pi[i] * p[i][j]
            df = abs(nv - pi[j])
            if df > maxd:
                maxd = df
            pi[j] = nv
            total += nv
        for j in range(s_tot):
            pi[j] /= total
        if maxd < 1e-14:
            break
    r1 = DtmcRes()
    r2 = DtmcRes()
    blk1 = 0.0
    blk2 = 0.0
    for n1 in range(k1 + 1):
        for n2 in range(k2 + 1):
            pv = pi[n1 * (k2 + 1) + n2]
            r1.l += n1 * pv
            r2.l += n2 * pv
            if n1 > 0:
                r1.u += pv
                r1.l_q += (n1 - 1) * pv
            if n2 > 0:
                r2.u += pv
                r2.l_q += (n2 - 1) * pv
            if n1 == k1:
                blk1 += pv * lam
            if n1 > 0 and n2 == k2:
                blk2 += pv * mu1 * (1.0 - mu2)
    r1.p_block = blk1
    r1.lam_eff = lam - blk1
    r1.w = r1.l / r1.lam_eff if r1.lam_eff > 0.0 else 0.0
    r1.w_q = r1.l_q / r1.lam_eff if r1.lam_eff > 0.0 else 0.0
    r2.p_block = blk2
    r2.lam_eff = r1.lam_eff - blk2
    r2.w = r2.l / r2.lam_eff if r2.lam_eff > 0.0 else 0.0
    r2.w_q = r2.l_q / r2.lam_eff if r2.lam_eff > 0.0 else 0.0
    return r1, r2, blk2 / (r1.lam_eff if r1.lam_eff > 0.0 else 1.0)
