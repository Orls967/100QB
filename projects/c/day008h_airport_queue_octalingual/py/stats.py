import math
import struct
from config import N_STN, SLOTS_PER_DAY, STN_SECURITY, default_rp
from dtmc import solve_dtmc_single, solve_dtmc_tandem, solve_geo_closed_form
from engine import sim_run
from network import Sim
from rng import FNV_OFFSET, MASK64, Xoshiro256, fnv1a_bytes, fnv1a_i64, fnv1a_u32pair, splitmix64


def run_single_sim(rng: Xoshiro256, c: int, k_cap: int, lam_bp: int, mu_bp: int, slots: int) -> tuple[float, float, float]:
    q = 0
    arr = 0
    blk = 0
    sum_l = 0
    sum_b = 0
    wstart = slots // 5
    for t in range(slots):
        a = rng.bernoulli(lam_bp)
        m = min(q, c)
        d = 0
        for _ in range(m):
            if rng.bernoulli(mu_bp) != 0:
                d += 1
        if q < k_cap:
            q += a
            if a != 0:
                arr += 1
        elif a != 0:
            blk += 1
        q -= d
        if t >= wstart:
            mb = min(q, c)
            sum_l += q
            sum_b += mb
    ns = float(slots - wstart)
    return sum_l / ns, sum_b / (ns * c), blk / float(arr + blk if arr + blk > 0 else 1)


def run_tandem_sim(rng: Xoshiro256, k1: int, k2: int, lam_bp: int, mu1_bp: int, mu2_bp: int, slots: int) -> tuple[float, float, float, float]:
    q1 = 0
    q2 = 0
    sl1 = 0
    sl2 = 0
    su1 = 0
    su2 = 0
    wstart = slots // 5
    for t in range(slots):
        a = rng.bernoulli(lam_bp)
        d1 = rng.bernoulli(mu1_bp) if q1 > 0 else 0
        d2 = rng.bernoulli(mu2_bp) if q2 > 0 else 0
        if q1 < k1:
            q1 += a
        q1 -= d1
        q2 -= d2
        if d1 != 0 and q2 < k2:
            q2 += 1
        if t >= wstart:
            sl1 += q1
            sl2 += q2
            if q1 > 0:
                su1 += 1
            if q2 > 0:
                su2 += 1
    ns = float(slots - wstart)
    return sl1 / ns, sl2 / ns, su1 / ns, su2 / ns


def write_s01_to_s13(c: dict, out_dir: str) -> None:
    hrng = FNV_OFFSET
    pairs = [(0, 0), (1, 2), (2, 1)]
    exp_hv = [0xBE40B757FE0B6697, 0x659394D6F2A43A95, 0x3ABDA5BAC0DD49CB]
    rng_ok = True
    for i in range(3):
        x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1])
        _, sd = splitmix64(x)
        st = Xoshiro256(sd)
        vals = [st.next_u64() for _ in range(4)]
        hv = fnv1a_bytes(struct.pack("<4Q", *vals))
        if hv != exp_hv[i]:
            rng_ok = False
        hrng = fnv1a_i64(hrng, hv)

    s0 = Sim(c, default_rp(c, 1))
    sim_run(s0)
    inv_ok = s0.inv_violations == 0
    with open(f"{out_dir}/s01_selftest.txt", "w", encoding="utf-8") as f:
        f.write(f"ok,h_rng,inv_ok,tot_arr,tot_done\n{1 if (rng_ok and inv_ok) else 0},{hrng:016x},{1 if inv_ok else 0},{s0.total_arr},{s0.total_done}\n")
    with open(f"{out_dir}/s02_config.csv", "w", encoding="utf-8") as f:
        f.write(
            f"key,value\nmaster_seed,{c['master_seed']}\nslots_per_day,{c['slots_per_day']}\n"
            f"parity_days,{c['parity_days']}\narr_scale_bp,{c['arr_scale_bp']}\nbelt_cap,{c['belt_cap']}\n"
        )
    with open(f"{out_dir}/s03_flights.csv", "w", encoding="utf-8") as f:
        f.write("flight_id,dep_slot,gate_close_slot,seats\n")
        for h in range(24):
            f.write(f"{h},{h * 360 + 300},{h * 360 + 300 - c['gate_close_lead']},180\n")
    with open(f"{out_dir}/s03_arrival_profile.csv", "w", encoding="utf-8") as f:
        f.write("hour,arr_bp\n")
        for h in range(24):
            f.write(f"{h},{s0.arr_prof[h]}\n")
    snm = ["CHECKIN", "SECURITY", "IMIGRASI", "GATE", "BAGGAGE_BELT", "CUSTOMS_ARRIVAL"]
    with open(f"{out_dir}/s04_stations.csv", "w", encoding="utf-8") as f:
        f.write("stn_id,name,c_max,mu_bp,cap,has_prio\n")
        for i in range(N_STN):
            f.write(f"{i},{snm[i]},{s0.stn[i].c_max},{c['mu_bp'][i]},{s0.stn[i].cap},{s0.stn[i].has_prio}\n")
    with open(f"{out_dir}/s04_roster.csv", "w", encoding="utf-8") as f:
        f.write("hour,checkin,security,imigrasi,gate,belt,customs\n")
        for h in range(24):
            f.write(f"{h}," + ",".join(str(v) for v in s0.roster[h]) + "\n")

    with open(f"{out_dir}/s05_queue_script.csv", "w", encoding="utf-8") as f:
        f.write("case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost\n")
        for k in range(12):
            q0 = (k % 4) + 1
            cap = 3 if k >= 6 else 10
            brk = 1 if (k % 3 == 2) else 0
            gc = 1 if k == 11 else 0
            acc = min(q0, cap)
            srv = 0 if gc != 0 else (acc - 1 if brk != 0 else acc)
            lost = acc if gc != 0 else (q0 - acc)
            f.write(f"{k},{k % 3},{q0},{cap},{k % 2},{brk},{gc},{acc},{srv},{lost}\n")

    with open(f"{out_dir}/s06_bag_script.csv", "w", encoding="utf-8") as f:
        f.write("case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed\n")
        for k in range(12):
            nb = k % 3
            bl = (k // 3) * 4
            cap = 10
            acc = nb if bl + nb <= cap else max(cap - bl, 0)
            f.write(f"{k},{nb},{bl},{cap},{c['bag_miss_bp']},{acc},{nb - acc},{1 if (k % 4 == 3 and acc > 0) else 0}\n")

    lam = 0.08
    mu = 0.12
    with open(f"{out_dir}/s07_dtmc_single.csv", "w", encoding="utf-8") as f:
        f.write("c,K,lam,mu,p_block,L,L_q,U,W,W_q\n")
        for cc in (1, 2, 3):
            for kk in (5, 10, 60):
                dr = solve_dtmc_single(cc, kk, lam, mu)
                f.write(f"{cc},{kk},{lam:.4f},{mu:.4f},{dr.p_block:.12f},{dr.l:.12f},{dr.l_q:.12f},{dr.u:.12f},{dr.w:.12f},{dr.w_q:.12f}\n")

    t1, t2, pbt = solve_dtmc_tandem(6, 6, 0.08, 0.14, 0.12)
    with open(f"{out_dir}/s08_dtmc_tandem.csv", "w", encoding="utf-8") as f:
        f.write(
            "K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n"
            f"6,6,0.0800,0.1400,0.1200,{t1.p_block:.12f},{t1.l:.12f},{t1.u:.12f},{t2.p_block:.12f},{t2.l:.12f},{t2.u:.12f},{pbt:.12f}\n"
        )
    cf = solve_geo_closed_form(lam, mu)
    d60 = solve_dtmc_single(1, 60, lam, mu)
    with open(f"{out_dir}/s09_geo_closed_form.csv", "w", encoding="utf-8") as f:
        f.write(
            "lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n"
            f"{lam:.4f},{mu:.4f},{cf.l:.12f},{cf.w:.12f},{cf.u:.12f},{d60.l:.12f},{abs(cf.l - d60.l) / cf.l:.12e}\n"
        )

    val_slots = 6000
    with open(f"{out_dir}/s10_single_validation.csv", "w", encoding="utf-8") as f:
        f.write("c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se\n")
        for ic, cc in enumerate((1, 2, 3)):
            kk = 10 if cc == 1 else 5
            dr = solve_dtmc_single(cc, kk, lam, mu)
            sl = sl2 = su = su2 = sp = sp2 = 0.0
            for r in range(50):
                x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(500 + ic, r)
                _, sd = splitmix64(x)
                rng = Xoshiro256(sd)
                vl, vu, vp = run_single_sim(rng, cc, kk, 800, 1200, val_slots)
                sl += vl
                sl2 += vl * vl
                su += vu
                su2 += vu * vu
                sp += vp
                sp2 += vp * vp
            for nm, dt, s1, s2 in (("L", dr.l, sl, sl2), ("U", dr.u, su, su2), ("p_block", dr.p_block, sp, sp2)):
                mn = s1 / 50.0
                se = math.sqrt(max((s2 / 50.0 - mn * mn) / 49.0, 1e-18))
                z = (mn - dt) / se
                f.write(f"{cc},{kk},{nm},{dt:.6f},{mn:.6f},{se:.6f},{z:.3f},{1 if abs(z) <= 3.0 else 0}\n")

    with open(f"{out_dir}/s11_tandem_validation.csv", "w", encoding="utf-8") as f:
        f.write("metric,dtmc,sim_mean,sim_se,z_score,within_3se\n")
        sl1 = sl12 = sl2 = sl22 = su1 = su12 = su2 = su22 = 0.0
        for r in range(50):
            x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(600, r)
            _, sd = splitmix64(x)
            rng = Xoshiro256(sd)
            l1, l2, u1, u2 = run_tandem_sim(rng, 6, 6, 800, 1400, 1200, val_slots)
            sl1 += l1
            sl12 += l1 * l1
            sl2 += l2
            sl22 += l2 * l2
            su1 += u1
            su12 += u1 * u1
            su2 += u2
            su22 += u2 * u2
        for nm, dt, s1, s2 in (("L1", t1.l, sl1, sl12), ("L2", t2.l, sl2, sl22), ("U1", t1.u, su1, su12), ("U2", t2.u, su2, su22)):
            mn = s1 / 50.0
            se = math.sqrt(max((s2 / 50.0 - mn * mn) / 49.0, 1e-18))
            z = (mn - dt) / se
            f.write(f"{nm},{dt:.6f},{mn:.6f},{se:.6f},{z:.3f},{1 if abs(z) <= 3.0 else 0}\n")

    rpw = default_rp(c, 2)
    rpw.audit_on = 1
    sw = Sim(c, rpw)
    sim_run(sw)
    n_blk = sw.audit_slots // 5
    blk_y = [sum(float(sw.q_series[(b * 5 + k) * N_STN + STN_SECURITY]) for k in range(5)) / 5.0 for b in range(n_blk)]
    best_d = 0
    best_mser = 1e30
    for d in range(0, n_blk // 2, 10):
        m = float(n_blk - d)
        mn = sum(blk_y[d:]) / m
        s2 = sum((y - mn) * (y - mn) for y in blk_y[d:])
        mser = s2 / (m * m)
        if mser < best_mser:
            best_mser = mser
            best_d = d * 5
    with open(f"{out_dir}/s12_warmup.csv", "w", encoding="utf-8") as f:
        f.write(f"warmup_slots,mser_min\n{best_d},{best_mser:.6f}\n")

    with open(f"{out_dir}/s13_calibration.csv", "w", encoding="utf-8") as f:
        f.write("arr_scale_bp,max_util,loss_frac,calibrated\n")
        for scale in (15000, 18000, 22000, 26000, 30000):
            rpc = default_rp(c, 1)
            rpc.arr_scale_bp = scale
            sc = Sim(c, rpc)
            sim_run(sc)
            max_u = max(sc.stn[st].busy_slots / float(SLOTS_PER_DAY * sc.stn[st].c_max) for st in range(4))
            lf = (sc.total_lost / float(sc.total_arr)) if sc.total_arr > 0 else 0.0
            cal = 1 if (0.55 <= max_u <= 0.92 and lf < 0.15) else 0
            f.write(f"{scale},{max_u:.4f},{lf:.4f},{cal}\n")
