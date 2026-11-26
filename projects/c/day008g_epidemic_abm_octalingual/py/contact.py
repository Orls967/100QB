from population import BP, ST_S, ST_IA, ST_IS


def vax_eff_bp(s, p, day: int) -> int:
    if p.vax_doses <= 0 or p.vax_day < 0:
        return 0
    if s.rp.vax_perfect != 0:
        return BP
    base = s.cfg["vax_eff1_bp"] if p.vax_doses == 1 else s.cfg["vax_eff2_bp"]
    elapsed = day - p.vax_day
    return base if elapsed <= 30 else max(base - (elapsed - 30) * 35, 2000)


def try_infect(s, src: int, tgt: int, lay: int, day: int, beta_eff: int):
    pt = s.pop[tgt]
    if pt.status != ST_S or s.seen_day[tgt] == day:
        return
    if pt.quar_until >= day:
        s.branches[5] += 1
        return
    p_bp = beta_eff
    if s.rp.interv_on != 0:
        m_day = s.rp.mask_ovr if s.rp.mask_ovr >= 0 else s.cfg["mask_day"]
        if day >= m_day:
            red = (s.cfg["mask_eff_bp"] * s.pop[src].compliance) // BP
            next_bp = (p_bp * (BP - red)) // BP
            if next_bp < p_bp:
                s.branches[4] += 1
            p_bp = next_bp
    ve = vax_eff_bp(s, pt, day)
    if ve > 0:
        p_bp = (p_bp * (BP - ve)) // BP
        if p_bp <= 0:
            s.branches[3] += 1
            return
    if s.rng[1].rand_int(1, BP) <= p_bp:
        s.seen_day[tgt] = day
        s.cand_tgt.append(tgt)
        s.cand_src.append(src)
        s.cand_lay.append(lay)
        s.cand_var.append(s.pop[src].variant)
    elif ve > 0:
        s.branches[3] += 1


def sim_contact_step(s, day: int):
    s.cand_tgt.clear()
    s.cand_src.clear()
    s.cand_lay.clear()
    s.cand_var.clear()
    if s.rp.contact_scale_bp == 0:
        return
    n = s.n
    sc_day = s.rp.school_ovr if s.rp.school_ovr >= 0 else s.cfg["school_close_day"]
    lk_day = s.rp.lock_ovr if s.rp.lock_ovr >= 0 else s.cfg["lockdown_day"]
    school_closed = s.rp.interv_on != 0 and day >= sc_day
    locked = s.rp.interv_on != 0 and day >= lk_day

    for i in range(n):
        pi = s.pop[i]
        st = pi.status
        if st != ST_IA and st != ST_IS:
            continue
        base_bp = s.rp.beta_bp
        if st == ST_IA:
            base_bp = (base_bp * s.cfg["asymp_rel_bp"]) // BP
        if pi.variant == 1:
            base_bp = (base_bp * s.cfg["variant_mult_bp"]) // BP
        if base_bp > BP:
            base_bp = BP

        h = pi.hh_id
        for idx in range(s.hh_off[h], s.hh_off[h + 1]):
            j = s.hh_mem[idx]
            if j != i:
                try_infect(s, i, j, 0, day, base_bp)
        if pi.quar_until >= day:
            s.branches[5] += 1
            continue

        if pi.sw_id >= 0 and not (pi.age_grp == 0 and school_closed):
            sw = pi.sw_id
            lo, hi = s.sw_off[sw], s.sw_off[sw + 1]
            sz = hi - lo
            if sz > 1:
                n_sw_c = 6 if pi.age_grp == 0 else 5
                if locked:
                    n_sw_c = (n_sw_c + 1) // 2
                for _ in range(n_sw_c):
                    j = s.sw_mem[lo + s.rng[1].rand_int(0, sz - 1)]
                    if j != i:
                        try_infect(s, i, j, 1, day, base_bp)

        n_comm = 12 if pi.super_flag != 0 else 3
        if locked:
            n_comm = (n_comm + 1) // 3
        for _ in range(n_comm):
            j = s.rng[1].rand_int(0, n - 1)
            if j != i:
                try_infect(s, i, j, 2, day, base_bp)
