def run_alloc_auction(s, wh, sku, stores, reqs, seqs, prio):
    n = len(stores)
    fills = [0] * n
    tot_req = sum(reqs)
    if tot_req <= 0: return fills
    avail = s.acc[wh]["on_hand"][sku]
    s.alloc_cases += 1; s.alloc_total_req += tot_req
    if avail <= 0: return fills
    if avail >= tot_req:
        s.alloc_total_fill += tot_req; s.alloc_nonzero += 1; s.branches[6] += 1
        return list(reqs)
    idx = sorted(range(n), key=lambda i: (-prio[i], seqs[i], stores[i]))
    rem = avail
    while rem > 0:
        max_short = max(reqs[i] - fills[i] for i in idx)
        if max_short <= 0: break
        prog = False
        for i in idx:
            if rem <= 0: break
            if reqs[i] - fills[i] == max_short:
                fills[i] += 1; rem -= 1; prog = True
        if not prog: break
    filled = avail - rem
    s.alloc_total_fill += filled
    if filled > 0: s.alloc_nonzero += 1
    s.branches[7] += 1
    return fills
