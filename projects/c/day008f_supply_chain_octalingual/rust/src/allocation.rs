use crate::entities::Sim;

pub fn run_alloc_auction(s: &mut Sim, wh: usize, sku: usize, stores: &[usize], reqs: &[i64], seqs: &[i64], prio: &[i32]) -> Vec<i64> {
    let n = stores.len();
    let mut fills = vec![0i64; n];
    let tot_req: i64 = reqs.iter().sum();
    if tot_req <= 0 { return fills; }
    let avail = s.acc[wh].on_hand[sku];
    s.alloc_cases += 1; s.alloc_total_req += tot_req;
    if avail <= 0 { return fills; }
    if avail >= tot_req {
        fills.copy_from_slice(reqs);
        s.alloc_total_fill += tot_req; s.alloc_nonzero += 1; s.branches[6] += 1;
        return fills;
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| prio[b].cmp(&prio[a]).then_with(|| seqs[a].cmp(&seqs[b])).then_with(|| stores[a].cmp(&stores[b])));
    let mut rem = avail;
    loop {
        let mut max_short = -1i64;
        for &i in &idx { max_short = max_short.max(reqs[i] - fills[i]); }
        if max_short <= 0 || rem <= 0 { break; }
        let mut prog = false;
        for &i in &idx {
            if rem <= 0 { break; }
            if reqs[i] - fills[i] == max_short { fills[i] += 1; rem -= 1; prog = true; }
        }
        if !prog { break; }
    }
    let filled = avail - rem;
    s.alloc_total_fill += filled;
    if filled > 0 { s.alloc_nonzero += 1; }
    s.branches[7] += 1;
    fills
}
