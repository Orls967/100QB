def run_alloc_auction(s, wh, sku, stores, reqs, seqs, prio)
  n = stores.length
  fills = Array.new(n, 0)
  tot_req = reqs.sum
  return fills if tot_req <= 0
  avail = s.acc[wh][:on_hand][sku]
  s.alloc_cases += 1; s.alloc_total_req += tot_req
  return fills if avail <= 0
  if avail >= tot_req
    s.alloc_total_fill += tot_req; s.alloc_nonzero += 1; s.branches[6] += 1
    return reqs.dup
  end
  idx = (0...n).to_a.sort do |a, b|
    if prio[a] != prio[b]
      prio[b] <=> prio[a]
    elsif seqs[a] != seqs[b]
      seqs[a] <=> seqs[b]
    else
      stores[a] <=> stores[b]
    end
  end
  rem = avail
  while rem > 0
    max_short = idx.map { |i| reqs[i] - fills[i] }.max
    break if max_short <= 0
    prog = false
    idx.each do |i|
      break if rem <= 0
      if reqs[i] - fills[i] == max_short
        fills[i] += 1; rem -= 1; prog = true
      end
    end
    break unless prog
  end
  filled = avail - rem
  s.alloc_total_fill += filled
  s.alloc_nonzero += 1 if filled > 0
  s.branches[7] += 1
  fills
end
