# frozen_string_literal: true

DtmcRes = Struct.new(:p_block, :l, :l_q, :u, :w, :w_q, :lam_eff, keyword_init: true)

def binom_probs(m, mu)
  out = [0.0] * 8
  out[0] = 1.0
  m.times do
    m.downto(1) { |k| out[k] = out[k] * (1.0 - mu) + out[k - 1] * mu }
    out[0] *= (1.0 - mu)
  end
  out
end

def solve_dtmc_single(c, k_cap, lam, mu)
  p = Array.new(64) { [0.0] * 64 }
  (0..k_cap).each do |n|
    m = [n, c].min
    bp = binom_probs(m, mu)
    (0..m).each do |d|
      if n < k_cap
        p[n][n - d] += (1.0 - lam) * bp[d]
        p[n][n + 1 - d] += lam * bp[d]
      else
        p[n][n - d] += bp[d]
      end
    end
  end
  pi = Array.new(64, 0.0)
  (0..k_cap).each { |i| pi[i] = 1.0 / (k_cap + 1) }
  25_000.times do
    maxd = 0.0
    sum = 0.0
    (0..k_cap).each do |j|
      nv = 0.0
      (0..k_cap).each { |i| nv += pi[i] * p[i][j] }
      df = (nv - pi[j]).abs
      maxd = df if df > maxd
      pi[j] = nv
      sum += nv
    end
    (0..k_cap).each { |j| pi[j] /= sum }
    break if maxd < 1e-14
  end
  r = DtmcRes.new(p_block: pi[k_cap], l: 0.0, l_q: 0.0, u: 0.0, w: 0.0, w_q: 0.0, lam_eff: lam * (1.0 - pi[k_cap]))
  (0..k_cap).each do |n|
    m = [n, c].min
    r.l += n * pi[n]
    r.l_q += (n - m) * pi[n]
    r.u += (m.to_f / c) * pi[n]
  end
  r.w = r.lam_eff > 0.0 ? r.l / r.lam_eff : 0.0
  r.w_q = r.lam_eff > 0.0 ? r.l_q / r.lam_eff : 0.0
  r
end

def solve_geo_closed_form(lam, mu)
  l = (lam * (1.0 - lam)) / (mu - lam)
  w = (1.0 - lam) / (mu - lam)
  wq = w - 1.0 / mu
  DtmcRes.new(p_block: 0.0, l: l, l_q: lam * wq, u: lam / mu, w: w, w_q: wq, lam_eff: lam)
end

def solve_dtmc_tandem(k1, k2, lam, mu1, mu2)
  s_tot = (k1 + 1) * (k2 + 1)
  p = Array.new(64) { [0.0] * 64 }
  (0..k1).each do |n1|
    (0..k2).each do |n2|
      s_idx = n1 * (k2 + 1) + n2
      [0, 1].each do |a|
        pa = a != 0 ? lam : (1.0 - lam)
        (0..(n1 > 0 ? 1 : 0)).each do |d1|
          pd1 = n1 == 0 ? 1.0 : (d1 != 0 ? mu1 : (1.0 - mu1))
          (0..(n2 > 0 ? 1 : 0)).each do |d2|
            pd2 = n2 == 0 ? 1.0 : (d2 != 0 ? mu2 : (1.0 - mu2))
            nn1 = n1 < k1 ? (n1 + a - d1) : (n1 - d1)
            n2_mid = n2 - d2
            nn2 = n2_mid < k2 ? (n2_mid + d1) : n2_mid
            p[s_idx][nn1 * (k2 + 1) + nn2] += pa * pd1 * pd2
          end
        end
      end
    end
  end
  pi = Array.new(64, 0.0)
  s_tot.times { |i| pi[i] = 1.0 / s_tot }
  25_000.times do
    maxd = 0.0
    sum = 0.0
    s_tot.times do |j|
      nv = 0.0
      s_tot.times { |i| nv += pi[i] * p[i][j] }
      df = (nv - pi[j]).abs
      maxd = df if df > maxd
      pi[j] = nv
      sum += nv
    end
    s_tot.times { |j| pi[j] /= sum }
    break if maxd < 1e-14
  end
  r1 = DtmcRes.new(p_block: 0.0, l: 0.0, l_q: 0.0, u: 0.0, w: 0.0, w_q: 0.0, lam_eff: 0.0)
  r2 = DtmcRes.new(p_block: 0.0, l: 0.0, l_q: 0.0, u: 0.0, w: 0.0, w_q: 0.0, lam_eff: 0.0)
  blk1 = 0.0
  blk2 = 0.0
  (0..k1).each do |n1|
    (0..k2).each do |n2|
      pv = pi[n1 * (k2 + 1) + n2]
      r1.l += n1 * pv; r2.l += n2 * pv
      if n1 > 0; r1.u += pv; r1.l_q += (n1 - 1) * pv; end
      if n2 > 0; r2.u += pv; r2.l_q += (n2 - 1) * pv; end
      blk1 += pv * lam if n1 == k1
      blk2 += pv * mu1 * (1.0 - mu2) if n1 > 0 && n2 == k2
    end
  end
  r1.p_block = blk1; r1.lam_eff = lam - blk1
  r1.w = r1.lam_eff > 0.0 ? r1.l / r1.lam_eff : 0.0; r1.w_q = r1.lam_eff > 0.0 ? r1.l_q / r1.lam_eff : 0.0
  r2.p_block = blk2; r2.lam_eff = r1.lam_eff - blk2
  r2.w = r2.lam_eff > 0.0 ? r2.l / r2.lam_eff : 0.0; r2.w_q = r2.lam_eff > 0.0 ? r2.l_q / r2.lam_eff : 0.0
  [r1, r2, blk2 / (r1.lam_eff > 0.0 ? r1.lam_eff : 1.0)]
end
