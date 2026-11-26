# frozen_string_literal: true
require "fileutils"

MASTER_SEED = 20261126
MASK64 = 0xFFFFFFFFFFFFFFFF
W_WALKS = 4000
W_RUIN = 2000
BENCH_STEPS = 200000
T_EXACT_MAX = 5000
TS_MSD = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096].freeze
TS_RET = [10, 100, 1000, 10000].freeze
TWO_M53 = 2.0**-53

def sm64(z)
  z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
  z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
  (z ^ (z >> 31)) & MASK64
end

class Rng
  def initialize(seed)
    st = seed & MASK64
    @s = Array.new(4) do
      st = (st + 0x9E3779B97F4A7C15) & MASK64
      sm64(st)
    end
  end

  def u64
    s0, s1, s2, s3 = @s
    x = (s0 + s3) & MASK64
    res = (((x << 23) | (x >> 41)) & MASK64) + s0
    res &= MASK64
    t = (s1 << 17) & MASK64
    s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t
    s3 = ((s3 << 45) | (s3 >> 19)) & MASK64
    @s[0] = s0; @s[1] = s1; @s[2] = s2; @s[3] = s3
    res
  end

  def u01; (u64 >> 11) * TWO_M53; end

  def rint(n)
    k = (u01 * n).floor
    k >= n ? n - 1 : k
  end
end

def task_seed(d, st, idx)
  h = 0x517CC1B727220A95
  [d, st, idx].each { |v| h = sm64(h ^ (v & MASK64)) }
  sm64(MASTER_SEED ^ h)
end

def fnv1a_i32(coords)
  h = 0xCBF29CE484222325
  coords.pack("l<*").each_byte { |b| h = ((h ^ b) * 0x100000001B3) & MASK64 }
  h
end

FileUtils.mkdir_p("out/ruby/figures")
$msd = Array.new(3) { Array.new(13, 0.0) }
$ret_sim = Array.new(3) { Array.new(4, 0.0) }
$ret_se = Array.new(3) { Array.new(4, 0.0) }
$ret_ex = Array.new(3) { Array.new(4, 0.0) }
$ruin_p = Array.new(6, 0.0)
$ruin_pe = Array.new(6, 0.0)
$s03_1024 = Array.new(3) { Array.new(W_WALKS, 0.0) }
$s04_10k = Array.new(3) { Array.new(W_WALKS, 0) }
$stage_sec = Array.new(12, 0.0)

def s01
  r = Rng.new(MASTER_SEED)
  x0, x1, x2 = r.u64, r.u64, r.u64
  ok = (x0 == 0xA9920E0E7DAB3EF0 && x1 == 0x033B9CEB4B3BA7A0 && x2 == 0xE7473A6A39AFCEA6)
  w = [1.0, 1.0 / 3.0, 0.0, 0.0, 0.0]
  (1..3).each { |n| w[n + 1] = ((10.0 * n * n + 10.0 * n + 3.0) * w[n] - n * n * w[n - 1]) / (9.0 * (n + 1) * (n + 1)) }
  eok = (w[2] - 15.0 / 81.0).abs < 1e-14 && (w[3] - 93.0 / 729.0).abs < 1e-14 && (w[4] - 639.0 / 6561.0).abs < 1e-14
  File.write("out/ruby/s01_selftest.txt", format("ok,x0,x1,x2,enum_ok\n%d,%016x,%016x,%016x,%d\n", (ok && eok) ? 1 : 0, x0, x1, x2, eok ? 1 : 0))
  exit(1) unless ok && eok
end

def s02
  File.open("out/ruby/s02_parity.csv", "w") do |f|
    f.puts "d,walk,x0,x1,x2,fnv1a_hex"
    (1..3).each do |d|
      coords = []; pos_all = []
      8.times do |w|
        r = Rng.new(task_seed(d, 2, w)); p = [0, 0, 0]
        200.times { k = r.rint(2 * d); p[k >> 1] += (k & 1) != 0 ? -1 : 1 }
        pos_all << p; d.times { |a| coords << p[a] }
      end
      h = fnv1a_i32(coords)
      8.times { |w| f.puts format("%d,%d,%d,%d,%d,%016x", d, w, pos_all[w][0], pos_all[w][1], pos_all[w][2], h) }
    end
  end
end

def s03
  File.open("out/ruby/s03_msd.csv", "w") do |f|
    f.puts "d,t,msd,se,exact"
    (1..3).each do |d|
      s1 = Array.new(13, 0.0); s2 = Array.new(13, 0.0)
      W_WALKS.times do |w|
        r = Rng.new(task_seed(d, 3, w)); p = [0, 0, 0]; ti = 0
        (1..4096).each do |t|
          k = r.rint(2 * d); p[k >> 1] += (k & 1) != 0 ? -1 : 1
          next unless t == TS_MSD[ti]
          r2 = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).to_f
          s1[ti] += r2; s2[ti] += r2 * r2
          $s03_1024[d - 1][w] = r2 if ti == 10
          ti += 1
        end
      end
      13.times do |ti|
        m = s1[ti] / W_WALKS; se = Math.sqrt([0.0, s2[ti] / W_WALKS - m * m].max / W_WALKS)
        $msd[d - 1][ti] = m
        f.puts format("%d,%d,%.8f,%.8f,%d", d, TS_MSD[ti], m, se, TS_MSD[ti])
      end
    end
  end
end

def s04
  File.open("out/ruby/s04_return_sim.csv", "w") do |f|
    f.puts "d,T,p_hat,se"
    (1..3).each do |d|
      hits = [0, 0, 0, 0]
      W_WALKS.times do |w|
        r = Rng.new(task_seed(d, 4, w)); p = [0, 0, 0]; ret_t = 0
        (1..10000).each do |t|
          k = r.rint(2 * d); p[k >> 1] += (k & 1) != 0 ? -1 : 1
          if (p[0] | p[1] | p[2]) == 0; ret_t = t; break; end
        end
        if ret_t > 0
          $s04_10k[d - 1][w] = 1
          4.times { |ti| hits[ti] += 1 if ret_t <= TS_RET[ti] }
        end
      end
      4.times do |ti|
        ph = hits[ti].to_f / W_WALKS; se = Math.sqrt([0.0, ph * (1.0 - ph)].max / W_WALKS)
        $ret_sim[d - 1][ti] = ph; $ret_se[d - 1][ti] = se
        f.puts format("%d,%d,%.8f,%.8f", d, TS_RET[ti], ph, se)
      end
    end
  end
end

def s05
  u1 = Array.new(T_EXACT_MAX + 1, 1.0); u2 = Array.new(T_EXACT_MAX + 1, 1.0)
  w3 = Array.new(T_EXACT_MAX + 1, 1.0); u3 = Array.new(T_EXACT_MAX + 1, 1.0); fv = Array.new(T_EXACT_MAX + 1, 0.0)
  w3[1] = 1.0 / 3.0
  (1...T_EXACT_MAX).each { |n| w3[n + 1] = ((10.0 * n * n + 10.0 * n + 3.0) * w3[n] - n * n * w3[n - 1]) / (9.0 * (n + 1) * (n + 1)) }
  (1..T_EXACT_MAX).each { |n| u1[n] = u1[n - 1] * (2.0 * n - 1.0) / (2.0 * n); u2[n] = u1[n] * u1[n]; u3[n] = u1[n] * w3[n] }
  us = [u1, u2, u3]
  File.open("out/ruby/s05_return_exact.csv", "w") do |f|
    f.puts "d,T,F_exact"
    (1..3).each do |d|
      u = us[d - 1]; acc = 0.0
      (1..T_EXACT_MAX).each do |n|
        s = u[n]; k = 1
        while k < n; s -= fv[k] * u[n - k]; k += 1; end
        fv[n] = s; acc += s
        4.times { |ti| if 2 * n == TS_RET[ti]; $ret_ex[d - 1][ti] = acc; f.puts format("%d,%d,%.10f", d, TS_RET[ti], acc); end }
      end
    end
  end
end

def s06
  File.open("out/ruby/s06_return_compare.csv", "w") do |f|
    f.puts "d,T,p_hat,se,F_exact,z_score,p_inf_est"
    p3_inf = $ret_sim[2][3] + ($ret_sim[2][3] - $ret_sim[2][2]) / (Math.sqrt(10.0) - 1.0)
    (1..3).each do |d|
      4.times do |ti|
        z = $ret_se[d - 1][ti] > 1e-12 ? ($ret_sim[d - 1][ti] - $ret_ex[d - 1][ti]) / $ret_se[d - 1][ti] : 0.0
        f.puts format("%d,%d,%.8f,%.8f,%.10f,%.4f,%.8f", d, TS_RET[ti], $ret_sim[d - 1][ti], $ret_se[d - 1][ti], $ret_ex[d - 1][ti], z, (d == 3 && ti == 3) ? p3_inf : 0.0)
      end
    end
  end
end

def s07
  ns = [10, 10, 10, 50, 50, 50]; xs = [2, 5, 8, 10, 25, 40]
  File.open("out/ruby/s07_ruin.csv", "w") do |f|
    f.puts "N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z"
    6.times do |ci|
      n, x0 = ns[ci], xs[ci]; wins = 0; st1 = 0.0; st2 = 0.0
      W_RUIN.times do |w|
        r = Rng.new(task_seed(1, 7, ci * 100000 + w)); x = x0; steps = 0
        while x > 0 && x < n; x += r.rint(2) == 0 ? 1 : -1; steps += 1; end
        wins += 1 if x == n; st1 += steps; st2 += steps * steps
      end
      ph = wins.to_f / W_RUIN; pse = Math.sqrt([0.0, ph * (1.0 - ph)].max / W_RUIN); pe = x0.to_f / n
      th = st1 / W_RUIN; tse = Math.sqrt([0.0, st2 / W_RUIN - th * th].max / W_RUIN); te = (x0 * (n - x0)).to_f
      $ruin_p[ci] = ph; $ruin_pe[ci] = pe
      f.puts format("%d,%d,%.8f,%.8f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f", n, x0, ph, pse, pe, (ph - pe) / pse, th, tse, te, (th - te) / tse)
    end
  end
end

def s08
  File.open("out/ruby/s08_endpoint.csv", "w") do |f|
    f.puts "d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count"
    hist = Array.new(40, 0); n_ep = 4000
    (1..3).each do |d|
      m1 = m2 = m4 = 0.0
      n_ep.times do |w|
        r = Rng.new(task_seed(d, 8, w)); p = [0, 0, 0]
        1024.times { k = r.rint(2 * d); p[k >> 1] += (k & 1) != 0 ? -1 : 1 }
        x = p[0].to_f; m1 += x; m2 += x * x; m4 += x * x * x * x
        if d == 1; b = [[((x + 128.0) / 6.4).floor, 0].max, 39].min; hist[b] += 1; end
      end
      m1 /= n_ep; v = m2 / n_ep - m1 * m1
      f.puts format("%d,%.6f,%.6f,%.6f,%.6f,-1,0,0,%d", d, v, 1024.0 / d, (m4 / n_ep) / (v * v), 3.0 + (d - 3.0) / 1024.0, n_ep)
    end
    40.times { |b| f.puts format("1,0,0,0,0,%d,%.2f,%.2f,%d", b, -128.0 + b * 6.4, -128.0 + (b + 1) * 6.4, hist[b]) }
  end
end

def s09
  File.open("out/ruby/s09_ci.csv", "w") do |f|
    f.puts "metric,d,estimate,boot_se,ci_lo,ci_hi"
    (1..3).each do |d|
      r1 = Rng.new(task_seed(d, 9, 0))
      boot = Array.new(200) { s = 0.0; W_WALKS.times { s += $s03_1024[d - 1][r1.rint(W_WALKS)] }; s / W_WALKS }.sort
      bm = boot.sum / 200.0; bv = boot.sum { |v| (v - bm)**2 } / 199.0
      f.puts format("msd_1024,%d,%.6f,%.6f,%.6f,%.6f", d, $msd[d - 1][10], Math.sqrt(bv), boot[5], boot[194])
      r2 = Rng.new(task_seed(d, 9, 1))
      boot2 = Array.new(200) { s = 0; W_WALKS.times { s += $s04_10k[d - 1][r2.rint(W_WALKS)] }; s.to_f / W_WALKS }.sort
      bm2 = boot2.sum / 200.0; bv2 = boot2.sum { |v| (v - bm2)**2 } / 199.0
      f.puts format("return_10000,%d,%.6f,%.6f,%.6f,%.6f", d, $ret_sim[d - 1][3], Math.sqrt(bv2), boot2[5], boot2[194])
    end
  end
end

def s10
  File.open("out/ruby/s10_fit.csv", "w") do |f|
    f.puts "d,slope,slope_se,intercept"
    (1..3).each do |d|
      xs = TS_MSD.map { |t| Math.log(t) }; ys = $msd[d - 1].map { |m| Math.log(m) }
      mx = xs.sum / 13.0; my = ys.sum / 13.0
      num = (0...13).sum { |i| (xs[i] - mx) * (ys[i] - my) }; den = xs.sum { |x| (x - mx)**2 }
      sl = num / den; ic = my - sl * mx; rss = (0...13).sum { |i| (ys[i] - (sl * xs[i] + ic))**2 }
      f.puts format("%d,%.6f,%.6f,%.6f", d, sl, Math.sqrt((rss / 11.0) / den), ic)
    end
  end
end

def s11
  File.open("out/ruby/s11_bench.csv", "w") do |f|
    f.puts "d,steps_per_rep,reps,ns_med,ns_mad"
    (1..3).each do |d|
      vals = Array.new(7) do |rep|
        r = Rng.new(task_seed(d, 11, rep)); p = [0, 0, 0]
        t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
        BENCH_STEPS.times { k = r.rint(2 * d); p[k >> 1] += (k & 1) != 0 ? -1 : 1 }
        (Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0) * 1e9 / BENCH_STEPS
      end.sort
      med = vals[3]; dev = vals.map { |v| (v - med).abs }.sort
      f.puts format("%d,%d,7,%.4f,%.4f", d, BENCH_STEPS, med, dev[3])
    end
  end
end

def s12
  cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
  s1 = +"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>"
  3.times do |d|
    pts = (0...13).map { |i| format("%.1f,%.1f", 60.0 + 500.0 * Math.log10(TS_MSD[i]) / Math.log10(4096), 310.0 - 280.0 * Math.log10([$msd[d][i], 1.0].max) / Math.log10(4096)) }.join(" ")
    s1 << "<polyline fill=\"none\" stroke=\"#{cols[d]}\" stroke-width=\"1.6\" points=\"#{pts}\"/>"
  end
  File.write("out/ruby/figures/msd.svg", s1 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby MSD Log-Log (y=t)</text></svg>\n")

  s2 = +"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>"
  3.times do |d|
    pts = (0...4).map { |i| format("%.1f,%.1f", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * $ret_ex[d][i]) }.join(" ")
    s2 << "<polyline fill=\"none\" stroke=\"#{cols[d]}\" stroke-width=\"1.6\" points=\"#{pts}\"/>"
    4.times { |i| s2 << format("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4\" fill=\"%s\"/>", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * $ret_sim[d][i], cols[d]) }
  end
  File.write("out/ruby/figures/return_prob.svg", s2 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby Return Prob F_d(T) Exact vs Sim</text></svg>\n")

  s3 = +"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>"
  6.times { |i| s3 << format("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"5\" fill=\"#1f77b4\"/>", 60.0 + 500.0 * $ruin_pe[i], 310.0 - 280.0 * $ruin_p[i]) }
  File.write("out/ruby/figures/ruin.svg", s3 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby Gambler's Ruin Exact vs Sim</text></svg>\n")

  File.write("out/ruby/s12_summary.json", format("{\"lang\":\"ruby\",\"W\":%d,\"ruin_W\":%d,\"bench_steps\":%d,\"toolchain\":\"ruby 4.0.9\",\"total_sec\":%.3f}\n", W_WALKS, W_RUIN, BENCH_STEPS, $stage_sec.sum))
end

cmd = ARGV[0] || "all"
stages = %w[s01 s02 s03 s04 s05 s06 s07 s08 s09 s10 s11 s12]
if cmd == "all"
  stages.each_with_index { |nm, i| puts "[ruby] tahap #{i + 1}/12 (#{nm})"; t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC); send(nm); $stage_sec[i] = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0 }
elsif stages.include?(cmd)
  send(cmd)
else
  exit(1)
end
