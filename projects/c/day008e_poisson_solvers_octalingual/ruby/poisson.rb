# frozen_string_literal: true
require "fileutils"

PI_VAL = 3.14159265358979323846
NS_ALL = [15, 31, 63, 127].freeze
MNAMES = %w[jacobi gs rbgs sor cg pcg_ssor mg].freeze
SIN_C = [1.0, -1.0/6.0, 1.0/120.0, -1.0/5040.0, 1.0/362880.0, -1.0/39916800.0, 1.0/6227020800.0, -1.0/1307674368000.0, 1.0/355687428096000.0, -1.0/121645100408832000.0].freeze

def sin_spec(x)
  sgn = 1.0
  if x < 0.0; x = -x; sgn = -sgn; end
  if x > PI_VAL; x -= PI_VAL; sgn = -sgn; end
  x = PI_VAL - x if x > 0.5 * PI_VAL
  x2 = x * x; p = SIN_C[9]
  8.downto(0) { |k| p = p * x2 + SIN_C[k] }
  sgn * (x * p)
end
def cos_spec(x); sin_spec(0.5 * PI_VAL - x); end
def omega_opt(n)
  { 15 => 1.673513677715992, 31 => 1.8214651907890225, 63 => 1.906454701582762, 127 => 1.952093233850055 }[n]
end
def fnv1a_f64(v)
  h = 0xCBF29CE484222325; mask = 0xFFFFFFFFFFFFFFFF
  v.pack("E*").each_byte { |b| h = ((h ^ b) * 0x100000001B3) & mask }
  h
end
def apply_a(n, u, out)
  n.times do |i|
    r = i * n
    n.times do |j|
      k = r + j
      w = j > 0 ? u[k - 1] : 0.0; e = j + 1 < n ? u[k + 1] : 0.0
      s = i > 0 ? u[k - n] : 0.0; nn = i + 1 < n ? u[k + n] : 0.0
      out[k] = (((4.0 * u[k] - w) - e) - s) - nn
    end
  end
end
def dot_seq(a, b)
  s = 0.0; i = 0; len = a.length
  while i < len; s += a[i] * b[i]; i += 1; end
  s
end
def norm2_seq(a); Math.sqrt(dot_seq(a, a)); end
def max_abs_diff(a, b)
  m = 0.0; i = 0; len = a.length
  while i < len; d = (a[i] - b[i]).abs; m = d if d > m; i += 1; end
  m
end
def build_prob(n)
  h = 1.0 / (n + 1); h2 = h * h
  sx = Array.new(n) { |i| sin_spec(PI_VAL * (i + 1) * h) }
  sy = Array.new(n) { |i| sin_spec(2.0 * PI_VAL * (i + 1) * h) }
  b = Array.new(n * n, 0.0); ue = Array.new(n * n, 0.0)
  n.times do |i|
    x = (i + 1) * h; xx = x * (1.0 - x)
    n.times do |j|
      y = (j + 1) * h; yy = y * (1.0 - y); sxy = sx[i] * sy[j]
      ue[i * n + j] = sxy + xx * yy
      b[i * n + j] = h2 * (5.0 * PI_VAL * PI_VAL * sxy + 2.0 * xx + 2.0 * yy)
    end
  end
  [b, ue]
end
def compute_res(n, u, b, r)
  apply_a(n, u, r); i = 0; n2 = n * n
  while i < n2; r[i] = b[i] - r[i]; i += 1; end
  norm2_seq(r)
end
def rbgs_sweep(n, u, b)
  2.times do |col|
    n.times do |i|
      row = i * n
      n.times do |j|
        next if ((i + j) & 1) != col
        k = row + j
        w = j > 0 ? u[k - 1] : 0.0; e = j + 1 < n ? u[k + 1] : 0.0
        s = i > 0 ? u[k - n] : 0.0; nn = i + 1 < n ? u[k + n] : 0.0
        u[k] = 0.25 * ((((b[k] + w) + e) + s) + nn)
      end
    end
  end
end
def restrict_fw(nf, rf)
  nc = (nf - 1) / 2; rc = Array.new(nc * nc, 0.0)
  nc.times do |ic|
    fi = 2 * ic + 1
    nc.times do |jc|
      fj = 2 * jc + 1; k = fi * nf + fj
      edge = ((rf[k - 1] + rf[k + 1]) + rf[k - nf]) + rf[k + nf]
      corn = ((rf[k - nf - 1] + rf[k - nf + 1]) + rf[k + nf - 1]) + rf[k + nf + 1]
      rc[ic * nc + jc] = 4.0 * (0.25 * rf[k] + 0.125 * edge + 0.0625 * corn)
    end
  end
  rc
end
def prolong_add(nf, uf, ec)
  nc = (nf - 1) / 2
  nc.times do |ic|
    fi = 2 * ic + 1
    nc.times do |jc|
      fj = 2 * jc + 1; k = fi * nf + fj; v = ec[ic * nc + jc]
      uf[k] += v
      uf[k - 1] += 0.5 * v; uf[k + 1] += 0.5 * v; uf[k - nf] += 0.5 * v; uf[k + nf] += 0.5 * v
      uf[k - nf - 1] += 0.25 * v; uf[k - nf + 1] += 0.25 * v; uf[k + nf - 1] += 0.25 * v; uf[k + nf + 1] += 0.25 * v
    end
  end
end
def solve_dense3(b, x)
  m = Array.new(9) { Array.new(10, 0.0) }
  3.times { |i| 3.times { |j|
    r = i * 3 + j; m[r][r] = 4.0
    m[r][r - 1] = -1.0 if j > 0; m[r][r + 1] = -1.0 if j < 2
    m[r][r - 3] = -1.0 if i > 0; m[r][r + 3] = -1.0 if i < 2
    m[r][9] = b[r]
  } }
  9.times do |c|
    piv = (c...9).max_by { |r| m[r][c].abs }
    m[c], m[piv] = m[piv], m[c] if piv != c
    ((c + 1)...9).each { |r| f = m[r][c] / m[c][c]; (c...10).each { |k| m[r][k] -= f * m[c][k] } }
  end
  8.downto(0) { |r| s = m[r][9]; ((r + 1)...9).each { |k| s -= m[r][k] * x[k] }; x[r] = s / m[r][r] }
end
def vcycle(n, u, b)
  if n == 3; solve_dense3(b, u); return; end
  rbgs_sweep(n, u, b); rbgs_sweep(n, u, b)
  r = Array.new(n * n, 0.0); compute_res(n, u, b, r)
  rc = restrict_fw(n, r); nc = (n - 1) / 2; ec = Array.new(nc * nc, 0.0)
  vcycle(nc, ec, rc); prolong_add(n, u, ec); rbgs_sweep(n, u, b)
end
def fmg_solve(n, u, b)
  if n == 3; solve_dense3(b, u); return; end
  bc = restrict_fw(n, b); nc = (n - 1) / 2; uc = Array.new(nc * nc, 0.0)
  fmg_solve(nc, uc, bc); u.fill(0.0); prolong_add(n, u, uc); vcycle(n, u, b)
end

FileUtils.mkdir_p("out/ruby/figures")
$res = Array.new(7) { Array.new(4) }
$hist63 = Array.new(7) { [] }
$disc_err = Array.new(4, 0.0)
$stage_s = Array.new(19, 0.0)

def run_solver(meth, n, save63)
  n2 = n * n; b, ue = build_prob(n)
  u = Array.new(n2, 0.0); r = Array.new(n2, 0.0); p = Array.new(n2, 0.0); ap = Array.new(n2, 0.0); z = Array.new(n2, 0.0)
  nb = norm2_seq(b); om = omega_opt(n)
  max_it = [100000, 50000, 50000, 10000, 5000, 3000, 100][meth]
  hist = []; rr = 0.0; rz = 0.0
  if meth == 4 || meth == 5; r = b.dup; rr = dot_seq(r, r); end
  p = r.dup if meth == 4
  if meth == 5
    n.times { |i| n.times { |j| k = i * n + j; z[k] = 0.25 * ((r[k] + (j > 0 ? z[k - 1] : 0.0)) + (i > 0 ? z[k - n] : 0.0)) } }
    (n - 1).downto(0) { |i| (n - 1).downto(0) { |j| k = i * n + j; p[k] = z[k] + 0.25 * ((j + 1 < n ? p[k + 1] : 0.0) + (i + 1 < n ? p[k + n] : 0.0)) } }
    z = p.dup; rz = dot_seq(r, z)
  end
  max_it.times do
    rel = 0.0
    if meth == 0
      n.times { |i| row = i * n; n.times { |j|
        k = row + j; w = j > 0 ? u[k - 1] : 0.0; e = j + 1 < n ? u[k + 1] : 0.0
        s = i > 0 ? u[k - n] : 0.0; nn = i + 1 < n ? u[k + n] : 0.0
        z[k] = 0.25 * ((((b[k] + w) + e) + s) + nn)
      } }
      u, z = z, u; rel = compute_res(n, u, b, r) / nb
    elsif meth == 1 || meth == 3
      n.times { |i| row = i * n; n.times { |j|
        k = row + j; w = j > 0 ? u[k - 1] : 0.0; e = j + 1 < n ? u[k + 1] : 0.0
        s = i > 0 ? u[k - n] : 0.0; nn = i + 1 < n ? u[k + n] : 0.0
        gs = 0.25 * ((((b[k] + w) + e) + s) + nn)
        u[k] = meth == 1 ? gs : ((1.0 - om) * u[k] + om * gs)
      } }
      rel = compute_res(n, u, b, r) / nb
    elsif meth == 2
      rbgs_sweep(n, u, b); rel = compute_res(n, u, b, r) / nb
    elsif meth == 4
      apply_a(n, p, ap); alpha = rr / dot_seq(p, ap); k = 0
      while k < n2; u[k] += alpha * p[k]; r[k] -= alpha * ap[k]; k += 1; end
      rr_new = dot_seq(r, r); rel = Math.sqrt(rr_new) / nb; hist << rel; break if rel <= 1e-8
      beta = rr_new / rr; k = 0; while k < n2; p[k] = r[k] + beta * p[k]; k += 1; end
      rr = rr_new; next
    elsif meth == 5
      apply_a(n, p, ap); alpha = rz / dot_seq(p, ap); k = 0
      while k < n2; u[k] += alpha * p[k]; r[k] -= alpha * ap[k]; k += 1; end
      rel = norm2_seq(r) / nb; hist << rel; break if rel <= 1e-8
      n.times { |i| n.times { |j| k = i * n + j; ap[k] = 0.25 * ((r[k] + (j > 0 ? ap[k - 1] : 0.0)) + (i > 0 ? ap[k - n] : 0.0)) } }
      (n - 1).downto(0) { |i| (n - 1).downto(0) { |j| k = i * n + j; z[k] = ap[k] + 0.25 * ((j + 1 < n ? z[k + 1] : 0.0) + (i + 1 < n ? z[k + n] : 0.0)) } }
      rz_new = dot_seq(r, z); beta = rz_new / rz; k = 0; while k < n2; p[k] = z[k] + beta * p[k]; k += 1; end
      rz = rz_new; next
    else
      vcycle(n, u, b); rel = compute_res(n, u, b, r) / nb
    end
    hist << rel; break if rel <= 1e-8
  end
  used = hist.length
  rho = used >= 51 ? (hist[-1] / hist[-51])**0.02 : (used >= 2 ? (hist[-1] / hist[0])**(1.0 / (used - 1)) : 0.0)
  $hist63[meth] = hist if save63 && n == 63
  { iters: used, rel_res: hist[-1], err_inf: max_abs_diff(u, ue), rho_emp: rho, hist_hash: fnv1a_f64(hist), sol_hash: fnv1a_f64(u) }
end

def s01
  u = (1..9).map(&:to_f); au = Array.new(9, 0.0); apply_a(3, u, au)
  h1 = fnv1a_f64(u); h2 = fnv1a_f64(au); h3 = fnv1a_f64([1.0, -2.5, 3.25])
  ok = max_abs_diff(au, [-2.0,-1.0,4.0,3.0,0.0,7.0,16.0,11.0,22.0]) == 0.0 && h1 == 0xc56de10c6671e6c2 && h2 == 0x59ed72ec82ee959c && h3 == 0xafc5dc1ca48a2366
  File.write("out/ruby/s01_selftest.txt", format("ok,h_u3,h_Au3,h_test\n%d,%016x,%016x,%016x\n", ok ? 1 : 0, h1, h2, h3))
  exit(1) unless ok
end
def s02
  File.open("out/ruby/s02_problem.csv", "w") do |f|
    f.puts "N,norm_b,h_b,h_ue"
    NS_ALL.each { |n| b, ue = build_prob(n); f.puts format("%d,%.17g,%016x,%016x", n, norm2_seq(b), fnv1a_f64(b), fnv1a_f64(ue)) }
  end
end
def s03
  File.open("out/ruby/s03_eigen.csv", "w") do |f|
    f.puts "N,j,k,lambda_exact,rel_err_inf"
    NS_ALL[0, 3].each do |n|
      h = 1.0 / (n + 1); v = Array.new(n * n, 0.0); av = Array.new(n * n, 0.0)
      [[1,1],[1,2],[3,5]].each do |j, k|
        lam = 4.0 - 2.0 * cos_spec(j * PI_VAL * h) - 2.0 * cos_spec(k * PI_VAL * h)
        n.times { |r| n.times { |c| v[r * n + c] = sin_spec(j * PI_VAL * (r + 1) * h) * sin_spec(k * PI_VAL * (c + 1) * h) } }
        apply_a(n, v, av); me = (0...(n * n)).map { |q| (av[q] - lam * v[q]).abs }.max
        f.puts format("%d,%d,%d,%.17g,%.3e", n, j, k, lam, me / lam)
      end
    end
  end
end
def run_m(m, fn, nmax)
  File.open(fn, "w") do |f|
    f.puts "N,iters,rel_res,err_inf,rho_emp"
    nmax.times { |i| r = run_solver(m, NS_ALL[i], NS_ALL[i] == 63); $res[m][i] = r; f.puts format("%d,%d,%.10e,%.10e,%.10f", NS_ALL[i], r[:iters], r[:rel_res], r[:err_inf], r[:rho_emp]) }
  end
end
def s04; run_m(0, "out/ruby/s04_jacobi.csv", 3); end
def s05; run_m(1, "out/ruby/s05_gs.csv", 3); end
def s06; run_m(2, "out/ruby/s06_rbgs.csv", 3); end
def s07; run_m(3, "out/ruby/s07_sor.csv", 4); end
def s08; run_m(4, "out/ruby/s08_cg.csv", 4); end
def s09; run_m(5, "out/ruby/s09_pcg_ssor.csv", 4); end
def s10; run_m(6, "out/ruby/s10_mg.csv", 4); end
def s11
  File.open("out/ruby/s11_fmg.csv", "w") do |f|
    f.puts "N,fmg_err_inf,disc_err_inf,ratio"
    NS_ALL.each_with_index { |n, i| b, ue = build_prob(n); u = Array.new(n * n, 0.0); fmg_solve(n, u, b); fe = max_abs_diff(u, ue); $disc_err[i] = $res[6][i][:err_inf]; f.puts format("%d,%.10e,%.10e,%.6f", n, fe, $disc_err[i], fe / $disc_err[i]) }
  end
end
def s12
  File.open("out/ruby/s12_order.csv", "w") do |f|
    f.puts "N,h,err_inf,slope"
    xs = NS_ALL.map { |n| Math.log(1.0 / (n + 1)) }; ys = $disc_err.map { |e| Math.log(e) }
    mx = xs.sum / 4.0; my = ys.sum / 4.0; sl = (0...4).sum { |i| (xs[i] - mx) * (ys[i] - my) } / xs.sum { |x| (x - mx)**2 }
    4.times { |i| f.puts format("%d,%.10e,%.10e,%.6f", NS_ALL[i], 1.0 / (NS_ALL[i] + 1), $disc_err[i], sl) }
  end
end
def s13
  File.open("out/ruby/s13_scaling.csv", "w") do |f|
    f.puts "method,N,iters,slope,slope_se"
    7.times do |m|
      cnt = m <= 2 ? 3 : 4
      xs = (0...cnt).map { |i| Math.log(NS_ALL[i]) }; ys = (0...cnt).map { |i| Math.log($res[m][i][:iters]) }
      mx = xs.sum / cnt; my = ys.sum / cnt; den = xs.sum { |x| (x - mx)**2 }
      sl = (0...cnt).sum { |i| (xs[i] - mx) * (ys[i] - my) } / den; ic = my - sl * mx
      rss = (0...cnt).sum { |i| (ys[i] - (sl * xs[i] + ic))**2 }; se = Math.sqrt((rss / (cnt - 2)) / den)
      cnt.times { |i| f.puts format("%s,%d,%d,%.6f,%.6f", MNAMES[m], NS_ALL[i], $res[m][i][:iters], sl, se) }
    end
  end
end
def s14
  File.open("out/ruby/s14_history.csv", "w") do |f|
    f.puts "method,iter,rel_res"
    7.times { |m| ln = $hist63[m].length; step = ln > 200 ? ln / 200 : 1; (0...ln).step(step) { |k| f.puts format("%s,%d,%.10e", MNAMES[m], k + 1, $hist63[m][k]) }; f.puts format("%s,%d,%.10e", MNAMES[m], ln, $hist63[m][-1]) if (ln - 1) % step != 0 }
  end
end
def s15
  File.open("out/ruby/s15_mg_factor.csv", "w") do |f|
    f.puts "N,iters,mg_factor_last5"
    [63, 127].each do |n|
      b, = build_prob(n); u = Array.new(n * n, 0.0); r = Array.new(n * n, 0.0); nb = norm2_seq(b); h = []
      (1..20).each { vcycle(n, u, b); rel = compute_res(n, u, b, r) / nb; h << rel; break if rel <= 1e-8 }
      f.puts format("%d,%d,%.8f", n, h.length, (h[-1] / h[-6])**0.2)
    end
  end
end
def s16
  flops = [9, 9, 9, 12, 21, 35, 36]
  File.open("out/ruby/s16_cost.csv", "w") do |f|
    f.puts "method,N,flops_per_unk_iter,ns_med,ns_mad"
    n = 63; n2 = n * n; b, = build_prob(n); u = Array.new(n2, 0.0); z = Array.new(n2, 0.0)
    7.times do |m|
      iters = m == 6 ? 3 : 10
      vals = Array.new(7) do
        u.fill(0.0); t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
        iters.times { if m == 6; vcycle(n, u, b); elsif m == 2; rbgs_sweep(n, u, b); else; apply_a(n, u, z); k = 0; while k < n2; u[k] += 0.25 * (b[k] - z[k]); k += 1; end; end }
        (Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0) * 1e9 / (iters * n2)
      end.sort
      med = vals[3]; dev = vals.map { |v| (v - med).abs }.sort
      f.puts format("%s,%d,%d,%.4f,%.4f", MNAMES[m], n, flops[m], med, dev[3])
    end
  end
end
def s17
  File.open("out/ruby/s17_parity.csv", "w") do |f|
    f.puts "method,N,iters,hist_fnv1a,sol_fnv1a"
    7.times { |m| 2.times { |ni| r = $res[m][ni]; f.puts format("%s,%d,%d,%016x,%016x", MNAMES[m], NS_ALL[ni], r[:iters], r[:hist_hash], r[:sol_hash]) } }
  end
end
def s18
  cols = %w[#1f77b4 #ff7f0e #2ca02c #d62728 #9467bd #8c564b #e377c2]
  s1 = +"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>"
  7.times { |m| ln = $hist63[m].length; step = ln > 150 ? ln / 150 : 1; pts = (0...ln).step(step).map { |k| format("%.1f,%.1f", 60.0 + 500.0 * Math.log10(k + 1.0) / Math.log10($hist63[0].length), 310.0 - 280.0 * [[-Math.log10([$hist63[m][k], 1e-12].max) / 9.0, 0.0].max, 1.0].min) }.join(" "); s1 << "<polyline fill=\"none\" stroke=\"#{cols[m]}\" stroke-width=\"1.4\" points=\"#{pts}\"/>" }
  File.write("out/ruby/figures/history.svg", s1 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby Residual History (N=63)</text></svg>\n")
  s2 = +"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>"
  7.times { |m| cnt = m <= 2 ? 3 : 4; pts = (0...cnt).map { |i| format("%.1f,%.1f", 60.0 + 500.0 * (Math.log10(NS_ALL[i]) - Math.log10(15)) / (Math.log10(127) - Math.log10(15)), 310.0 - 280.0 * Math.log10($res[m][i][:iters]) / 5.0) }.join(" "); s2 << "<polyline fill=\"none\" stroke=\"#{cols[m]}\" stroke-width=\"1.5\" points=\"#{pts}\"/>" }
  File.write("out/ruby/figures/scaling.svg", s2 + "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby Iterations vs N</text></svg>\n")
  pts3 = (0...4).map { |i| format("%.1f,%.1f", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * (Math.log10($disc_err[i]) + 5.0) / 4.0) }.join(" ")
  File.write("out/ruby/figures/order.svg", "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"1.8\" points=\"#{pts3}\"/><text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">Ruby Discretization Error vs h</text></svg>\n")
end
def s19
  File.write("out/ruby/s19_summary.json", format("{\"lang\":\"ruby\",\"N_max\":127,\"N_stat_max\":63,\"toolchain\":\"ruby 4.0.9\",\"total_sec\":%.3f}\n", $stage_s.sum))
end

cmd = ARGV[0] || "all"
stages = (1..19).map { |i| format("s%02d", i) }
if cmd == "all"
  stages.each_with_index { |nm, i| puts "[ruby] tahap #{i + 1}/19 (#{nm})"; t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC); send(nm); $stage_s[i] = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0 }
elsif stages.include?(cmd); send(cmd); else exit(1); end
