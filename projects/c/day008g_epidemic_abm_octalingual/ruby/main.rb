# frozen_string_literal: true

require "fileutils"
require_relative "rng"
require_relative "population"
require_relative "contact"
require_relative "disease"
require_relative "interventions"
require_relative "engine"
require_relative "analysis"
require_relative "svg"

def manifest_scan(dir, prefix, out)
  Dir.children(dir).sort.each do |nm|
    next if nm.start_with?(".")

    p = File.join(dir, nm)
    if File.directory?(p)
      manifest_scan(p, "#{prefix}#{nm}/", out)
    elsif nm != "files_manifest.csv"
      b = File.binread(p)
      lc = b.count("\n")
      out.puts(format("%s%s,%d,%d,%016x", prefix, nm, lc, b.bytesize, fnv1a_bytes(b)))
    end
  end
end

out_dir = "out/ruby"
FileUtils.mkdir_p("#{out_dir}/figures")
cfg = load_config("config.kv")
if ARGV[0] == "s01"
  write_s01_to_s09(cfg, out_dir)
  exit 0
end

tot = 0.0
(1..9).each { |st| puts format("[ruby] tahap %d/39 (s%02d)", st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s01_to_s09(cfg, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

n_main = 3000
beta_cal = 400
puts "[ruby] tahap 10/39 (s10)"
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
s_no = Sim.new(cfg, default_rp(n_main, cfg["main_days"], beta_cal, 0), true)
sim_run(s_no)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts "[ruby] tahap 11/39 (s11)"
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
s_main = Sim.new(cfg, default_rp(n_main, cfg["main_days"], beta_cal, 1), true)
sim_run(s_main)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

(12..34).each { |st| puts format("[ruby] tahap %d/39 (s%02d)", st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s10_to_s37(s_no, s_main, cfg, out_dir, beta_cal)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts "[ruby] tahap 35/39 (s35)"
vals = 5.times.map do
  sb = Sim.new(cfg, default_rp(2000, 60, beta_cal, 1), false)
  tb = Process.clock_gettime(Process::CLOCK_MONOTONIC)
  sim_run(sb)
  dt = [Process.clock_gettime(Process::CLOCK_MONOTONIC) - tb, 1e-9].max
  120_000.0 / dt
end.sort
med = vals[2]
dev = vals.map { |v| (v - med).abs }.sort
File.write("#{out_dir}/s35_bench.csv", format("person_days,pd_per_sec_med,pd_per_sec_mad\n120000,%.2f,%.2f\n", med, dev[2]))

puts "[ruby] tahap 36/39 (s36)\n[ruby] tahap 37/39 (s37)\n[ruby] tahap 38/39 (s38)"
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s38_svgs(s_no, s_main, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts "[ruby] tahap 39/39 (s39)"
File.write("#{out_dir}/s39_summary.json", format("{\"lang\":\"ruby\",\"N\":%d,\"days\":%d,\"N_sw\":500,\"toolchain\":\"ruby 4.0.9\",\"total_sec\":%.3f}\n", n_main, cfg["main_days"], tot))
File.open("#{out_dir}/files_manifest.csv", "w") do |mf|
  mf.puts("file,lines,bytes,fnv1a")
  manifest_scan(out_dir, "", mf)
end
