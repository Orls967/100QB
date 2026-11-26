# frozen_string_literal: true

require "fileutils"
require_relative "rng"
require_relative "config"
require_relative "dtmc"
require_relative "station"
require_relative "passenger"
require_relative "network"
require_relative "engine"
require_relative "stats"
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
      out.puts(format("%s%s,%d,%d,%016x", prefix, nm, b.count("\n"), b.bytesize, fnv1a_bytes(b)))
    end
  end
end

out_dir = "out/ruby"
FileUtils.mkdir_p("#{out_dir}/figures")
cfg = load_config("config.kv")
if ARGV[0] == "s01"
  write_s01_to_s13(cfg, out_dir)
  exit 0
end

tot = 0.0
(1..13).each { |st| puts format("[ruby] tahap %d/48 (s%02d)", st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s01_to_s13(cfg, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

main_days = 7
puts "[ruby] tahap 14/48 (s14)"
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
rpm = default_rp(cfg, main_days); rpm.audit_on = 1
s_main = Sim.new(cfg, rpm); sim_run(s_main)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

(15..41).each { |st| puts format("[ruby] tahap %d/48 (s%02d)", st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s14_to_s45(cfg, out_dir, s_main)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts "[ruby] tahap 42/48 (s42)"
vals = 5.times.map do
  sb = Sim.new(cfg, default_rp(cfg, 3))
  tb = Process.clock_gettime(Process::CLOCK_MONOTONIC)
  sim_run(sb)
  dt = [Process.clock_gettime(Process::CLOCK_MONOTONIC) - tb, 1e-9].max
  sb.sys_sum_l.to_f / dt
end.sort
med = vals[2]
dev = vals.map { |v| (v - med).abs }.sort
File.write("#{out_dir}/s42_bench.csv", format("slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,%.2f,%.2f\n", med, dev[2]))

puts "[ruby] tahap 43/48 (s43)\n[ruby] tahap 44/48 (s44)\n[ruby] tahap 45/48 (s45)\n[ruby] tahap 46/48 (s46)"
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s46_svgs(s_main, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts "[ruby] tahap 47/48 (s47)\n[ruby] tahap 48/48 (s48)"
File.write("#{out_dir}/s48_summary.json", format("{\"lang\":\"ruby\",\"days\":%d,\"sw_days\":1,\"toolchain\":\"ruby 4.0.9\",\"total_sec\":%.3f}\n", main_days, tot))
File.open("#{out_dir}/files_manifest.csv", "w") { |mf| mf.puts("file,lines,bytes,fnv1a"); manifest_scan(out_dir, "", mf) }
