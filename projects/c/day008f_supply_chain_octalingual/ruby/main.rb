require 'fileutils'
%w[rng inventory allocation entities engine clearing analysis svg].each { |m| require_relative m }

out_dir = 'out/ruby'
FileUtils.mkdir_p("#{out_dir}/figures")
cfg = load_config('config.kv')
cmd = ARGV[0] || 'all'
if cmd == 's01'
  write_s01_to_s09(cfg, out_dir)
  exit 0
end

tot = 0.0
(1..9).each { |st| puts format('[ruby] tahap %d/40 (s%02d)', st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s01_to_s09(cfg, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts '[ruby] tahap 10/40 (s10)'
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
sim = Sim.new(cfg, true, out_dir)
sim_run_days(sim, cfg['num_days'], 8000, 0, -1, 1)
sim.close_logs
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

(11..38).each { |st| puts format('[ruby] tahap %d/40 (s%02d)', st, st) }
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s11_to_s38(sim, cfg, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts '[ruby] tahap 39/40 (s39)'
t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
write_s39_svgs(sim, out_dir)
tot += Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

puts '[ruby] tahap 40/40 (s40)'
File.write("#{out_dir}/s40_summary.json", format("{\"lang\":\"ruby\",\"days\":5,\"ev_per_day\":8000,\"ev_total\":%d,\"toolchain\":\"ruby 4.0.9\",\"total_sec\":%.3f}\n", cfg['num_days'] * 8000, tot))
mf = ['file,lines,bytes,fnv1a']
Dir.glob("#{out_dir}/**/*").sort.each do |p|
  next if File.directory?(p) || File.basename(p) == 'files_manifest.csv'
  rel = p.sub("#{out_dir}/", '')
  bytes = File.binread(p)
  mf << format('%s,%d,%d,%016x', rel, bytes.count("\n"), bytes.bytesize, fnv1a_bytes(bytes))
end
File.write("#{out_dir}/files_manifest.csv", mf.join("\n") + "\n")
