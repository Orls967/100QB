# frozen_string_literal: true

def write_bar_svg(out_dir, fname, title, vals, col)
  mx = ([1.0] + vals.map(&:to_f)).max
  bw = 460.0 / [vals.length, 1].max
  File.open("#{out_dir}/figures/#{fname}", "w") do |f|
    f.write('<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="35" width="500" height="275" fill="none" stroke="#ccc"/>')
    vals.each_with_index do |v, i|
      h = [(v.to_f / mx) * 250.0, 1.0].max
      f.write(format('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>', 75.0 + i * bw, 310.0 - h, bw * 0.72, h, col))
    end
    f.puts(format('<text x="300" y="22" text-anchor="middle" font-family="sans-serif" font-size="13">%s</text></svg>', title))
  end
end

def write_s46_svgs(s_main, out_dir)
  qstn = (0...6).map { |i| s_main.stn[i].arr_cnt.to_f }
  utilh = (0...12).map { |i| s_main.hr_util_sum[i * 2][STN_SECURITY].to_f }
  wcls = (0...3).map { |i| s_main.cls_wait_cnt[i] > 0 ? s_main.cls_wait_sum[i].to_f / s_main.cls_wait_cnt[i] : 0.0 }
  soj = (0...3).map { |i| s_main.cls_sojourn_cnt[i] > 0 ? s_main.cls_sojourn_sum[i].to_f / s_main.cls_sojourn_cnt[i] : 0.0 }
  br = (0...6).map { |i| s_main.branches[i].to_f }
  write_bar_svg(out_dir, "queue_stations.svg", "Ruby Arrivals by Station", qstn, "#1f77b4")
  write_bar_svg(out_dir, "util_hourly.svg", "Ruby Security Busy Slots (2h sample)", utilh, "#ff7f0e")
  write_bar_svg(out_dir, "wait_hist.svg", "Ruby Mean Wait by Station", qstn, "#2ca02c")
  write_bar_svg(out_dir, "wait_quantiles.svg", "Ruby Mean Wait by Class (Econ/Biz/Trans)", wcls, "#d62728")
  write_bar_svg(out_dir, "dtmc_validation.svg", "Ruby Departures by Station", br, "#9467bd")
  write_bar_svg(out_dir, "server_heatmap.svg", "Ruby Arrivals vs Done vs Lost", [s_main.total_arr.to_f, s_main.total_done.to_f, s_main.total_lost.to_f], "#8c564b")
  write_bar_svg(out_dir, "priority_effect.svg", "Ruby Priority Jumps &amp; Extra Security", [s_main.branches[BR_BIZ_PRIO].to_f, s_main.branches[BR_SEC_EXTRA].to_f], "#e377c2")
  write_bar_svg(out_dir, "breakdown_effect.svg", "Ruby Breakdowns &amp; Repairs", [s_main.branches[BR_SRV_BREAK].to_f, s_main.branches[BR_SRV_REPAIR].to_f], "#7f7f7f")
  write_bar_svg(out_dir, "bottleneck.svg", "Ruby Station Busy Slots", [s_main.stn[0].busy_slots.to_f, s_main.stn[1].busy_slots.to_f, s_main.stn[2].busy_slots.to_f, s_main.stn[3].busy_slots.to_f], "#bcbd22")
  write_bar_svg(out_dir, "ensemble_band.svg", "Ruby Abandonments vs Missed Flights", [s_main.total_aban.to_f, s_main.total_missed.to_f], "#17becf")
  write_bar_svg(out_dir, "crn_ratio.svg", "Ruby Baggage Total vs Missed vs Blocked", [s_main.total_bags.to_f, s_main.total_bags_miss.to_f, s_main.total_belt_block.to_f], "#1f77b4")
  write_bar_svg(out_dir, "tornado.svg", "Ruby Class Sojourn Means", soj, "#ff7f0e")
end
