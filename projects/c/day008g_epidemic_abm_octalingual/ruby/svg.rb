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

def write_s38_svgs(s_no, s_main, out_dir)
  d_cnt = s_main.rp.days
  step = [d_cnt / 12, 1].max
  epi = Array.new(12, 0.0)
  inc = Array.new(12, 0.0)
  rt = Array.new(12, 0.0)
  hosp = Array.new(12, 0.0)
  ens = Array.new(12, 0.0)
  12.times do |i|
    d = [i * step, d_cnt - 1].min
    epi[i] = (s_main.daily_cnt[d][ST_IA] + s_main.daily_cnt[d][ST_IS]).to_f
    inc[i] = s_main.daily_inc_age[d].sum.to_f
    rt[i] = (s_no.daily_cnt[d][ST_IA] + s_no.daily_cnt[d][ST_IS]).to_f
    hosp[i] = s_main.daily_hosp[d][0].to_f
    ens[i] = s_main.daily_var[d][1].to_f
  end
  write_bar_svg(out_dir, "epi_curves.svg", "Ruby Active Infectious (Interv)", epi, "#1f77b4")
  write_bar_svg(out_dir, "incidence_age.svg", "Ruby Daily Incidence Sample", inc, "#ff7f0e")
  write_bar_svg(out_dir, "rt_series.svg", "Ruby Active Infectious (No Interv)", rt, "#d62728")
  write_bar_svg(out_dir, "hospital_cap.svg", "Ruby Hospital Occupancy Sample", hosp, "#9467bd")
  write_bar_svg(out_dir, "layer_attack.svg", "Ruby Infections by Layer (HH/SW/Comm)", s_main.layer_inf.map(&:to_f), "#2ca02c")
  write_bar_svg(out_dir, "offspring_dist.svg", "Ruby Infections by Age Group", s_main.age_inf.map(&:to_f), "#8c564b")
  write_bar_svg(out_dir, "sweep_heatmap.svg", "Ruby No-Interv vs Interv Total Infections", [s_no.total_trans.to_f, s_main.total_trans.to_f], "#e377c2")
  write_bar_svg(out_dir, "vaccine_effect.svg", "Ruby Vaccine Doses &amp; Blocked Transmissions", [s_main.branches[18].to_f, s_main.branches[3].to_f], "#17becf")
  write_bar_svg(out_dir, "variant_share.svg", "Ruby Variant 1 Daily Incidence", ens, "#bcbd22")
  write_bar_svg(out_dir, "ensemble_band.svg", "Ruby Final Recovered vs Deaths", [s_main.daily_cnt[d_cnt - 1][ST_R].to_f, s_main.daily_cnt[d_cnt - 1][ST_D].to_f], "#7f7f7f")
end
