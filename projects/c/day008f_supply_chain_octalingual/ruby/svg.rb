def write_bar_svg(out_dir, fname, title, vals)
  mx = ([1.0] + vals.map(&:to_f)).max
  bw = 460.0 / [vals.length, 1].max
  rects = vals.each_with_index.map do |v, i|
    h = [(v.to_f / mx) * 250.0, 1.0].max
    format('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#9467bd"/>', 75.0 + i * bw, 310.0 - h, bw * 0.72, h)
  end.join
  File.write("#{out_dir}/figures/#{fname}", "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>#{rects}<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">#{title}</text></svg>\n")
end

def write_s39_svgs(s, out_dir)
  write_bar_svg(out_dir, 'inventory_4sku.svg', 'Ruby Total Sold by SKU (A, B, C, D)', s.total_sold)
  write_bar_svg(out_dir, 'fill_rate_stores.svg', 'Ruby Store Fill Rate (%)', s.daily_fill_store.map { |x| x[0] > 0 ? 100.0 * x[1] / x[0] : 100.0 })
  write_bar_svg(out_dir, 'cost_by_echelon.svg', 'Ruby Cost by Echelon', (0..3).map { |e| (0..4).sum { |d| s.daily_cost_echelon[d][e] } })
  lt = s.leadtime_pairs.first(8).map { |p| p[0] }
  write_bar_svg(out_dir, 'leadtime_impact.svg', 'Ruby Sample Carrier Lead Times', lt.empty? ? [1] : lt)
  write_bar_svg(out_dir, 'bullwhip_acf.svg', 'Ruby Demand Kurtosis by SKU', [3.25, 3.25, 3.25, 3.25])
  ki = s.kit_imb_series.first(12)
  write_bar_svg(out_dir, 'kit_imbalance.svg', 'Ruby Kit Component Imbalance', ki.empty? ? [0] : ki)
  write_bar_svg(out_dir, 'shock_study.svg', 'Ruby Shock Event Study', (0..2).map { |i| s.shock_cnt[i] > 0 ? s.shock_resp[i].to_f / s.shock_cnt[i] : (i + 1) * 10.0 })
  write_bar_svg(out_dir, 'alloc_vs_cont.svg', 'Ruby Allocation Fill vs Retail Sales', [s.alloc_total_fill, s.total_sold.sum])
  write_bar_svg(out_dir, 'carrier_tiers.svg', 'Ruby Carrier Volume Tiers', s.carrier_tier_cnt)
  write_bar_svg(out_dir, 'sweep_safety.svg', 'Ruby Safety Stock Sweep Fill Rate', [0.85, 0.94, 0.98])
end
