# frozen_string_literal: true

module Longevity
  module Svg
    module_function

    def fmt(val, prec = 2)
      format("%<val>.#{prec}f", val: val)
    end

    def svg_open(width, height, title)
      <<~SVG
        <svg xmlns="http://www.w3.org/2000/svg" width="#{width}" height="#{height}" viewBox="0 0 #{width} #{height}">
        <rect width="100%" height="100%" fill="#ffffff"/>
        <style>text{font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;fill:#1f2937;}.title{font-size:15px;font-weight:700;}.axis{stroke:#4b5563;stroke-width:1.2;}.lbl{font-size:11px;fill:#4b5563;}</style>
        <text x="#{width / 2}" y="26" text-anchor="middle" class="title">#{title}</text>
      SVG
    end

    def axes_box(px0, px1, py0, py1)
      "<line x1=\"#{px0}\" y1=\"#{py1}\" x2=\"#{px1}\" y2=\"#{py1}\" class=\"axis\"/>\n" \
        "<line x1=\"#{px0}\" y1=\"#{py0}\" x2=\"#{px0}\" y2=\"#{py1}\" class=\"axis\"/>\n"
    end

    def line_el(x1_pt, y1_pt, x2_pt, y2_pt, stroke:, width: 1.8, dash: nil)
      d = dash ? " stroke-dasharray=\"#{dash}\"" : ''
      "<line x1=\"#{fmt(x1_pt)}\" y1=\"#{fmt(y1_pt)}\" x2=\"#{fmt(x2_pt)}\" y2=\"#{fmt(y2_pt)}\" " \
        "stroke=\"#{stroke}\"#{d} stroke-width=\"#{width}\"/>\n"
    end

    def poly_el(points, stroke:, width: 2.2)
      pts = points.map { |x, y| "#{fmt(x)},#{fmt(y)}" }.join(' ')
      "<polyline fill=\"none\" stroke=\"#{stroke}\" stroke-width=\"#{width}\" points=\"#{pts}\"/>\n"
    end

    def circle_el(center_x, center_y, fill:, rad: 4.5)
      "<circle cx=\"#{fmt(center_x)}\" cy=\"#{fmt(center_y)}\" r=\"#{rad}\" fill=\"#{fill}\"/>\n"
    end

    def footer_el(text_x, label)
      "<text x=\"#{text_x}\" y=\"398\" text-anchor=\"middle\" class=\"lbl\">#{label}</text>\n</svg>\n"
    end

    def write_mortality_surface(path, ages, years, ln_m)
      out = svg_open(760, 420, 'Lee-Carter Synthetic Log-Mortality Surface ln m(x,t)')
      out << axes_box(75.0, 710.0, 55.0, 365.0)
      mx = ->(idx) { 75.0 + ((idx.to_f / (ages.size - 1)) * 635.0) }
      my = ->(val) { 365.0 - (((val + 8.5) / 9.0) * 310.0) }

      [[0, '#dc2626'], [13, '#d97706'], [26, '#059669'], [39, '#2563eb']].each do |t_idx, col|
        pts = ages.each_index.map { |i| [mx.call(i), my.call(ln_m[i][t_idx])] }
        out << poly_el(pts, stroke: col, width: 2.0)
      end

      out << footer_el(392, "Age x (30..100) | Red: #{years[0]} -&gt; Blue: #{years[-1]}")
      File.write(path, out)
    end

    def write_life_expectancy_fan(path, fan_data)
      out = svg_open(760, 420, 'Cohort 2020 Survival Probability Fan Chart (Age 65, 5%/50%/95%)')
      out << axes_box(75.0, 710.0, 55.0, 365.0)
      surv_fan = fan_data[:survival_fan]
      mx = ->(idx) { 75.0 + ((idx.to_f / (surv_fan.size - 1)) * 635.0) }
      my = ->(val) { 365.0 - (val * 310.0) }

      upper = surv_fan.each_with_index.map { |r, i| "#{fmt(mx.call(i))},#{fmt(my.call(r[:q95]))}" }
      lower = surv_fan.each_with_index.reverse_each.map { |r, i| "#{fmt(mx.call(i))},#{fmt(my.call(r[:q05]))}" }
      out << "<polygon fill=\"#93c5fd\" fill-opacity=\"0.45\" points=\"#{(upper + lower).join(' ')}\"/>\n"

      med = surv_fan.each_with_index.map { |r, i| [mx.call(i), my.call(r[:q50])] }
      out << poly_el(med, stroke: '#1d4ed8', width: 2.4)
      lbl = "Age 66..100 | Median e65 = #{fmt(fan_data[:e65_q50])} yrs " \
            "(#{fmt(fan_data[:e65_q05])}..#{fmt(fan_data[:e65_q95])})"
      out << footer_el(392, lbl)
      File.write(path, out)
    end

    def write_k_moments_qq(path, k_paths, rw_params)
      out = svg_open(760, 420, 'Forecast k(T+40) Empirical vs Theoretical Normal Quantiles')
      out << axes_box(80.0, 710.0, 55.0, 365.0)
      vals = k_paths.map { |p| p[39] }.sort
      mu = rw_params[:k_last] + (40.0 * rw_params[:d_hat])
      sd = Math.sqrt(40.0) * rw_params[:sigma_hat]
      probs = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
      zs = [-2.3263, -1.6449, -1.2816, -0.6745, 0.0, 0.6745, 1.2816, 1.6449, 2.3263]
      theo = zs.map { |z| mu + (sd * z) }
      emp = probs.map { |p| Stats.quantile_sorted(vals, p) }
      vmin = [theo.min, emp.min].min - 1.0
      vmax = [theo.max, emp.max].max + 1.0
      mx = ->(v) { 80.0 + (((v - vmin) / (vmax - vmin)) * 630.0) }
      my = ->(v) { 365.0 - (((v - vmin) / (vmax - vmin)) * 310.0) }

      out << line_el(mx.call(vmin), my.call(vmin), mx.call(vmax), my.call(vmax), stroke: '#dc2626', dash: '5,4')
      theo.each_index { |i| out << circle_el(mx.call(theo[i]), my.call(emp[i]), fill: '#2563eb') }
      out << footer_el(395, 'Theoretical Normal Quantile (h=40) vs Empirical')
      File.write(path, out)
    end

    def write_cv_vs_n(path, cv_rows, sim_rows)
      out = svg_open(760, 420, 'Portfolio Coefficient of Variation CV(N) — Undiversifiable Longevity Floor')
      out << axes_box(80.0, 710.0, 55.0, 365.0)
      finite = cv_rows.reject { |r| r[:n] == 'inf' }
      inf_cv = cv_rows.find { |r| r[:n] == 'inf' }[:cv]
      ymax = finite.map { |r| r[:cv] }.max * 1.15
      mx = ->(i) { 80.0 + ((i.to_f / (finite.size - 1)) * 630.0) }
      my = ->(v) { 365.0 - ((v / ymax) * 310.0) }

      out << line_el(80.0, my.call(inf_cv), 710.0, my.call(inf_cv), stroke: '#059669', dash: '5,4')
      out << poly_el(finite.each_with_index.map { |r, i| [mx.call(i), my.call(r[:cv])] }, stroke: '#2563eb')
      finite.each_with_index do |r, i|
        out << circle_el(mx.call(i), my.call(r[:cv]), fill: '#2563eb')
        out << "<text x=\"#{fmt(mx.call(i))}\" y=\"381\" text-anchor=\"middle\" class=\"lbl\">N=#{r[:n]}</text>\n"
      end
      draw_sim_squares(out, finite, sim_rows, mx, my)
      out << footer_el(395, "Blue: Exact LTV | Red: Full Sim | Green: N=inf (#{fmt(inf_cv * 100.0)}%)")
      File.write(path, out)
    end

    def draw_sim_squares(out, finite, sim_rows, map_x, map_y)
      sim_rows.each do |sr|
        idx = finite.index { |r| r[:n] == sr[:n] }
        next unless idx

        x = fmt(map_x.call(idx) - 4.0)
        y = fmt(map_y.call(sr[:cv]) - 4.0)
        out << "<rect x=\"#{x}\" y=\"#{y}\" width=\"8\" height=\"8\" fill=\"#dc2626\"/>\n"
      end
    end

    def write_liability_hist_var(path, samples, mean_val, var_val)
      out = svg_open(760, 420, 'Portfolio Liability Distribution (N=100, S=5,000) with Mean &amp; VaR 99.5%')
      out << axes_box(75.0, 710.0, 55.0, 365.0)
      smin, smax = samples.minmax
      draw_hist_bars(out, samples, smin, smax)
      mx = ->(v) { 75.0 + (((v - smin) / (smax - smin)) * 635.0) }
      out << line_el(mx.call(mean_val), 55.0, mx.call(mean_val), 365.0, stroke: '#111827', width: 2.0)
      out << line_el(mx.call(var_val), 55.0, mx.call(var_val), 365.0, stroke: '#dc2626', dash: '5,4', width: 2.2)
      out << footer_el(392, "Black: Mean (#{fmt(mean_val, 1)}) | Red: VaR 99.5% (#{fmt(var_val, 1)})")
      File.write(path, out)
    end

    def draw_hist_bars(out, samples, smin, smax)
      bins = 32
      bw = (smax - smin) / bins.to_f
      counts = Array.new(bins, 0)
      samples.each { |v| counts[((v - smin) / bw).floor.clamp(0, bins - 1)] += 1 }
      cmax = counts.max.to_f
      bar_w = 635.0 / bins.to_f
      counts.each_with_index do |c, i|
        bh = (c / cmax) * 310.0
        x = fmt(75.0 + (i * bar_w))
        out << "<rect x=\"#{x}\" y=\"#{fmt(365.0 - bh)}\" width=\"#{fmt(bar_w - 1.0)}\" " \
               "height=\"#{fmt(bh)}\" fill=\"#60a5fa\"/>\n"
      end
    end

    def write_capital_vs_n(path, base_caps, param_caps)
      out = svg_open(760, 420, 'Solvency Capital (VaR 99.5% - Mean) as % of Mean vs N')
      out << axes_box(80.0, 710.0, 55.0, 365.0)
      ymax = [base_caps.map { |r| r[:capital_pct] }.max, param_caps.map { |r| r[:capital_pct] }.max].max * 1.15
      mx = ->(i) { 80.0 + ((i.to_f / (base_caps.size - 1)) * 630.0) }
      my = ->(v) { 365.0 - ((v / ymax) * 310.0) }

      out << poly_el(base_caps.each_with_index.map { |r, i| [mx.call(i), my.call(r[:capital_pct])] }, stroke: '#2563eb')
      out << poly_el(param_caps.each_with_index.map { |r, i|
        [mx.call(i), my.call(r[:capital_pct])]
      }, stroke: '#dc2626')
      base_caps.each_with_index do |r, i|
        out << circle_el(mx.call(i), my.call(r[:capital_pct]), fill: '#2563eb')
        out << circle_el(mx.call(i), my.call(param_caps[i][:capital_pct]), fill: '#dc2626')
        out << "<text x=\"#{fmt(mx.call(i))}\" y=\"381\" text-anchor=\"middle\" class=\"lbl\">N=#{r[:n]}</text>\n"
      end
      out << footer_el(395, 'Blue: Process Risk Only | Red: With Parameter Uncertainty')
      File.write(path, out)
    end
  end
end
