# frozen_string_literal: true

require 'fileutils'

module BacktestOverfitting
  module Svg
    WIDTH = 860
    HEIGHT = 500
    PAD_L = 82
    PAD_R = 30
    PAD_T = 48
    PAD_B = 64

    module_function

    def write_line_chart(path, title:, x_label:, y_label:, series:)
      FileUtils.mkdir_p(File.dirname(path))
      xmin, xmax = bounds(series.flat_map { |s| s[:x] })
      ymin, ymax = bounds(series.flat_map { |s| s[:y] })
      lines = [header(title, x_label, y_label, xmin, xmax, ymin, ymax)]
      series.each_with_index { |s, idx| append_series(lines, s, idx, [xmin, xmax, ymin, ymax]) }
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def write_bar_chart(path, title:, y_label:, groups:)
      FileUtils.mkdir_p(File.dirname(path))
      ymax = [groups.flat_map { |g| g[:bars].map { |b| b[:value] } }.max * 1.12, 0.1].max
      lines = [header(title, '', y_label, 0.0, groups.length.to_f, 0.0, ymax)]
      draw_bar_groups(lines, groups, ymax)
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def write_histogram(path, title:, x_label:, y_label:, series_vals:, bins: 22)
      FileUtils.mkdir_p(File.dirname(path))
      all_v = series_vals.flat_map { |s| s[:values] }
      lo, hi = bounds(all_v)
      hists = series_vals.map { |s| compute_hist(s[:values], lo, hi, bins) }
      ymax = [hists.flatten.max * 1.15, 0.05].max
      bw = (hi - lo) / bins
      xs = Array.new(bins) { |i| lo + ((i + 0.5) * bw) }
      lines = [header(title, x_label, y_label, lo, hi, 0.0, ymax)]
      series_vals.each_with_index do |s, idx|
        append_series(lines, s.merge(x: xs, y: hists[idx]), idx, [lo, hi, 0.0, ymax])
      end
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def compute_hist(vals, lo, hi, bins)
      counts = Array.new(bins, 0.0)
      span = [hi - lo, 1.0e-12].max
      vals.each do |v|
        bi = (((v - lo) / span) * bins).floor.clamp(0, bins - 1)
        counts[bi] += 1.0
      end
      tot = [vals.length, 1].max
      counts.map { |c| c / tot }
    end
    private_class_method :compute_hist

    def append_series(lines, series, idx, lims)
      xmin, xmax, ymin, ymax = lims
      pts = series[:x].each_with_index.map do |xv, i|
        "#{sx(xv, xmin, xmax).round(2)},#{sy(series[:y][i], ymin, ymax).round(2)}"
      end.join(' ')
      dash = series[:dashed] ? ' stroke-dasharray="6,4"' : ''
      lines << %(<polyline fill="none" stroke="#{series[:color]}" stroke-width="2.3"#{dash} points="#{pts}"/>)
      lx = PAD_L + 14 + (idx * 210)
      lines << format('<line x1="%<x1>d" y1="30" x2="%<x2>d" y2="30" stroke="%<c>s" stroke-width="2.5"%<d>s/>',
                      x1: lx, x2: lx + 22, c: series[:color], d: dash)
      lines << %(<text x="#{lx + 28}" y="34" font-size="11" fill="#1f2937">#{series[:label]}</text>)
    end
    private_class_method :append_series

    def draw_bar_groups(lines, groups, ymax)
      pw = WIDTH - PAD_L - PAD_R
      gw = pw.to_f / groups.length
      groups.each_with_index do |g, gi|
        gx = PAD_L + (gi * gw) + (0.12 * gw)
        bw = (0.76 * gw) / g[:bars].length
        g[:bars].each_with_index do |b, bi|
          y0 = sy(0.0, 0.0, ymax)
          yt = sy(b[:value], 0.0, ymax)
          lines << format('<rect x="%<x>.1f" y="%<y>.1f" width="%<w>.1f" height="%<h>.1f" fill="%<c>s"/>',
                          x: gx + (bi * bw), y: yt, w: [bw - 3, 4.0].max, h: y0 - yt, c: b[:color])
        end
        tx = (PAD_L + ((gi + 0.5) * gw)).round(1)
        lines << (%(<text x="#{tx}" y="#{HEIGHT - 42}" text-anchor="middle" font-size="10" ) +
                  %(fill="#374151">#{g[:label]}</text>))
      end
    end
    private_class_method :draw_bar_groups

    def bounds(vals)
      lo, hi = vals.compact.minmax
      return [lo - 1.0, hi + 1.0] if (hi - lo).abs < 1.0e-12

      pad = 0.06 * (hi - lo)
      [lo - pad, hi + pad]
    end
    private_class_method :bounds

    def sx(val, xmin, xmax) = PAD_L + ((val - xmin) / (xmax - xmin) * (WIDTH - PAD_L - PAD_R))
    private_class_method :sx

    def sy(val, ymin, ymax) = HEIGHT - PAD_B - ((val - ymin) / (ymax - ymin) * (HEIGHT - PAD_T - PAD_B))
    private_class_method :sy

    def header(title, x_label, y_label, _xmin, _xmax, ymin, ymax)
      pw = WIDTH - PAD_L - PAD_R
      ph = HEIGHT - PAD_T - PAD_B
      ymid = PAD_T + (ph / 2)
      [
        %(<svg xmlns="http://www.w3.org/2000/svg" width="#{WIDTH}" height="#{HEIGHT}" ) +
          %(viewBox="0 0 #{WIDTH} #{HEIGHT}">),
        %(<rect width="#{WIDTH}" height="#{HEIGHT}" fill="#ffffff"/>),
        %(<text x="#{WIDTH / 2}" y="20" text-anchor="middle" font-size="14" ) +
          %(font-weight="bold" fill="#111827">#{title}</text>),
        %(<rect x="#{PAD_L}" y="#{PAD_T}" width="#{pw}" height="#{ph}" fill="#f9fafb" stroke="#d1d5db"/>),
        %(<text x="#{PAD_L + (pw / 2)}" y="#{HEIGHT - 16}" text-anchor="middle" ) +
          %(font-size="12" fill="#374151">#{x_label}</text>),
        %(<text x="18" y="#{ymid}" text-anchor="middle" font-size="12" fill="#374151" ) +
          %(transform="rotate(-90 18 #{ymid})">#{y_label}</text>),
        %(<text x="#{PAD_L - 6}" y="#{HEIGHT - PAD_B}" text-anchor="end" ) +
          %(font-size="10" fill="#6b7280">#{format('%.3g', ymin)}</text>),
        %(<text x="#{PAD_L - 6}" y="#{PAD_T + 8}" text-anchor="end" ) +
          %(font-size="10" fill="#6b7280">#{format('%.3g', ymax)}</text>)
      ].join("\n")
    end
    private_class_method :header
  end
end
