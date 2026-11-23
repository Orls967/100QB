# frozen_string_literal: true

require 'fileutils'

module FactorPortfolioLab
  # Pure-Ruby SVG chart renderer for line/multi-series, scatter, and bar charts.
  module Svg
    PALETTE = %w[#1f77b4 #d62728 #2ca02c #ff7f0e #9467bd #8c564b #e377c2 #17becf].freeze

    module_function

    def line_chart(path:, title:, x_label:, y_label:, series:, width: 760, height: 420)
      xs = series.flat_map { |s| s[:x] }
      ys = series.flat_map { |s| s[:y] }
      x0, x1 = padded_bounds(xs)
      y0, y1 = padded_bounds(ys)
      lines = svg_header(width, height, title, x_label, y_label)
      add_grid_and_ticks(lines, x0, x1, y0, y1, width, height)

      series.each_with_index do |s, idx|
        color = s[:color] || PALETTE[idx % PALETTE.length]
        pts = s[:x].each_with_index.map do |xv, i|
          "#{map_x(xv, x0, x1, width).round(1)},#{map_y(s[:y][i], y0, y1, height).round(1)}"
        end.join(' ')
        dash = s[:dashed] ? ' stroke-dasharray="6,4"' : ''
        lines << "<polyline fill=\"none\" stroke=\"#{color}\" stroke-width=\"2.2\"#{dash} points=\"#{pts}\"/>"
        s[:x].each_with_index do |xv, i|
          cx = map_x(xv, x0, x1, width).round(1)
          cy = map_y(s[:y][i], y0, y1, height).round(1)
          lines << "<circle cx=\"#{cx}\" cy=\"#{cy}\" r=\"2.8\" fill=\"#{color}\"/>"
        end
      end

      add_legend(lines, series, width)
      write_svg(path, lines)
    end

    def bar_chart(path:, title:, y_label:, labels:, values:, color: '#1f77b4', width: 780, height: 420)
      y0 = [values.min, 0.0].min
      y1 = [values.max * 1.12, 1.0e-6].max
      lines = svg_header(width, height, title, '', y_label)
      add_y_ticks(lines, y0, y1, width, height)

      plot_w = width - 110
      slot = plot_w.to_f / [labels.length, 1].max
      bar_w = [slot * 0.62, 54.0].min
      y_zero = map_y(0.0, y0, y1, height)

      labels.each_with_index do |lbl, i|
        bx = 75 + (i * slot) + ((slot - bar_w) * 0.5)
        yv = map_y(values[i], y0, y1, height)
        top = [yv, y_zero].min
        bh = [(yv - y_zero).abs, 1.0].max
        lines << "<rect x=\"#{bx.round(1)}\" y=\"#{top.round(1)}\" width=\"#{bar_w.round(1)}\" " \
                 "height=\"#{bh.round(1)}\" fill=\"#{color}\" opacity=\"0.88\"/>"
        tx = (bx + (bar_w * 0.5)).round(1)
        lines << "<text x=\"#{tx}\" y=\"#{height - 22}\" font-size=\"10\" text-anchor=\"middle\" " \
                 "fill=\"#333\">#{esc(lbl.to_s[0, 14])}</text>"
      end

      write_svg(path, lines)
    end

    def svg_header(width, height, title, x_label, y_label)
      [
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"#{width}\" height=\"#{height}\" " \
        "viewBox=\"0 0 #{width} #{height}\">",
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        "<text x=\"#{width / 2}\" y=\"26\" font-size=\"15\" font-weight=\"bold\" " \
        "text-anchor=\"middle\" fill=\"#111\">#{esc(title)}</text>",
        "<text x=\"#{width / 2}\" y=\"#{height - 6}\" font-size=\"11\" text-anchor=\"middle\" " \
        "fill=\"#444\">#{esc(x_label)}</text>",
        "<text x=\"16\" y=\"#{height / 2}\" font-size=\"11\" text-anchor=\"middle\" fill=\"#444\" " \
        "transform=\"rotate(-90 16 #{height / 2})\">#{esc(y_label)}</text>"
      ]
    end

    def add_grid_and_ticks(lines, x0, x1, y0, y1, width, height)
      add_y_ticks(lines, y0, y1, width, height)
      5.times do |i|
        xv = x0 + ((x1 - x0) * i / 4.0)
        px = map_x(xv, x0, x1, width).round(1)
        lines << "<line x1=\"#{px}\" y1=\"45\" x2=\"#{px}\" y2=\"#{height - 45}\" stroke=\"#eee\"/>"
        lines << "<text x=\"#{px}\" y=\"#{height - 28}\" font-size=\"10\" text-anchor=\"middle\" " \
                 "fill=\"#555\">#{format_tick(xv)}</text>"
      end
      lines << "<rect x=\"75\" y=\"45\" width=\"#{width - 110}\" height=\"#{height - 90}\" " \
               'fill="none" stroke="#999" stroke-width="1"/>'
    end

    def add_y_ticks(lines, y0, y1, width, height)
      5.times do |i|
        yv = y0 + ((y1 - y0) * i / 4.0)
        py = map_y(yv, y0, y1, height).round(1)
        lines << "<line x1=\"75\" y1=\"#{py}\" x2=\"#{width - 35}\" y2=\"#{py}\" stroke=\"#eee\"/>"
        lines << "<text x=\"70\" y=\"#{(py + 3).round(1)}\" font-size=\"10\" text-anchor=\"end\" " \
                 "fill=\"#555\">#{format_tick(yv)}</text>"
      end
      lines << "<rect x=\"75\" y=\"45\" width=\"#{width - 110}\" height=\"#{height - 90}\" " \
               'fill="none" stroke="#999" stroke-width="1"/>'
    end

    def add_legend(lines, series, width)
      return if series.length <= 1

      series.each_with_index do |s, idx|
        color = s[:color] || PALETTE[idx % PALETTE.length]
        y = 58 + (idx * 16)
        x = width - 190
        lines << "<rect x=\"#{x}\" y=\"#{y - 9}\" width=\"12\" height=\"10\" fill=\"#{color}\"/>"
        lines << "<text x=\"#{x + 16}\" y=\"#{y}\" font-size=\"10\" fill=\"#222\">#{esc(s[:name])}</text>"
      end
    end

    def padded_bounds(vals)
      lo = vals.min
      hi = vals.max
      span = hi - lo
      pad = span <= 1e-12 ? 0.5 : span * 0.08
      [lo - pad, hi + pad]
    end

    def map_x(val, x0, x1, width) = 75.0 + (((val - x0) / [x1 - x0, 1.0e-15].max) * (width - 110.0))

    def map_y(val, y0, y1, height) = (height - 45.0) - (((val - y0) / [y1 - y0, 1.0e-15].max) * (height - 90.0))

    def format_tick(val)
      return '0' if val.abs < 1e-12
      return format('%.3g', val) if val.abs < 0.01 || val.abs >= 1000.0

      format('%.3f', val).sub(/0+\z/, '').sub(/\.\z/, '')
    end

    def esc(str) = str.to_s.gsub('&', '&amp;').gsub('<', '&lt;').gsub('>', '&gt;')

    def write_svg(path, lines)
      FileUtils.mkdir_p(File.dirname(path))
      File.write(path, "#{lines.join("\n")}\n</svg>\n")
    end
  end
end
