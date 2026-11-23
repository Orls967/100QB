# frozen_string_literal: true

require 'fileutils'

module CdoTranche
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
      xs = series.flat_map { |s| s[:x] }
      ys = series.flat_map { |s| s[:y] }
      xmin, xmax = bounds(xs)
      ymin, ymax = bounds(ys)

      lines = [header(title, x_label, y_label, xmin, xmax, ymin, ymax)]
      series.each_with_index do |s, idx|
        pts = s[:x].each_with_index.map do |xv, i|
          "#{sx(xv, xmin, xmax).round(2)},#{sy(s[:y][i], ymin, ymax).round(2)}"
        end.join(' ')
        dash = s[:dashed] ? ' stroke-dasharray="6,4"' : ''
        lines << %(<polyline fill="none" stroke="#{s[:color]}" stroke-width="2.3"#{dash} points="#{pts}"/>)
        lx = PAD_L + 14 + (idx * 200)
        lines << %(<line x1="#{lx}" y1="30" x2="#{lx + 22}" y2="30" stroke="#{s[:color]}" stroke-width="2.5"#{dash}/>)
        lines << %(<text x="#{lx + 28}" y="34" font-size="11" fill="#1f2937">#{s[:label]}</text>)
      end
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def bounds(vals)
      clean = vals.compact
      lo = clean.min
      hi = clean.max
      return [lo - 1.0, hi + 1.0] if (hi - lo).abs < 1.0e-12

      pad = 0.06 * (hi - lo)
      [lo - pad, hi + pad]
    end
    private_class_method :bounds

    def sx(val, xmin, xmax) = PAD_L + ((val - xmin) / (xmax - xmin) * (WIDTH - PAD_L - PAD_R))
    private_class_method :sx

    def sy(val, ymin, ymax) = HEIGHT - PAD_B - ((val - ymin) / (ymax - ymin) * (HEIGHT - PAD_T - PAD_B))
    private_class_method :sy

    def header(title, x_label, y_label, xmin, xmax, ymin, ymax)
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
        %(<text x="#{PAD_L}" y="#{HEIGHT - 42}" font-size="10" fill="#6b7280">#{format('%.3g', xmin)}</text>),
        %(<text x="#{WIDTH - PAD_R}" y="#{HEIGHT - 42}" text-anchor="end" ) +
          %(font-size="10" fill="#6b7280">#{format('%.3g', xmax)}</text>),
        %(<text x="#{PAD_L - 6}" y="#{HEIGHT - PAD_B}" text-anchor="end" ) +
          %(font-size="10" fill="#6b7280">#{format('%.3g', ymin)}</text>),
        %(<text x="#{PAD_L - 6}" y="#{PAD_T + 8}" text-anchor="end" ) +
          %(font-size="10" fill="#6b7280">#{format('%.3g', ymax)}</text>)
      ].join("\n")
    end
    private_class_method :header
  end
end
