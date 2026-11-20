# frozen_string_literal: true

module AsianOptions
  module Svg
    COLORS = ['#1f6feb', '#d1495b', '#2a9d8f'].freeze
    SERIES_NAMES = ['Naive MC', 'Antithetic', 'Control Variate'].freeze
    X_TICKS = [1_000, 10_000, 100_000].freeze

    class ConvergenceChart
      DIMS = { w: 900, h: 520, pad_l: 85, pad_r: 30, pad_t: 55, pad_b: 65 }.freeze

      def initialize(points, slopes)
        @points = points
        @slopes = slopes
        @plot_w = DIMS[:w] - DIMS[:pad_l] - DIMS[:pad_r]
        @plot_h = DIMS[:h] - DIMS[:pad_t] - DIMS[:pad_b]
        @bounds = compute_bounds
      end

      def render
        lines = [
          Svg.header(DIMS[:w], DIMS[:h]),
          %(<rect width="#{DIMS[:w]}" height="#{DIMS[:h]}" fill="#ffffff"/>),
          Svg.title(DIMS[:w] / 2, 30, 'Konvergensi Standard Error vs Jumlah Path (Log-Log, Ruby)')
        ]
        lines.concat(x_grid_lines, y_grid_lines)
        lines << Svg.frame_rect(DIMS[:pad_l], DIMS[:pad_t], @plot_w, @plot_h)
        lines.concat(series_lines, legend_lines, axis_labels)
        lines << '</svg>'
        "#{lines.join("\n")}\n"
      end

      private

      def compute_bounds
        {
          min_x: Math.log10(900.0), max_x: Math.log10(230_000.0),
          min_y: Math.log10(@points.map(&:se_control_variate).min * 0.7),
          max_y: Math.log10(@points.map(&:se_naive).max * 1.35)
        }
      end

      def map_x(n_val)
        DIMS[:pad_l] + (((Math.log10(n_val.to_f) - @bounds[:min_x]) / (@bounds[:max_x] - @bounds[:min_x])) * @plot_w)
      end

      def map_y(se_val)
        ratio = (Math.log10(se_val) - @bounds[:min_y]) / (@bounds[:max_y] - @bounds[:min_y])
        DIMS[:pad_t] + ((1.0 - ratio) * @plot_h)
      end

      def x_grid_lines
        X_TICKS.flat_map do |tick|
          x = map_x(tick)
          [
            format('<line x1="%<x>.1f" y1="%<y1>d" x2="%<x>.1f" y2="%<y2>d" stroke="#e5e7eb" stroke-dasharray="4,4"/>',
                   x: x, y1: DIMS[:pad_t], y2: DIMS[:pad_t] + @plot_h),
            format('<text x="%<x>.1f" y="%<y>d" text-anchor="middle" font-size="12" fill="#374151">%<t>d</text>',
                   x: x, y: DIMS[:pad_t] + @plot_h + 22, t: tick)
          ]
        end
      end

      def y_grid_lines
        (@bounds[:min_y].ceil..@bounds[:max_y].floor).flat_map do |exp_i|
          val = 10.0**exp_i
          y = map_y(val)
          [
            format('<line x1="%<x1>d" y1="%<y>.1f" x2="%<x2>d" y2="%<y>.1f" stroke="#e5e7eb" stroke-dasharray="4,4"/>',
                   x1: DIMS[:pad_l], x2: DIMS[:pad_l] + @plot_w, y: y),
            format('<text x="%<x>d" y="%<y>.1f" text-anchor="end" font-size="12" fill="#374151">%<v>.0e</text>',
                   x: DIMS[:pad_l] - 10, y: y + 4.0, v: val)
          ]
        end
      end

      def series_lines
        %i[se_naive se_antithetic se_control_variate].each_with_index.flat_map do |attr, idx|
          color = COLORS[idx]
          coords = @points.map { |pt| [map_x(pt.n), map_y(pt.public_send(attr))] }
          pts_str = coords.map { |x, y| format('%<x>.1f,%<y>.1f', x: x, y: y) }.join(' ')
          [%(<polyline fill="none" stroke="#{color}" stroke-width="2.5" points="#{pts_str}"/>)] +
            coords.map { |x, y| format('<circle cx="%<x>.1f" cy="%<y>.1f" r="4" fill="%<c>s"/>', x: x, y: y, c: color) }
        end
      end

      def legend_lines
        slope_vals = [@slopes.naive, @slopes.antithetic, @slopes.control_variate]
        box_x = DIMS[:pad_l] + @plot_w - 252
        box_y = DIMS[:pad_t] + 8
        lines = [
          "<rect x=\"#{box_x}\" y=\"#{box_y}\" width=\"244\" height=\"72\" rx=\"6\" " \
          'fill="#ffffff" fill-opacity="0.9" stroke="#d1d5db"/>'
        ]
        SERIES_NAMES.each_with_index do |name, idx|
          ly = DIMS[:pad_t] + 22 + (idx * 22)
          lx = DIMS[:pad_l] + @plot_w - 240
          col = COLORS[idx]
          lines << %(<line x1="#{lx}" y1="#{ly}" x2="#{lx + 24}" y2="#{ly}" stroke="#{col}" stroke-width="2.5"/>)
          lines << format('<text x="%<x>d" y="%<y>d" font-size="12" fill="#111827">%<n>s (slope=%<s>.3f)</text>',
                          x: lx + 32, y: ly + 4, n: name, s: slope_vals[idx])
        end
        lines
      end

      def axis_labels
        mid_x = DIMS[:pad_l] + (@plot_w / 2)
        mid_y = DIMS[:pad_t] + (@plot_h / 2)
        [
          "<text x=\"#{mid_x}\" y=\"#{DIMS[:h] - 16}\" text-anchor=\"middle\" " \
          'font-size="13" fill="#111827">Jumlah Path N (skala log)</text>',
          "<text x=\"22\" y=\"#{mid_y}\" text-anchor=\"middle\" font-size=\"13\" " \
          "fill=\"#111827\" transform=\"rotate(-90 22 #{mid_y})\">Standard Error (SE)</text>"
        ]
      end
    end
  end
end
