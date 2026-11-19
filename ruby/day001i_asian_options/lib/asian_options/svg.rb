# frozen_string_literal: true

module AsianOptions
  module Svg
    module_function

    def write_convergence_svg(path, points, slopes)
      raise ArgumentError, 'points tidak boleh kosong' if points.empty?

      File.write(path, ConvergenceChart.new(points, slopes).render)
    end

    def write_estimators_svg(path, estimators)
      raise ArgumentError, 'estimators tidak boleh kosong' if estimators.empty?

      File.write(path, EstimatorsChart.new(estimators).render)
    end

    def header(width, height)
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"#{width}\" height=\"#{height}\" " \
        "viewBox=\"0 0 #{width} #{height}\" font-family=\"system-ui, -apple-system, sans-serif\">"
    end

    def title(x_pos, y_pos, text)
      "<text x=\"#{x_pos}\" y=\"#{y_pos}\" text-anchor=\"middle\" " \
        "font-size=\"18\" font-weight=\"bold\" fill=\"#111827\">#{text}</text>"
    end

    def frame_rect(x_pos, y_pos, width, height)
      %(<rect x="#{x_pos}" y="#{y_pos}" width="#{width}" height="#{height}" fill="none" stroke="#9ca3af"/>)
    end

    class EstimatorsChart
      DIMS = { w: 900, h: 440, pad_l: 210, pad_r: 60, pad_t: 60, pad_b: 60 }.freeze

      def initialize(estimators)
        @estimators = estimators
        @plot_w = DIMS[:w] - DIMS[:pad_l] - DIMS[:pad_r]
        @plot_h = DIMS[:h] - DIMS[:pad_t] - DIMS[:pad_b]
        @min_x, @max_x = compute_bounds
      end

      def render
        lines = [
          Svg.header(DIMS[:w], DIMS[:h]),
          %(<rect width="#{DIMS[:w]}" height="#{DIMS[:h]}" fill="#ffffff"/>),
          Svg.title(DIMS[:w] / 2, 32, 'Perbandingan Estimasi Harga + CI 95% &amp; VRF (200k Path, Ruby)'),
          Svg.frame_rect(DIMS[:pad_l], DIMS[:pad_t], @plot_w, @plot_h)
        ]
        @estimators.each_with_index { |est, idx| lines.concat(row_lines(est, idx)) }
        lines << x_axis_label
        lines << '</svg>'
        "#{lines.join("\n")}\n"
      end

      private

      def compute_bounds
        min_x = @estimators.map(&:ci_low).min
        max_x = @estimators.map(&:ci_high).max
        margin = (max_x - min_x) * 0.15
        [min_x - margin, max_x + margin]
      end

      def map_x(val)
        DIMS[:pad_l] + (((val - @min_x) / (@max_x - @min_x)) * @plot_w)
      end

      def x_axis_label
        mid_x = DIMS[:pad_l] + (@plot_w / 2)
        "<text x=\"#{mid_x}\" y=\"#{DIMS[:h] - 18}\" text-anchor=\"middle\" " \
          'font-size="13" fill="#111827">Harga Opsi Asia</text>'
      end

      def row_lines(est, idx)
        cy = DIMS[:pad_t] + (((idx + 0.5) / @estimators.length) * @plot_h)
        x_lo = map_x(est.ci_low)
        x_hi = map_x(est.ci_high)
        x_mid = map_x(est.price)
        color = COLORS[idx % COLORS.length]

        [
          format('<text x="%<x>d" y="%<y>.1f" text-anchor="end" font-size="12" fill="#111827">%<n>s</text>',
                 x: DIMS[:pad_l] - 12, y: cy + 4.0, n: est.name),
          format('<line x1="%<x1>.1f" y1="%<y>.1f" x2="%<x2>.1f" y2="%<y>.1f" stroke="%<c>s" stroke-width="3"/>',
                 x1: x_lo, x2: x_hi, y: cy, c: color),
          cap_line(x_lo, cy, color),
          cap_line(x_hi, cy, color),
          format('<circle cx="%<x>.1f" cy="%<y>.1f" r="5" fill="%<c>s"/>', x: x_mid, y: cy, c: color),
          format('<text x="%<x>.1f" y="%<y>.1f" font-size="11" fill="#374151">%<p>.4f (VRF %<v>.1fx)</text>',
                 x: x_hi + 8.0, y: cy + 4.0, p: est.price, v: est.vrf)
        ]
      end

      def cap_line(x_pos, y_pos, color)
        format('<line x1="%<x>.1f" y1="%<y1>.1f" x2="%<x>.1f" y2="%<y2>.1f" stroke="%<c>s" stroke-width="2"/>',
               x: x_pos, y1: y_pos - 6.0, y2: y_pos + 6.0, c: color)
      end
    end
  end
end
