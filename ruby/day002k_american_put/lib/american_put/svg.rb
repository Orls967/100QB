# frozen_string_literal: true

module AmericanPut
  module Svg
    PAD = { left: 82, right: 32, top: 54, bottom: 64 }.freeze
    WIDTH = 860
    HEIGHT = 500
    PLOT_W = WIDTH - PAD[:left] - PAD[:right]
    PLOT_H = HEIGHT - PAD[:top] - PAD[:bottom]

    module_function

    def write_convergence_svg(path, crr_sweep, psor_sweep, lsm_sweep, ref_price)
      bounds = convergence_bounds(crr_sweep, psor_sweep, lsm_sweep, ref_price)
      lines = svg_header('Konvergensi Harga Put Amerika: CRR, PSOR PDE & LSM')
      lines += draw_grid_and_axes(bounds, 'Harga Put Amerika (USD)')
      lines << draw_ref_line(ref_price, bounds)
      lines += draw_series(crr_sweep.map(&:american_price), bounds, '#1d4ed8')
      lines += draw_series(psor_sweep.map(&:american_price), bounds, '#059669')
      lines += draw_lsm_series(lsm_sweep, bounds, '#d97706')
      lines += draw_convergence_legend(PAD[:left] + 16, PAD[:top] + 18)
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def write_spot_profile_svg(path, spot_points)
      all_vals = spot_points.flat_map { |pt| [pt.intrinsic, pt.bs_european_put, pt.crr_american, pt.psor_american] }
      bounds = [0.0, all_vals.max * 1.08]
      lines = svg_header('Profil Harga Put Amerika vs Eropa & Nilai Intrinsik (S0 = 80..120)')
      lines += draw_grid_and_axes(bounds, 'Nilai Opsi (USD)')
      lines += draw_spot_x_labels(spot_points)
      lines += draw_series(spot_points.map(&:intrinsic), bounds, '#64748b', dashed: true)
      lines += draw_series(spot_points.map(&:bs_european_put), bounds, '#dc2626')
      lines += draw_series(spot_points.map(&:crr_american), bounds, '#1d4ed8')
      lines += draw_series(spot_points.map(&:psor_american), bounds, '#059669', dashed: true)
      lines += draw_spot_legend(PAD[:left] + PLOT_W - 210, PAD[:top] + 18)
      lines << '</svg>'
      File.write(path, "#{lines.join("\n")}\n")
    end

    def convergence_bounds(crr_sweep, psor_sweep, lsm_sweep, ref_price)
      vals = crr_sweep.map(&:american_price) + psor_sweep.map(&:american_price) +
             lsm_sweep.map(&:ci_low) + lsm_sweep.map(&:ci_high) + [ref_price]
      [vals.min - 0.02, vals.max + 0.02]
    end
    private_class_method :convergence_bounds

    def svg_header(title)
      [
        format('<svg xmlns="http://www.w3.org/2000/svg" width="%<w>d" height="%<h>d" viewBox="0 0 %<w>d %<h>d">',
               w: WIDTH, h: HEIGHT),
        %(  <rect width="#{WIDTH}" height="#{HEIGHT}" fill="#ffffff"/>),
        format('  <text x="%<x>d" y="30" text-anchor="middle" font-family="sans-serif" ' \
               'font-size="16" font-weight="bold" fill="#0f172a">%<t>s</text>', x: WIDTH / 2, t: title)
      ]
    end
    private_class_method :svg_header

    def draw_grid_and_axes(bounds, y_label)
      out = [format('  <rect x="%<x>d" y="%<y>d" width="%<w>d" height="%<h>d" fill="#f8fafc" stroke="#cbd5e1"/>',
                    x: PAD[:left], y: PAD[:top], w: PLOT_W, h: PLOT_H)]
      5.times { |idx| out.concat(grid_tick(idx, bounds)) }
      yc = PAD[:top] + (PLOT_H / 2)
      out << format('  <text x="20" y="%<yc>d" text-anchor="middle" transform="rotate(-90 20 %<yc>d)" ' \
                    'font-family="sans-serif" font-size="12" fill="#334155">%<lbl>s</text>', yc: yc, lbl: y_label)
      out
    end
    private_class_method :draw_grid_and_axes

    def grid_tick(idx, bounds)
      frac = idx / 4.0
      y_px = (PAD[:top] + (PLOT_H * (1.0 - frac))).round(2)
      val = bounds[0] + ((bounds[1] - bounds[0]) * frac)
      [
        %(  <line x1="#{PAD[:left]}" y1="#{y_px}" x2="#{PAD[:left] + PLOT_W}" ) +
          %(y2="#{y_px}" stroke="#e2e8f0" stroke-dasharray="4,4"/>),
        %(  <text x="#{PAD[:left] - 8}" y="#{(y_px + 4).round(2)}" text-anchor="end" ) +
          %(font-family="monospace" font-size="11" fill="#475569">#{format('%.3f', val)}</text>)
      ]
    end
    private_class_method :grid_tick

    def draw_ref_line(ref_price, bounds)
      y_px = map_y(ref_price, bounds).round(2)
      %(  <line x1="#{PAD[:left]}" y1="#{y_px}" x2="#{PAD[:left] + PLOT_W}" ) +
        %(y2="#{y_px}" stroke="#0f172a" stroke-width="1.2" stroke-dasharray="6,4"/>)
    end
    private_class_method :draw_ref_line

    def draw_series(values, bounds, color, dashed: false)
      return [] if values.empty?

      coords = values.each_with_index.map do |val, idx|
        x_px = PAD[:left] + (PLOT_W * idx / [values.length - 1, 1].max.to_f)
        [x_px.round(2), map_y(val, bounds).round(2)]
      end
      pts = coords.map { |x_px, y_px| "#{x_px},#{y_px}" }.join(' ')
      dash = dashed ? ' stroke-dasharray="6,3"' : ''
      out = [%(  <polyline fill="none" stroke="#{color}" stroke-width="2.2"#{dash} points="#{pts}"/>)]
      coords.each { |x_px, y_px| out << %(  <circle cx="#{x_px}" cy="#{y_px}" r="3.5" fill="#{color}"/>) }
      out
    end
    private_class_method :draw_series

    def draw_lsm_series(lsm_sweep, bounds, color)
      out = draw_series(lsm_sweep.map(&:price), bounds, color)
      lsm_sweep.each_with_index do |pt, idx|
        x_px = (PAD[:left] + (PLOT_W * idx / [lsm_sweep.length - 1, 1].max.to_f)).round(2)
        y_lo = map_y(pt.ci_low, bounds).round(2)
        y_hi = map_y(pt.ci_high, bounds).round(2)
        out << %(  <line x1="#{x_px}" y1="#{y_lo}" x2="#{x_px}" y2="#{y_hi}" stroke="#{color}" stroke-width="1.5"/>)
      end
      out
    end
    private_class_method :draw_lsm_series

    def draw_spot_x_labels(spot_points)
      y_px = PAD[:top] + PLOT_H + 20
      spot_points.each_with_index.map do |pt, idx|
        x_px = (PAD[:left] + (PLOT_W * idx / [spot_points.length - 1, 1].max.to_f)).round(2)
        %(  <text x="#{x_px}" y="#{y_px}" text-anchor="middle" font-family="monospace" ) +
          %(font-size="11" fill="#334155">S0=#{pt.s0.to_i}</text>)
      end
    end
    private_class_method :draw_spot_x_labels

    def draw_convergence_legend(pos_x, pos_y)
      legend_box(pos_x, pos_y, 210, [
                   ['#1d4ed8', '● CRR Binomial (n=50..4000)'], ['#059669', '● PSOR PDE (100..800)'],
                   ['#d97706', '● LSM MC (5k..100k, CI 95%)'], ['#0f172a', '-- Referensi CRR-4000']
                 ])
    end
    private_class_method :draw_convergence_legend

    def draw_spot_legend(pos_x, pos_y)
      legend_box(pos_x, pos_y, 205, [
                   ['#1d4ed8', '● CRR Amerika'], ['#059669', '-- PSOR Amerika'],
                   ['#dc2626', '● BS Put Eropa'], ['#64748b', '-- Intrinsik max(K-S, 0)']
                 ])
    end
    private_class_method :draw_spot_legend

    def legend_box(pos_x, pos_y, box_w, items)
      out = [format('  <rect x="%<x>d" y="%<y>d" width="%<w>d" height="78" ' \
                    'fill="#ffffff" fill-opacity="0.9" stroke="#cbd5e1"/>',
                    x: pos_x - 8, y: pos_y - 12, w: box_w)]
      items.each_with_index do |(color, text), idx|
        out << format('  <text x="%<x>d" y="%<y>d" font-family="sans-serif" font-size="12" fill="%<c>s">%<t>s</text>',
                      x: pos_x, y: pos_y + 4 + (idx * 18), c: color, t: text)
      end
      out
    end
    private_class_method :legend_box

    def map_y(val, bounds)
      span = bounds[1] - bounds[0]
      frac = span.zero? ? 0.5 : (val - bounds[0]) / span
      PAD[:top] + (PLOT_H * (1.0 - frac))
    end
    private_class_method :map_y
  end
end
