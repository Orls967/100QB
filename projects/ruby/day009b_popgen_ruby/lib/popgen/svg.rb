# frozen_string_literal: true

module PopGen
  # Hand-written pure-Ruby XML SVG generator (no external plotting libraries).
  module SVG
    module_function

    def esc(str)
      str.to_s.gsub('&', '&amp;').gsub('<', '&lt;').gsub('>', '&gt;').gsub('"', '&quot;')
    end

    def header(w, h, title, subtitle)
      lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        %(<svg xmlns="http://www.w3.org/2000/svg" width="#{w}" height="#{h}" viewBox="0 0 #{w} #{h}">),
        %(  <rect x="0" y="0" width="#{w}" height="#{h}" fill="#fafbfc"/>),
        %(  <rect x="1" y="1" width="#{w - 2}" height="#{h - 2}" fill="none" stroke="#d8dee9" stroke-width="1"/>),
        %(  <text x="#{w / 2}" y="28" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="15" font-weight="bold" fill="#1f2937">#{esc(title)}</text>)
      ]
      if subtitle && !subtitle.empty?
        lines << %(  <text x="#{w / 2}" y="46" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="11" fill="#4b5563">#{esc(subtitle)}</text>)
      end
      lines
    end

    def safe_tx(v, is_log)
      return v unless is_log

      Math.log10(v > 1e-30 ? v : 1e-30)
    end

    # Multi-series XY plot with optional SE error bars, reference line, and log axes.
    # series: Array of { label:, color:, x:, y:, yerr: nil, line: true, points: true, dashed: false }
    def write_xy(path, title:, subtitle: '', xlabel:, ylabel:, series:,
                 width: 760, height: 460, log_x: false, log_y: false,
                 ref_y: nil, ref_label: nil)
      ml = 78
      mr = 32
      mt = 62
      mb = 60
      pw = width - ml - mr
      ph = height - mt - mb

      xs_all = []
      ys_all = []
      series.each do |sr|
        sr[:x].each_with_index do |xv, i|
          yv = sr[:y][i]
          next unless xv.finite? && yv.finite?
          next if log_x && xv <= 0.0
          next if log_y && yv <= 0.0

          xs_all << safe_tx(xv, log_x)
          ty = safe_tx(yv, log_y)
          ys_all << ty
          if sr[:yerr] && !log_y && sr[:yerr][i] && sr[:yerr][i].finite?
            e = sr[:yerr][i].abs
            ys_all << (ty - e)
            ys_all << (ty + e)
          end
        end
      end
      ys_all << safe_tx(ref_y, log_y) if ref_y && (!log_y || ref_y > 0.0)

      xmin, xmax = xs_all.empty? ? [0.0, 1.0] : [xs_all.min, xs_all.max]
      ymin, ymax = ys_all.empty? ? [0.0, 1.0] : [ys_all.min, ys_all.max]
      if xmin >= xmax
        xmin -= 1.0
        xmax += 1.0
      end
      if ymin >= ymax
        ymin -= 1.0
        ymax += 1.0
      end
      xpad = 0.04 * (xmax - xmin)
      ypad = 0.08 * (ymax - ymin)
      xmin -= xpad
      xmax += xpad
      ymin -= ypad
      ymax += ypad

      out = header(width, height, title, subtitle)
      out << %(  <rect x="#{ml}" y="#{mt}" width="#{pw}" height="#{ph}" fill="#ffffff" stroke="#9ca3af" stroke-width="1"/>)

      6.times do |t|
        frac = t / 5.0
        gx = ml + frac * pw
        vx = xmin + frac * (xmax - xmin)
        dx = log_x ? (10.0**vx) : vx
        out << format('  <line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e5e7eb" stroke-width="1"/>', gx, mt, gx, mt + ph)
        out << format('  <text x="%.1f" y="%d" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%.3g</text>', gx, mt + ph + 16, dx)

        gy = mt + ph - frac * ph
        vy = ymin + frac * (ymax - ymin)
        dy = log_y ? (10.0**vy) : vy
        out << format('  <line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e7eb" stroke-width="1"/>', ml, gy, ml + pw, gy)
        out << format('  <text x="%d" y="%.1f" text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%.3g</text>', ml - 8, gy + 3.5, dy)
      end

      out << %(  <text x="#{ml + pw / 2}" y="#{height - 14}" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#111827">#{esc(xlabel)}</text>)
      out << %(  <text x="18" y="#{mt + ph / 2}" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#111827" transform="rotate(-90 18 #{mt + ph / 2})">#{esc(ylabel)}</text>)

      if ref_y && (!log_y || ref_y > 0.0)
        py = mt + ph - ((safe_tx(ref_y, log_y) - ymin) / (ymax - ymin)) * ph
        out << format('  <line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#6b7280" stroke-width="1.2" stroke-dasharray="5,4"/>', ml, py, ml + pw, py)
        if ref_label && !ref_label.empty?
          out << format('  <text x="%d" y="%.1f" text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#4b5563">%s</text>', ml + pw - 6, py - 4.0, esc(ref_label))
        end
      end

      series.each do |sr|
        col = sr[:color] || '#2563eb'
        if sr.fetch(:line, true) && sr[:x].length >= 2
          pts = []
          sr[:x].each_with_index do |xv, i|
            yv = sr[:y][i]
            next unless xv.finite? && yv.finite?
            next if log_x && xv <= 0.0
            next if log_y && yv <= 0.0

            px = ml + ((safe_tx(xv, log_x) - xmin) / (xmax - xmin)) * pw
            py = mt + ph - ((safe_tx(yv, log_y) - ymin) / (ymax - ymin)) * ph
            pts << format('%.1f,%.1f', px, py)
          end
          dash = sr[:dashed] ? ' stroke-dasharray="6,4"' : ''
          out << %(  <polyline fill="none" stroke="#{col}" stroke-width="2"#{dash} points="#{pts.join(' ')}"/>) unless pts.empty?
        end

        if sr[:yerr] && !log_y
          sr[:x].each_with_index do |xv, i|
            yv = sr[:y][i]
            ev = sr[:yerr][i]&.abs
            next unless xv.finite? && yv.finite? && ev && ev.finite? && ev > 0.0
            next if log_x && xv <= 0.0

            px = ml + ((safe_tx(xv, log_x) - xmin) / (xmax - xmin)) * pw
            py1 = mt + ph - ((yv - ev - ymin) / (ymax - ymin)) * ph
            py2 = mt + ph - ((yv + ev - ymin) / (ymax - ymin)) * ph
            out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>', px, py1, px, py2, col)
            out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>', px - 3.0, py1, px + 3.0, py1, col)
            out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.2"/>', px - 3.0, py2, px + 3.0, py2, col)
          end
        end

        next unless sr.fetch(:points, true)

        sr[:x].each_with_index do |xv, i|
          yv = sr[:y][i]
          next unless xv.finite? && yv.finite?
          next if log_x && xv <= 0.0
          next if log_y && yv <= 0.0

          px = ml + ((safe_tx(xv, log_x) - xmin) / (xmax - xmin)) * pw
          py = mt + ph - ((safe_tx(yv, log_y) - ymin) / (ymax - ymin)) * ph
          out << format('  <circle cx="%.1f" cy="%.1f" r="3.2" fill="%s" stroke="#ffffff" stroke-width="0.8"/>', px, py, col)
        end
      end

      unless series.empty?
        lx = ml + pw - 206
        ly = mt + 10
        lh = 16 * series.length + 10
        out << %(  <rect x="#{lx}" y="#{ly}" width="198" height="#{lh}" rx="4" fill="#ffffff" fill-opacity="0.9" stroke="#d1d5db"/>)
        series.each_with_index do |sr, s_idx|
          col = sr[:color] || '#2563eb'
          cy = ly + 15 + s_idx * 16
          dash = sr[:dashed] ? ' stroke-dasharray="4,3"' : ''
          out << %(  <line x1="#{lx + 8}" y1="#{cy - 3}" x2="#{lx + 26}" y2="#{cy - 3}" stroke="#{col}" stroke-width="2.2"#{dash}/>)
          out << %(  <text x="#{lx + 32}" y="#{cy}" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#1f2937">#{esc(sr[:label])}</text>)
        end
      end

      out << '</svg>'
      File.write(path, "#{out.join("\n")}\n")
      true
    end

    def write_bars(path, title:, subtitle: '', ylabel:, labels:, vals:, errs: nil, colors: nil, log_y: false)
      w = 760
      h = 460
      ml = 80
      mr = 30
      mt = 62
      mb = 74
      pw = w - ml - mr
      ph = h - mt - mb

      if log_y
        lvs = vals.map { |v| safe_tx(v > 1e-12 ? v : 1e-12, true) }
        ymin = lvs.min - 0.3
        ymax = lvs.max + 0.3
      else
        lows = vals.each_with_index.map { |v, i| v - (errs ? errs[i].to_f.abs : 0.0) }
        highs = vals.each_with_index.map { |v, i| v + (errs ? errs[i].to_f.abs : 0.0) }
        ymin = [0.0, lows.min].min
        ymax = highs.max
        ymax = ymin + 1.0 if ymin >= ymax
        ymax += 0.14 * (ymax - ymin)
      end

      out = header(w, h, title, subtitle)
      out << %(  <rect x="#{ml}" y="#{mt}" width="#{pw}" height="#{ph}" fill="#ffffff" stroke="#9ca3af" stroke-width="1"/>)

      6.times do |t|
        frac = t / 5.0
        gy = mt + ph - frac * ph
        vy = ymin + frac * (ymax - ymin)
        dy = log_y ? (10.0**vy) : vy
        out << format('  <line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e7eb"/>', ml, gy, ml + pw, gy)
        out << format('  <text x="%d" y="%.1f" text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%.3g</text>', ml - 8, gy + 3.5, dy)
      end

      n_bars = [labels.length, 1].max
      slot = pw.to_f / n_bars
      bw = slot * 0.56
      y0_px = log_y ? (mt + ph) : [mt + ph - ((0.0 - ymin) / (ymax - ymin)) * ph, mt + ph].min

      labels.each_with_index do |lbl, i|
        cx = ml + (i + 0.5) * slot
        tv = log_y ? safe_tx(vals[i] > 1e-12 ? vals[i] : 1e-12, true) : vals[i]
        py = mt + ph - ((tv - ymin) / (ymax - ymin)) * ph
        top = [py, y0_px].min
        bh = (y0_px - py).abs
        col = colors ? colors[i] : '#2563eb'
        out << format('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" fill-opacity="0.85" stroke="#1f2937" stroke-width="0.8"/>', cx - 0.5 * bw, top, bw, bh, col)

        if errs && errs[i] && errs[i] > 0.0 && !log_y
          py1 = mt + ph - ((vals[i] - errs[i] - ymin) / (ymax - ymin)) * ph
          py2 = mt + ph - ((vals[i] + errs[i] - ymin) / (ymax - ymin)) * ph
          out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#111827" stroke-width="1.4"/>', cx, py1, cx, py2)
          out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#111827" stroke-width="1.4"/>', cx - 5.0, py1, cx + 5.0, py1)
          out << format('  <line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#111827" stroke-width="1.4"/>', cx - 5.0, py2, cx + 5.0, py2)
        end

        out << format('  <text x="%.1f" y="%.1f" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="9.5" fill="#111827">%.3g</text>', cx, top - 5.0, vals[i])
        out << format('  <text x="%.1f" y="%d" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%s</text>', cx, mt + ph + 18, esc(lbl))
      end

      out << %(  <text x="18" y="#{mt + ph / 2}" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#111827" transform="rotate(-90 18 #{mt + ph / 2})">#{esc(ylabel)}</text>)
      out << '</svg>'
      File.write(path, "#{out.join("\n")}\n")
      true
    end

    def write_hist(path, title:, subtitle: '', xlabel:, ylabel:, centers:, counts:, theory_y: nil, theory_label: nil)
      w = 760
      h = 460
      ml = 76
      mr = 32
      mt = 62
      mb = 60
      pw = w - ml - mr
      ph = h - mt - mb

      dx = centers.length > 1 ? (centers[1] - centers[0]) : 1.0
      xmin = centers.first - 0.5 * dx
      xmax = centers.last + 0.5 * dx
      ymax = ([counts.max || 0.0, theory_y ? (theory_y.max || 0.0) : 0.0].max)
      ymax = 1.0 if ymax <= 0.0
      ymax *= 1.14

      out = header(w, h, title, subtitle)
      out << %(  <rect x="#{ml}" y="#{mt}" width="#{pw}" height="#{ph}" fill="#ffffff" stroke="#9ca3af" stroke-width="1"/>)

      6.times do |t|
        frac = t / 5.0
        gx = ml + frac * pw
        gy = mt + ph - frac * ph
        out << format('  <line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e5e7eb"/>', gx, mt, gx, mt + ph)
        out << format('  <line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e5e7eb"/>', ml, gy, ml + pw, gy)
        out << format('  <text x="%.1f" y="%d" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%.3g</text>', gx, mt + ph + 16, xmin + frac * (xmax - xmin))
        out << format('  <text x="%d" y="%.1f" text-anchor="end" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#374151">%.3g</text>', ml - 8, gy + 3.5, frac * ymax)
      end

      bw = (dx / (xmax - xmin)) * pw * 0.90
      centers.each_with_index do |c, i|
        cx = ml + ((c - xmin) / (xmax - xmin)) * pw
        bh = (counts[i] / ymax) * ph
        out << format('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#60a5fa" fill-opacity="0.72" stroke="#2563eb" stroke-width="0.6"/>', cx - 0.5 * bw, mt + ph - bh, bw, bh)
      end

      if theory_y
        pts = centers.each_with_index.map do |c, i|
          cx = ml + ((c - xmin) / (xmax - xmin)) * pw
          py = mt + ph - (theory_y[i] / ymax) * ph
          format('%.1f,%.1f', cx, py)
        end
        out << %(  <polyline fill="none" stroke="#dc2626" stroke-width="2.2" points="#{pts.join(' ')}"/>)
        if theory_label
          out << %(  <rect x="#{ml + pw - 190}" y="#{mt + 10}" width="182" height="26" rx="4" fill="#ffffff" fill-opacity="0.9" stroke="#d1d5db"/>)
          out << %(  <line x1="#{ml + pw - 180}" y1="#{mt + 23}" x2="#{ml + pw - 160}" y2="#{mt + 23}" stroke="#dc2626" stroke-width="2.2"/>)
          out << %(  <text x="#{ml + pw - 154}" y="#{mt + 26}" font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#1f2937">#{esc(theory_label)}</text>)
        end
      end

      out << %(  <text x="#{ml + pw / 2}" y="#{h - 14}" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#111827">#{esc(xlabel)}</text>)
      out << %(  <text x="18" y="#{mt + ph / 2}" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="12" fill="#111827" transform="rotate(-90 18 #{mt + ph / 2})">#{esc(ylabel)}</text>)
      out << '</svg>'
      File.write(path, "#{out.join("\n")}\n")
      true
    end
  end
end
