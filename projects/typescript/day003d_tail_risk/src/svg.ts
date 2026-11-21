/**
 * Pure TypeScript SVG chart generator (no external plotting dependencies)
 * producing 7 publication-quality vector figures for day003d.
 */

export interface LineSeries {
  label: string;
  color: string;
  x: number[];
  y: number[];
  dashed?: boolean;
  errorBars?: number[];
}

export interface BarGroup {
  category: string;
  values: Array<{ label: string; value: number; color: string }>;
}

function escapeXml(s: string): string {
  return s
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function niceTicks(minVal: number, maxVal: number, count = 6): number[] {
  const span = maxVal - minVal;
  if (!(span > 0.0) || !Number.isFinite(span)) {
    return [minVal - 1, minVal, minVal + 1];
  }
  const rawStep = span / Math.max(count - 1, 1);
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const niceMult = norm <= 1.2 ? 1 : norm <= 2.5 ? 2 : norm <= 5.5 ? 5 : 10;
  const step = niceMult * mag;
  const start = Math.floor(minVal / step) * step;
  const ticks: number[] = [];
  for (let t = start; t <= maxVal + 0.5 * step; t += step) {
    if (t >= minVal - 0.1 * step && t <= maxVal + 0.1 * step) {
      ticks.push(Math.abs(t) < 1e-12 * step ? 0.0 : t);
    }
  }
  if (ticks.length < 2) {
    return [minVal, 0.5 * (minVal + maxVal), maxVal];
  }
  return ticks;
}

function fmtTick(v: number): string {
  if (Math.abs(v) < 1e-12) return "0";
  if (Math.abs(v) >= 1000 || Math.abs(v) < 1e-2) return v.toExponential(1);
  if (Math.abs(v - Math.round(v)) < 1e-9) return String(Math.round(v));
  return v.toFixed(Math.abs(v) < 1 ? 3 : 2);
}

/** Generic clean multi-series line chart with optional horizontal reference line & error bars. */
export function renderLineChartSvg(opts: {
  title: string;
  subtitle?: string;
  xLabel: string;
  yLabel: string;
  series: LineSeries[];
  hLine?: { y: number; label: string; color: string };
  width?: number;
  height?: number;
}): string {
  const width = opts.width ?? 820;
  const height = opts.height ?? 500;
  const padL = 78;
  const padR = 32;
  const padT = 62;
  const padB = 64;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  let xMin = Infinity;
  let xMax = -Infinity;
  let yMin = opts.hLine ? opts.hLine.y : Infinity;
  let yMax = opts.hLine ? opts.hLine.y : -Infinity;

  for (const s of opts.series) {
    for (let i = 0; i < s.x.length; i += 1) {
      const xv = s.x[i] ?? 0.0;
      const yv = s.y[i] ?? 0.0;
      const eb = s.errorBars?.[i] ?? 0.0;
      if (xv < xMin) xMin = xv;
      if (xv > xMax) xMax = xv;
      if (yv - eb < yMin) yMin = yv - eb;
      if (yv + eb > yMax) yMax = yv + eb;
    }
  }
  if (!(xMax > xMin)) {
    xMin -= 1;
    xMax += 1;
  }
  const yPad = 0.08 * Math.max(yMax - yMin, 1e-6);
  yMin -= yPad;
  yMax += yPad;

  const mapX = (x: number): number => padL + ((x - xMin) / (xMax - xMin)) * plotW;
  const mapY = (y: number): number => padT + plotH - ((y - yMin) / (yMax - yMin)) * plotH;

  const xTicks = niceTicks(xMin, xMax, 6);
  const yTicks = niceTicks(yMin, yMax, 6);

  const parts: string[] = [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
    `<rect width="100%" height="100%" fill="#ffffff"/>`,
    `<text x="${width / 2}" y="26" text-anchor="middle" font-family="system-ui, -apple-system, sans-serif" font-size="16" font-weight="700" fill="#111827">${escapeXml(opts.title)}</text>`,
  ];
  if (opts.subtitle) {
    parts.push(
      `<text x="${width / 2}" y="46" text-anchor="middle" font-family="system-ui, -apple-system, sans-serif" font-size="12" fill="#4b5563">${escapeXml(opts.subtitle)}</text>`,
    );
  }

  for (const xt of xTicks) {
    const px = mapX(xt);
    parts.push(
      `<line x1="${px.toFixed(1)}" y1="${padT}" x2="${px.toFixed(1)}" y2="${padT + plotH}" stroke="#e5e7eb" stroke-width="1"/>`,
      `<text x="${px.toFixed(1)}" y="${padT + plotH + 18}" text-anchor="middle" font-family="system-ui, sans-serif" font-size="11" fill="#374151">${fmtTick(xt)}</text>`,
    );
  }
  for (const yt of yTicks) {
    const py = mapY(yt);
    parts.push(
      `<line x1="${padL}" y1="${py.toFixed(1)}" x2="${padL + plotW}" y2="${py.toFixed(1)}" stroke="#e5e7eb" stroke-width="1"/>`,
      `<text x="${padL - 8}" y="${(py + 4).toFixed(1)}" text-anchor="end" font-family="system-ui, sans-serif" font-size="11" fill="#374151">${fmtTick(yt)}</text>`,
    );
  }

  parts.push(
    `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#374151" stroke-width="1.2"/>`,
  );

  if (opts.hLine) {
    const hy = mapY(opts.hLine.y);
    parts.push(
      `<line x1="${padL}" y1="${hy.toFixed(1)}" x2="${padL + plotW}" y2="${hy.toFixed(1)}" stroke="${opts.hLine.color}" stroke-width="1.5" stroke-dasharray="6,4"/>`,
      `<text x="${padL + 8}" y="${(hy - 6).toFixed(1)}" font-family="system-ui, sans-serif" font-size="11" font-weight="600" fill="${opts.hLine.color}">${escapeXml(opts.hLine.label)}</text>`,
    );
  }

  for (const s of opts.series) {
    const coords = s.x
      .map((xv, i) => `${mapX(xv).toFixed(2)},${mapY(s.y[i] ?? 0.0).toFixed(2)}`)
      .join(" ");
    const dash = s.dashed ? ` stroke-dasharray="5,4"` : "";
    parts.push(
      `<polyline fill="none" stroke="${s.color}" stroke-width="2.2"${dash} points="${coords}"/>`,
    );
    if (s.errorBars) {
      for (let i = 0; i < s.x.length; i += 1) {
        const px = mapX(s.x[i] ?? 0.0);
        const yv = s.y[i] ?? 0.0;
        const eb = s.errorBars[i] ?? 0.0;
        const pyTop = mapY(yv + eb);
        const pyBot = mapY(yv - eb);
        parts.push(
          `<line x1="${px.toFixed(1)}" y1="${pyTop.toFixed(1)}" x2="${px.toFixed(1)}" y2="${pyBot.toFixed(1)}" stroke="${s.color}" stroke-width="1.4"/>`,
          `<circle cx="${px.toFixed(1)}" cy="${mapY(yv).toFixed(1)}" r="3.5" fill="${s.color}"/>`,
        );
      }
    }
  }

  // Legend top-right
  let ly = padT + 18;
  const lx = padL + plotW - 220;
  parts.push(
    `<rect x="${lx - 8}" y="${padT + 6}" width="218" height="${opts.series.length * 18 + 10}" fill="#ffffff" fill-opacity="0.9" stroke="#d1d5db" rx="4"/>`,
  );
  for (const s of opts.series) {
    const dash = s.dashed ? ` stroke-dasharray="4,3"` : "";
    parts.push(
      `<line x1="${lx}" y1="${ly}" x2="${lx + 22}" y2="${ly}" stroke="${s.color}" stroke-width="2.2"${dash}/>`,
      `<text x="${lx + 28}" y="${ly + 4}" font-family="system-ui, sans-serif" font-size="11" fill="#111827">${escapeXml(s.label)}</text>`,
    );
    ly += 18;
  }

  parts.push(
    `<text x="${padL + plotW / 2}" y="${height - 16}" text-anchor="middle" font-family="system-ui, sans-serif" font-size="12" fill="#111827">${escapeXml(opts.xLabel)}</text>`,
    `<text x="18" y="${padT + plotH / 2}" text-anchor="middle" transform="rotate(-90 18 ${padT + plotH / 2})" font-family="system-ui, sans-serif" font-size="12" fill="#111827">${escapeXml(opts.yLabel)}</text>`,
    `</svg>`,
  );
  return parts.join("\n");
}

/** Grouped bar chart (supports negative and positive values around a zero baseline). */
export function renderGroupedBarSvg(opts: {
  title: string;
  subtitle?: string;
  yLabel: string;
  groups: BarGroup[];
  hLine?: { y: number; label: string; color: string };
  width?: number;
  height?: number;
}): string {
  const width = opts.width ?? 820;
  const height = opts.height ?? 500;
  const padL = 78;
  const padR = 32;
  const padT = 76;
  const padB = 64;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  let yMin = 0.0;
  let yMax = opts.hLine ? opts.hLine.y : 0.0;
  for (const g of opts.groups) {
    for (const v of g.values) {
      if (v.value < yMin) yMin = v.value;
      if (v.value > yMax) yMax = v.value;
    }
  }
  const pad = 0.12 * Math.max(yMax - yMin, 1e-4);
  if (yMin < 0.0) yMin -= pad;
  yMax += pad;

  const mapY = (y: number): number => padT + plotH - ((y - yMin) / (yMax - yMin)) * plotH;
  const zeroY = mapY(0.0);
  const yTicks = niceTicks(yMin, yMax, 6);

  const parts: string[] = [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
    `<rect width="100%" height="100%" fill="#ffffff"/>`,
    `<text x="${width / 2}" y="24" text-anchor="middle" font-family="system-ui, sans-serif" font-size="16" font-weight="700" fill="#111827">${escapeXml(opts.title)}</text>`,
  ];
  if (opts.subtitle) {
    parts.push(
      `<text x="${width / 2}" y="42" text-anchor="middle" font-family="system-ui, sans-serif" font-size="12" fill="#4b5563">${escapeXml(opts.subtitle)}</text>`,
    );
  }

  for (const yt of yTicks) {
    const py = mapY(yt);
    parts.push(
      `<line x1="${padL}" y1="${py.toFixed(1)}" x2="${padL + plotW}" y2="${py.toFixed(1)}" stroke="#e5e7eb" stroke-width="1"/>`,
      `<text x="${padL - 8}" y="${(py + 4).toFixed(1)}" text-anchor="end" font-family="system-ui, sans-serif" font-size="11" fill="#374151">${fmtTick(yt)}</text>`,
    );
  }
  parts.push(
    `<line x1="${padL}" y1="${zeroY.toFixed(1)}" x2="${padL + plotW}" y2="${zeroY.toFixed(1)}" stroke="#374151" stroke-width="1.2"/>`,
    `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#374151" stroke-width="1.2"/>`,
  );

  if (opts.hLine) {
    const hy = mapY(opts.hLine.y);
    parts.push(
      `<line x1="${padL}" y1="${hy.toFixed(1)}" x2="${padL + plotW}" y2="${hy.toFixed(1)}" stroke="${opts.hLine.color}" stroke-width="1.5" stroke-dasharray="6,4"/>`,
      `<text x="${padL + 6}" y="${(hy - 5).toFixed(1)}" font-family="system-ui, sans-serif" font-size="11" font-weight="600" fill="${opts.hLine.color}">${escapeXml(opts.hLine.label)}</text>`,
    );
  }

  const nGroups = Math.max(opts.groups.length, 1);
  const groupW = plotW / nGroups;
  const legendSeen = new Map<string, string>();

  opts.groups.forEach((g, gi) => {
    const nBars = Math.max(g.values.length, 1);
    const innerW = 0.76 * groupW;
    const barW = innerW / nBars;
    const xStart = padL + gi * groupW + 0.12 * groupW;
    g.values.forEach((v, bi) => {
      legendSeen.set(v.label, v.color);
      const bx = xStart + bi * barW;
      const by = v.value >= 0.0 ? mapY(v.value) : zeroY;
      const bh = Math.max(1.0, Math.abs(mapY(v.value) - zeroY));
      parts.push(
        `<rect x="${(bx + 1).toFixed(1)}" y="${by.toFixed(1)}" width="${Math.max(2, barW - 2).toFixed(1)}" height="${bh.toFixed(1)}" fill="${v.color}" fill-opacity="0.88"/>`,
      );
    });
    parts.push(
      `<text x="${(padL + (gi + 0.5) * groupW).toFixed(1)}" y="${padT + plotH + 20}" text-anchor="middle" font-family="system-ui, sans-serif" font-size="11" fill="#111827">${escapeXml(g.category)}</text>`,
    );
  });

  // Top legend bar
  let lx = padL;
  for (const [label, color] of legendSeen.entries()) {
    parts.push(
      `<rect x="${lx}" y="53" width="12" height="12" fill="${color}"/>`,
      `<text x="${lx + 16}" y="63" font-family="system-ui, sans-serif" font-size="11" fill="#111827">${escapeXml(label)}</text>`,
    );
    lx += 16 + label.length * 7.2 + 18;
  }

  parts.push(
    `<text x="18" y="${padT + plotH / 2}" text-anchor="middle" transform="rotate(-90 18 ${padT + plotH / 2})" font-family="system-ui, sans-serif" font-size="12" fill="#111827">${escapeXml(opts.yLabel)}</text>`,
    `</svg>`,
  );
  return parts.join("\n");
}

/** Histogram comparison figure for two parameter distributions (e.g. GARCH t-MLE vs Gauss QMLE). */
export function renderHistogramSvg(opts: {
  title: string;
  subtitle?: string;
  xLabel: string;
  edges: number[];
  series: Array<{ label: string; density: number[]; color: string }>;
  trueValue?: number;
  width?: number;
  height?: number;
}): string {
  const width = opts.width ?? 820;
  const height = opts.height ?? 480;
  const padL = 74;
  const padR = 30;
  const padT = 62;
  const padB = 60;
  const plotW = width - padL - padR;
  const plotH = height - padT - padB;

  const xMin = opts.edges[0] ?? 0.0;
  const xMax = opts.edges[opts.edges.length - 1] ?? 1.0;
  let yMax = 0.0;
  for (const s of opts.series) {
    for (const d of s.density) {
      if (d > yMax) yMax = d;
    }
  }
  yMax *= 1.12;
  const mapX = (x: number): number => padL + ((x - xMin) / (xMax - xMin)) * plotW;
  const mapY = (y: number): number => padT + plotH - (y / Math.max(yMax, 1e-12)) * plotH;

  const parts: string[] = [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
    `<rect width="100%" height="100%" fill="#ffffff"/>`,
    `<text x="${width / 2}" y="25" text-anchor="middle" font-family="system-ui, sans-serif" font-size="16" font-weight="700" fill="#111827">${escapeXml(opts.title)}</text>`,
  ];
  if (opts.subtitle) {
    parts.push(
      `<text x="${width / 2}" y="44" text-anchor="middle" font-family="system-ui, sans-serif" font-size="12" fill="#4b5563">${escapeXml(opts.subtitle)}</text>`,
    );
  }

  for (const xt of niceTicks(xMin, xMax, 6)) {
    const px = mapX(xt);
    parts.push(
      `<line x1="${px.toFixed(1)}" y1="${padT}" x2="${px.toFixed(1)}" y2="${padT + plotH}" stroke="#e5e7eb"/>`,
      `<text x="${px.toFixed(1)}" y="${padT + plotH + 18}" text-anchor="middle" font-family="system-ui, sans-serif" font-size="11" fill="#374151">${fmtTick(xt)}</text>`,
    );
  }
  for (const yt of niceTicks(0, yMax, 5)) {
    const py = mapY(yt);
    parts.push(
      `<line x1="${padL}" y1="${py.toFixed(1)}" x2="${padL + plotW}" y2="${py.toFixed(1)}" stroke="#e5e7eb"/>`,
      `<text x="${padL - 8}" y="${(py + 4).toFixed(1)}" text-anchor="end" font-family="system-ui, sans-serif" font-size="11" fill="#374151">${fmtTick(yt)}</text>`,
    );
  }

  for (const s of opts.series) {
    for (let i = 0; i < s.density.length; i += 1) {
      const x0 = mapX(opts.edges[i] ?? xMin);
      const x1 = mapX(opts.edges[i + 1] ?? xMax);
      const yTop = mapY(s.density[i] ?? 0.0);
      parts.push(
        `<rect x="${x0.toFixed(1)}" y="${yTop.toFixed(1)}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="${(padT + plotH - yTop).toFixed(1)}" fill="${s.color}" fill-opacity="0.38" stroke="${s.color}" stroke-width="1"/>`,
      );
    }
  }

  if (opts.trueValue !== undefined) {
    const tx = mapX(opts.trueValue);
    parts.push(
      `<line x1="${tx.toFixed(1)}" y1="${padT}" x2="${tx.toFixed(1)}" y2="${padT + plotH}" stroke="#111827" stroke-width="2" stroke-dasharray="6,4"/>`,
      `<text x="${(tx + 5).toFixed(1)}" y="${padT + 16}" font-family="system-ui, sans-serif" font-size="11" font-weight="700" fill="#111827">Benar=${fmtTick(opts.trueValue)}</text>`,
    );
  }

  parts.push(
    `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#374151" stroke-width="1.2"/>`,
    `<text x="${padL + plotW / 2}" y="${height - 14}" text-anchor="middle" font-family="system-ui, sans-serif" font-size="12" fill="#111827">${escapeXml(opts.xLabel)}</text>`,
    `<text x="18" y="${padT + plotH / 2}" text-anchor="middle" transform="rotate(-90 18 ${padT + plotH / 2})" font-family="system-ui, sans-serif" font-size="12" fill="#111827">Kepadatan</text>`,
  );

  let ly = padT + 18;
  const lx = padL + plotW - 190;
  for (const s of opts.series) {
    parts.push(
      `<rect x="${lx}" y="${ly - 9}" width="12" height="12" fill="${s.color}" fill-opacity="0.5" stroke="${s.color}"/>`,
      `<text x="${lx + 18}" y="${ly + 1}" font-family="system-ui, sans-serif" font-size="11" fill="#111827">${escapeXml(s.label)}</text>`,
    );
    ly += 18;
  }

  parts.push(`</svg>`);
  return parts.join("\n");
}
