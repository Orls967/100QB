/** Pure TypeScript SVG chart generator for convergence and estimator comparison. */

export interface ConvergencePoint {
  n: number;
  se_naive: number;
  se_antithetic: number;
  se_control_variate: number;
}

export interface EstimatorSummary {
  name: string;
  price: number;
  se: number;
  ci_low: number;
  ci_high: number;
  vrf: number;
}

export interface Slopes {
  naive: number;
  antithetic: number;
  control_variate: number;
}

function escapeXml(text: string): string {
  return text
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

export function renderConvergenceSvg(points: readonly ConvergencePoint[], slopes: Slopes): string {
  const w = 900;
  const h = 520;
  const padL = 85;
  const padR = 30;
  const padT = 55;
  const padB = 65;
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;

  const minLogX = Math.log10(900);
  const maxLogX = Math.log10(230_000);
  let minSe = Number.POSITIVE_INFINITY;
  let maxSe = 0.0;
  for (const pt of points) {
    minSe = Math.min(minSe, pt.se_control_variate);
    maxSe = Math.max(maxSe, pt.se_naive);
  }
  const minLogY = Math.log10(minSe * 0.7);
  const maxLogY = Math.log10(maxSe * 1.35);

  const mapX = (n: number): number =>
    padL + ((Math.log10(n) - minLogX) / (maxLogX - minLogX)) * plotW;
  const mapY = (se: number): number =>
    padT + (1.0 - (Math.log10(se) - minLogY) / (maxLogY - minLogY)) * plotH;

  const gridLines: string[] = [];
  for (const nTick of [1_000, 10_000, 100_000]) {
    const x = mapX(nTick);
    gridLines.push(
      `<line x1="${x.toFixed(1)}" y1="${padT}" x2="${x.toFixed(1)}" y2="${padT + plotH}" stroke="#e5e7eb" stroke-dasharray="4,4"/>`,
      `<text x="${x.toFixed(1)}" y="${padT + plotH + 22}" text-anchor="middle" font-size="12" fill="#374151">${nTick.toLocaleString("en-US")}</text>`,
    );
  }
  for (let exp = Math.ceil(minLogY); exp <= Math.floor(maxLogY); exp += 1) {
    const val = 10 ** exp;
    const y = mapY(val);
    gridLines.push(
      `<line x1="${padL}" y1="${y.toFixed(1)}" x2="${padL + plotW}" y2="${y.toFixed(1)}" stroke="#e5e7eb" stroke-dasharray="4,4"/>`,
      `<text x="${padL - 10}" y="${(y + 4).toFixed(1)}" text-anchor="end" font-size="12" fill="#374151">${val.toExponential(0)}</text>`,
    );
  }

  const series = [
    {
      label: `Naive MC (slope=${slopes.naive.toFixed(3)})`,
      color: "#1f6feb",
      getY: (pt: ConvergencePoint) => pt.se_naive,
    },
    {
      label: `Antithetic (slope=${slopes.antithetic.toFixed(3)})`,
      color: "#d1495b",
      getY: (pt: ConvergencePoint) => pt.se_antithetic,
    },
    {
      label: `Control Variate (slope=${slopes.control_variate.toFixed(3)})`,
      color: "#2a9d8f",
      getY: (pt: ConvergencePoint) => pt.se_control_variate,
    },
  ];

  const pathsSvg: string[] = [];
  const legendSvg: string[] = [];
  series.forEach((s, idx) => {
    const coords = points.map((pt) => `${mapX(pt.n).toFixed(1)},${mapY(s.getY(pt)).toFixed(1)}`);
    pathsSvg.push(
      `<polyline fill="none" stroke="${s.color}" stroke-width="2.5" points="${coords.join(" ")}"/>`,
    );
    for (const pt of points) {
      pathsSvg.push(
        `<circle cx="${mapX(pt.n).toFixed(1)}" cy="${mapY(s.getY(pt)).toFixed(1)}" r="4" fill="${s.color}"/>`,
      );
    }
    const ly = padT + 22 + idx * 22;
    const lx = padL + plotW - 240;
    legendSvg.push(
      `<line x1="${lx}" y1="${ly}" x2="${lx + 24}" y2="${ly}" stroke="${s.color}" stroke-width="2.5"/>`,
      `<text x="${lx + 32}" y="${ly + 4}" font-size="12" fill="#111827">${escapeXml(s.label)}</text>`,
    );
  });

  return [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="system-ui, -apple-system, sans-serif">`,
    `<rect width="${w}" height="${h}" fill="#ffffff"/>`,
    `<text x="${w / 2}" y="30" text-anchor="middle" font-size="18" font-weight="bold" fill="#111827">Konvergensi Standard Error vs Jumlah Path (Log-Log, TypeScript)</text>`,
    ...gridLines,
    `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#9ca3af"/>`,
    ...pathsSvg,
    `<rect x="${padL + plotW - 252}" y="${padT + 8}" width="244" height="72" rx="6" fill="#ffffff" fill-opacity="0.9" stroke="#d1d5db"/>`,
    ...legendSvg,
    `<text x="${padL + plotW / 2}" y="${h - 16}" text-anchor="middle" font-size="13" fill="#111827">Jumlah Path N (skala log)</text>`,
    `<text x="22" y="${padT + plotH / 2}" text-anchor="middle" font-size="13" fill="#111827" transform="rotate(-90 22 ${padT + plotH / 2})">Standard Error (SE)</text>`,
    `</svg>`,
  ].join("\n");
}

export function renderEstimatorsSvg(estimators: readonly EstimatorSummary[]): string {
  const w = 900;
  const h = 440;
  const padL = 210;
  const padR = 60;
  const padT = 60;
  const padB = 60;
  const plotW = w - padL - padR;
  const plotH = h - padT - padB;

  let minX = Number.POSITIVE_INFINITY;
  let maxX = Number.NEGATIVE_INFINITY;
  for (const e of estimators) {
    minX = Math.min(minX, e.ci_low);
    maxX = Math.max(maxX, e.ci_high);
  }
  const margin = (maxX - minX) * 0.15;
  minX -= margin;
  maxX += margin;

  const mapX = (val: number): number => padL + ((val - minX) / (maxX - minX)) * plotW;
  const colors = ["#1f6feb", "#d1495b", "#2a9d8f"];

  const rowsSvg: string[] = [];
  estimators.forEach((e, idx) => {
    const cy = padT + ((idx + 0.5) / estimators.length) * plotH;
    const xLo = mapX(e.ci_low);
    const xHi = mapX(e.ci_high);
    const xMid = mapX(e.price);
    const col = colors[idx % colors.length] ?? "#1f6feb";
    rowsSvg.push(
      `<text x="${padL - 12}" y="${cy + 4}" text-anchor="end" font-size="12" fill="#111827">${escapeXml(e.name)}</text>`,
      `<line x1="${xLo.toFixed(1)}" y1="${cy.toFixed(1)}" x2="${xHi.toFixed(1)}" y2="${cy.toFixed(1)}" stroke="${col}" stroke-width="3"/>`,
      `<line x1="${xLo.toFixed(1)}" y1="${(cy - 6).toFixed(1)}" x2="${xLo.toFixed(1)}" y2="${(cy + 6).toFixed(1)}" stroke="${col}" stroke-width="2"/>`,
      `<line x1="${xHi.toFixed(1)}" y1="${(cy - 6).toFixed(1)}" x2="${xHi.toFixed(1)}" y2="${(cy + 6).toFixed(1)}" stroke="${col}" stroke-width="2"/>`,
      `<circle cx="${xMid.toFixed(1)}" cy="${cy.toFixed(1)}" r="5" fill="${col}"/>`,
      `<text x="${(xHi + 8).toFixed(1)}" y="${cy + 4}" font-size="11" fill="#374151">${e.price.toFixed(4)} (VRF ${e.vrf.toFixed(1)}x)</text>`,
    );
  });

  return [
    `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="system-ui, -apple-system, sans-serif">`,
    `<rect width="${w}" height="${h}" fill="#ffffff"/>`,
    `<text x="${w / 2}" y="32" text-anchor="middle" font-size="18" font-weight="bold" fill="#111827">Perbandingan Estimasi Harga + CI 95% &amp; VRF (200k Path, TypeScript)</text>`,
    `<rect x="${padL}" y="${padT}" width="${plotW}" height="${plotH}" fill="none" stroke="#9ca3af"/>`,
    ...rowsSvg,
    `<text x="${padL + plotW / 2}" y="${h - 18}" text-anchor="middle" font-size="13" fill="#111827">Harga Opsi Asia</text>`,
    `</svg>`,
  ].join("\n");
}
