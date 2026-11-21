/**
 * Strict deterministic 17-significant-digit JSON serializer matching project conventions.
 * Rejects NaN and +-Infinity to catch numerical failures immediately.
 */

function sanitizeValue(val: unknown, path: string): unknown {
  if (typeof val === "number") {
    if (!Number.isFinite(val)) {
      throw new Error(`Non-finite float at ${path}: ${val}`);
    }
    // Preserve full IEEE-754 precision via shortest round-trip or 17 sig figs
    return Number(val.toPrecision(17));
  }
  if (Array.isArray(val)) {
    return val.map((item, idx) => sanitizeValue(item, `${path}[${idx}]`));
  }
  if (val !== null && typeof val === "object") {
    const out: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(val)) {
      if (v !== undefined) {
        out[k] = sanitizeValue(v, `${path}.${k}`);
      }
    }
    return out;
  }
  return val;
}

export function serializeJson(data: unknown): string {
  return JSON.stringify(sanitizeValue(data, "$"), null, 2);
}
