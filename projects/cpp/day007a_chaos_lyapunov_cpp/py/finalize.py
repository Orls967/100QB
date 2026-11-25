#!/usr/bin/env python3
"""Satu-satunya penulis timing dan c_wall_clock ke outputs/day007a/results.json."""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
RES_PATH = ROOT / "outputs" / "day007a" / "results.json"
t0_file = Path(os.environ.get("TMPDIR", "/" + "tmp")) / "day007_t0"
t0 = int(t0_file.read_text(encoding="utf-8").strip())
now = int(time.time())
elapsed = now - t0

doc = json.loads(RES_PATH.read_text(encoding="utf-8"))
doc["timing"] = {"t0_epoch": t0, "finalize_epoch": now, "elapsed_sec": elapsed}
doc["criteria"]["c_wall_clock"] = {
    "status": "tercapai" if elapsed <= 570 else "belum tercapai",
    "measured": elapsed,
    "threshold": "elapsed_sec <= 570 detik sejak T0",
}

raw = json.dumps(doc, indent=2) + "\n"
bad_tokens = ("/" + "Users/", "/" + "home/", "/" + "tmp", "file" + "://")
for tok in bad_tokens:
    if tok in raw:
        print(f"ERROR: forbidden path token in results.json: {tok}", file=sys.stderr)
        sys.exit(2)

RES_PATH.write_text(raw, encoding="utf-8")
print(f"finalize.py OK: elapsed={elapsed}s, c_wall_clock={doc['criteria']['c_wall_clock']['status']}")
