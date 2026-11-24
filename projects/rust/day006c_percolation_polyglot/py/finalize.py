#!/usr/bin/env python3
"""Single wall-clock finalizer for day006c results.json."""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[4]


def main() -> int:
    t0_path = Path(os.environ.get("TMPDIR", "/tmp")) / "day006_t0"
    t0 = int(t0_path.read_text(encoding="utf-8").strip())
    now = int(time.time())
    elapsed = now - t0

    res_path = REPO_ROOT / "outputs" / "day006c" / "results.json"
    doc = json.loads(res_path.read_text(encoding="utf-8"))
    doc["timing"] = {"t0_epoch": t0, "finalize_epoch": now, "elapsed_sec": elapsed}
    doc["criteria"]["c_wall_clock"] = {
        "status": "tercapai" if elapsed <= 570 else "belum tercapai",
        "measured": elapsed,
        "threshold": "elapsed_sec <= 570",
    }
    serialized = json.dumps(doc, indent=2) + "\n"
    for forbidden in ("/Users/", "/home/", "/tmp", "file://"):
        if forbidden in serialized:
            print(f"Forbidden path substring found: {forbidden}", file=sys.stderr)
            return 1
    res_path.write_text(serialized, encoding="utf-8")
    print(f"finalize.py OK: elapsed_sec={elapsed}, status={doc['criteria']['c_wall_clock']['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
