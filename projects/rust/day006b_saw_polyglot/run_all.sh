#!/usr/bin/env bash
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$PROJ/../../.." && pwd)"
PY="$ROOT/.venv/bin/python"
[[ -x "$PY" ]] || PY="python3"
T0_FILE="${T0_FILE:-${TMPDIR:-/tmp}/day006_t0}"
[[ -f "$T0_FILE" ]] || date +%s > "$T0_FILE"

cargo build --release --manifest-path "$PROJ/rust/Cargo.toml"
(cd "$PROJ/go" && go build -o "$PROJ/rust/target/release/orchestrate" .)
timeout 60 "$PROJ/rust/target/release/orchestrate" \
  --saw-bin "$PROJ/rust/target/release/saw" \
  --nmax 18 --k 7 --heatmap-n 14 --seed 20261124 \
  --out "$PROJ/rust/target/go_out.json"
"$PY" "$PROJ/py/saw_analysis.py" \
  --go-json "$PROJ/rust/target/go_out.json" \
  --out-json "$ROOT/outputs/day006b/results.json" \
  --fig-dir "$ROOT/outputs/figures" \
  --t0-file "$T0_FILE"
