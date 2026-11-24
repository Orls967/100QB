#!/usr/bin/env bash
set -euo pipefail
PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$PROJ/../../.." && pwd)"
PY="$ROOT/.venv/bin/python"
[[ -x "$PY" ]] || PY="python3"

cargo build --release --manifest-path "$PROJ/rust/Cargo.toml"
(cd "$PROJ/go" && go build -o "$PROJ/rust/target/release/orchestrate" .)
"$PROJ/rust/target/release/orchestrate" \
  --bin "$PROJ/rust/target/release/percolate" \
  --seed 20261124 --timeout 55 \
  --out "$PROJ/rust/target/go_perc.json"
"$PROJ/rust/target/release/percolate" sample 96 0.5 20261124 > "$PROJ/rust/target/sample_96.json"
"$PY" "$PROJ/py/percolation.py" \
  --go-json "$PROJ/rust/target/go_perc.json" \
  --sample-json "$PROJ/rust/target/sample_96.json" \
  --out-json "$ROOT/outputs/day006c/results.json" \
  --fig-dir "$ROOT/outputs/figures"
