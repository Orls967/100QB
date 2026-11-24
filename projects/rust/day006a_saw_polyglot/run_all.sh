#!/usr/bin/env bash
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$PROJ_DIR/../../.." && pwd)"

PY_BIN="python3"
if [[ -x "$REPO_ROOT/.venv/bin/python" ]]; then
  PY_BIN="$REPO_ROOT/.venv/bin/python"
fi

NMAX="${NMAX:-24}"
PREFIX_K="${PREFIX_K:-8}"
SEED="${SEED:-20261124}"
T0_EPOCH="${T0_EPOCH:-$(date +%s)}"
T0_STR="${T0_STR:-$(TZ=Asia/Makassar date)}"

echo "[1/3] Building Rust core (release)..."
cargo build --release --manifest-path "$PROJ_DIR/rust/Cargo.toml"
SAW_BIN="$PROJ_DIR/rust/target/release/saw"

echo "[2/3] Building and running Go concurrent orchestrator (Nmax=$NMAX, K=$PREFIX_K)..."
ORCH_BIN="$PROJ_DIR/rust/target/release/orchestrate"
(cd "$PROJ_DIR/go" && go build -o "$ORCH_BIN" .)
"$ORCH_BIN" \
  --saw-bin "$SAW_BIN" \
  --nmax "$NMAX" \
  --k "$PREFIX_K" \
  --heatmap-n 16 \
  --seed "$SEED" \
  --timeout-sec 90 \
  --out "$PROJ_DIR/rust/target/go_orchestrated.json"

echo "[3/3] Running Python reference, critical exponent regressions, and visualizations..."
"$PY_BIN" "$PROJ_DIR/py/saw_analysis.py" \
  --go-json "$PROJ_DIR/rust/target/go_orchestrated.json" \
  --out-json "$REPO_ROOT/outputs/day006a/results.json" \
  --fig-dir "$REPO_ROOT/outputs/figures" \
  --t0-epoch "$T0_EPOCH" \
  --t0-str "$T0_STR"

echo "Done -> $REPO_ROOT/outputs/day006a/results.json"
