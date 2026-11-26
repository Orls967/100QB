#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

echo "=== Build ==="
make all

echo "=== [1/8] C ==="
./build/walk_c all

echo "=== [2/8] C++ ==="
./build/walk_cpp all

echo "=== [3/8] Rust ==="
./rust/target/release/walk_rust all

echo "=== [4/8] Go ==="
./build/walk_go all

echo "=== [5/8] Ruby ==="
ruby ruby/walk.rb all

echo "=== [6/8] Python ==="
../../../.venv/bin/python py/walk.py all

echo "=== [7/8] TypeScript (node --experimental-strip-types) ==="
node --experimental-strip-types ts/walk.ts all

echo "=== [8/8] Notebook (.ipynb) ==="
../../../.venv/bin/python - << 'EOF'
import nbformat
from nbclient import NotebookClient
p = "../../../days/day008d_sim_lattice-walks-octalingual.ipynb"
nb = nbformat.read(p, as_version=4)
NotebookClient(nb, timeout=600, kernel_name="python3").execute()
nbformat.write(nb, p)
EOF

echo "=== Compare & Aggregate ==="
../../../.venv/bin/python compare.py
