#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
make all

echo "=== [1/8] C11 ==="
./build/sc_c all
echo "=== [2/8] C++20 ==="
./build/sc_cpp all
echo "=== [3/8] Rust ==="
./rust/target/release/sc_rust all
echo "=== [4/8] Go ==="
./build/sc_go all
echo "=== [5/8] Ruby ==="
ruby ruby/main.rb all
echo "=== [6/8] Python ==="
../../../.venv/bin/python py/main.py all
echo "=== [7/8] TypeScript ==="
node --experimental-strip-types ts/main.ts all
echo "=== [8/8] Jupyter Notebook ==="
../../../.venv/bin/python -c "
import nbformat
from nbclient import NotebookClient
p = '../../../days/day008f_sim_supply-chain-octalingual.ipynb'
nb = nbformat.read(p, as_version=4)
NotebookClient(nb, timeout=300, kernel_name='python3').execute(cwd='../../../days')
for i, c in enumerate(nb['cells']):
    c['id'] = f'cell-{i:02d}'
    c.get('metadata', {}).pop('execution', None)
nbformat.write(nb, p)
"
echo "=== Compare & Figures ==="
../../../.venv/bin/python compare.py
