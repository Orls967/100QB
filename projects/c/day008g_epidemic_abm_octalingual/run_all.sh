#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

make build
./build/epi_c all --config config.kv --out out/c
./build/epi_cpp all --config config.kv --out out/cpp
./rust/target/release/epi_rust all
./build/epi_go all
ruby -wc ruby/*.rb
ruby ruby/main.rb all
../../../.venv/bin/python py/main.py all
node --experimental-strip-types ts/main.ts all
