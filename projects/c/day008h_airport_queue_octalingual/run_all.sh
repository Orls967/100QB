#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

make build
./build/aq_c all --config config.kv --out out/c
./build/aq_cpp all --config config.kv --out out/cpp
./rust/target/release/aq_rust all
./build/aq_go all
ruby -wc ruby/*.rb
ruby ruby/main.rb all
../../../.venv/bin/python py/main.py all
node --experimental-strip-types ts/main.ts all
