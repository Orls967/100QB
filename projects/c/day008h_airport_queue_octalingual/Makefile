CC ?= cc
CXX ?= c++
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off

C_SRCS = c/rng.c c/config.c c/dtmc.c c/station.c c/passenger.c c/network.c c/engine.c c/stats.c c/analysis.c c/svg.c c/main.c
CPP_SRCS = cpp/rng.cpp cpp/config.cpp cpp/dtmc.cpp cpp/station.cpp cpp/passenger.cpp cpp/network.cpp cpp/engine.cpp cpp/stats.cpp cpp/analysis.cpp cpp/svg.cpp cpp/main.cpp

.PHONY: all build build-c build-cpp build-rust build-go test clean

build: build-c build-cpp build-rust build-go

build-c:
	@mkdir -p build
	$(CC) $(CFLAGS) $(C_SRCS) -o build/aq_c -lm

build-cpp:
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(CPP_SRCS) -o build/aq_cpp

build-rust:
	cd rust && cargo build --release --offline

build-go:
	@mkdir -p build
	cd go && go vet ./... && go build -o ../build/aq_go .

all: build
	./run_all.sh

test: build
	./build/aq_c s01 --config config.kv --out out/c
	./build/aq_cpp s01 --config config.kv --out out/cpp
	./rust/target/release/aq_rust s01
	./build/aq_go s01
	ruby -wc ruby/*.rb && ruby ruby/main.rb s01
	../../../.venv/bin/python py/main.py s01
	node --experimental-strip-types ts/main.ts s01
	for f in s01_selftest.txt; do for lg in cpp rust go ruby py ts; do cmp -s out/c/$$f out/$$lg/$$f || exit 1; done; done
	@echo "day008h selftest 7-CLI parity OK"

clean:
	rm -rf build rust/target
