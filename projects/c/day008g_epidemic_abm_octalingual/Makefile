CC ?= cc
CXX ?= c++
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off

C_SRCS = c/rng.c c/population.c c/contact.c c/disease.c c/interventions.c c/engine.c c/analysis.c c/svg.c c/main.c
CPP_SRCS = cpp/rng.cpp cpp/population.cpp cpp/contact.cpp cpp/disease.cpp cpp/interventions.cpp cpp/engine.cpp cpp/analysis.cpp cpp/svg.cpp cpp/main.cpp

.PHONY: all build build-c build-cpp build-rust build-go test clean

build: build-c build-cpp build-rust build-go

build-c:
	@mkdir -p build
	$(CC) $(CFLAGS) $(C_SRCS) -o build/epi_c -lm

build-cpp:
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(CPP_SRCS) -o build/epi_cpp

build-rust:
	cd rust && cargo build --release --offline

build-go:
	@mkdir -p build
	cd go && go vet ./... && go build -o ../build/epi_go .

all: build
	./run_all.sh

test: build
	./build/epi_c s01 --config config.kv --out out/c
	./build/epi_cpp s01 --config config.kv --out out/cpp
	./rust/target/release/epi_rust s01
	./build/epi_go s01
	ruby -wc ruby/*.rb && ruby ruby/main.rb s01
	../../../.venv/bin/python py/main.py s01
	node --experimental-strip-types ts/main.ts s01
	for f in s01_selftest.txt; do for lg in cpp rust go ruby py ts; do cmp -s out/c/$$f out/$$lg/$$f || exit 1; done; done
	@echo "day008g selftest 7-CLI parity OK"

clean:
	rm -rf build rust/target
