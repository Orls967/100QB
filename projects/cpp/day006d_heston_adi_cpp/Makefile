CXX ?= c++
BASE_FLAGS = -O3 -std=c++20 -fno-fast-math -ffp-contract=off -Wall -Wextra -Wpedantic -pthread
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),x86_64)
  NATIVE_FLAG = -march=native
else
  NATIVE_FLAG := $(shell echo 'int main(){}' | $(CXX) -x c++ -mcpu=native - -o /dev/null 2>/dev/null && echo -mcpu=native || (echo 'int main(){}' | $(CXX) -x c++ -mcpu=apple-m4 - -o /dev/null 2>/dev/null && echo -mcpu=apple-m4 || echo -mcpu=apple-m3))
endif

VARIANTS = $(wildcard src/v*.cpp)
TESTS = $(wildcard tests/*.cpp)
OUT_DIR = ../../../outputs/day006d

.PHONY: all test bench vecreport clean

all: build/heston_tests $(if $(wildcard src/bench.cpp),build/heston_bench build/heston_bench_native)

build/heston_tests: $(TESTS) $(VARIANTS) src/common.hpp src/reference.hpp src/variants.hpp tests/check.hpp
	@mkdir -p build
	$(CXX) $(BASE_FLAGS) -Isrc -Itests $(TESTS) $(VARIANTS) -o $@

build/heston_bench: src/bench.cpp $(VARIANTS) src/common.hpp src/reference.hpp src/variants.hpp
	@mkdir -p build
	$(CXX) $(BASE_FLAGS) -Isrc src/bench.cpp $(VARIANTS) -o $@

build/heston_bench_native: src/bench.cpp $(VARIANTS) src/common.hpp src/reference.hpp src/variants.hpp
	@mkdir -p build
	$(CXX) $(BASE_FLAGS) $(NATIVE_FLAG) -Isrc src/bench.cpp $(VARIANTS) -o $@

test: build/heston_tests
	./build/heston_tests

bench: build/heston_bench build/heston_bench_native
	@mkdir -p $(OUT_DIR)
	./build/heston_bench all $(OUT_DIR)
	./build/heston_bench_native native $(OUT_DIR)

vecreport:
	@mkdir -p build $(OUT_DIR)
	$(CXX) $(BASE_FLAGS) -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Isrc -c src/v3_batched.cpp -o build/v3.o 2> $(OUT_DIR)/vectorize_remarks.txt
	$(CXX) $(BASE_FLAGS) -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Isrc -c src/v4_fused.cpp -o build/v4.o 2>> $(OUT_DIR)/vectorize_remarks.txt

clean:
	rm -rf build
