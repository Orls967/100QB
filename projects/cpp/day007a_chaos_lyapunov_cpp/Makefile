CXX ?= c++
CXXFLAGS = -O3 -std=c++20 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
OUT_DIR = ../../../outputs/day007a
FIG_DIR = ../../../outputs/figures

.PHONY: all test run clean

all: build/chaos build/test_chaos

build/chaos: src/chaos.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/chaos.cpp -o $@

build/test_chaos: tests/test_chaos.cpp src/chaos.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) tests/test_chaos.cpp -o $@

test: build/test_chaos
	./build/test_chaos

run: build/chaos
	@mkdir -p $(OUT_DIR) $(FIG_DIR)
	./build/chaos $(OUT_DIR) $(FIG_DIR)

clean:
	rm -rf build
