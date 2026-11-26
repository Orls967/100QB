CC       := cc
CXX      := c++
CFLAGS   := -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
BUILD    := build

C_SRCS   := c/rng.c c/inventory.c c/allocation.c c/entities.c c/engine.c c/clearing.c c/analysis.c c/svg.c c/main.c
CPP_SRCS := cpp/rng.cpp cpp/inventory.cpp cpp/allocation.cpp cpp/entities.cpp cpp/engine.cpp cpp/clearing.cpp cpp/analysis.cpp cpp/svg.cpp cpp/main.cpp

.PHONY: all build-c build-cpp build-rust build-go test run clean

all: build-c build-cpp build-rust build-go

$(BUILD):
	mkdir -p $(BUILD)

build-c: $(C_SRCS) c/common.h | $(BUILD)
	$(CC) $(CFLAGS) $(C_SRCS) -o $(BUILD)/sc_c -lm

build-cpp: $(CPP_SRCS) cpp/common.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(CPP_SRCS) -o $(BUILD)/sc_cpp

build-rust: rust/src/main.rs rust/Cargo.toml
	cd rust && cargo build --release --offline

build-go: go/main.go go/go.mod | $(BUILD)
	cd go && go vet ./... && go build -o ../$(BUILD)/sc_go .

test: all
	./$(BUILD)/sc_c s01
	./$(BUILD)/sc_cpp s01
	./rust/target/release/sc_rust s01
	./$(BUILD)/sc_go s01
	ruby -wc ruby/*.rb && ruby ruby/main.rb s01
	../../../.venv/bin/python py/main.py s01
	node --experimental-strip-types ts/main.ts s01

run: all
	bash run_all.sh

clean:
	rm -rf $(BUILD) rust/target
