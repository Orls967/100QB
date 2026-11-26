CC       := cc
CXX      := c++
CFLAGS   := -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
BUILD    := build

.PHONY: all build-c build-cpp build-rust build-go test run clean

all: build-c build-cpp build-rust build-go

$(BUILD):
	mkdir -p $(BUILD)

build-c: c/poisson.c | $(BUILD)
	$(CC) $(CFLAGS) c/poisson.c -o $(BUILD)/poisson_c -lm

build-cpp: cpp/poisson.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) cpp/poisson.cpp -o $(BUILD)/poisson_cpp

build-rust: rust/src/main.rs rust/Cargo.toml
	cd rust && cargo build --release --offline

build-go: go/main.go go/go.mod | $(BUILD)
	cd go && go vet ./... && go build -o ../$(BUILD)/poisson_go .

test: all
	./$(BUILD)/poisson_c s01
	./$(BUILD)/poisson_cpp s01
	./rust/target/release/poisson_rust s01
	./$(BUILD)/poisson_go s01
	ruby -wc ruby/poisson.rb && ruby ruby/poisson.rb s01
	../../../.venv/bin/python py/poisson.py s01
	node --experimental-strip-types ts/poisson.ts s01

run: all
	bash run_all.sh

clean:
	rm -rf $(BUILD) rust/target
