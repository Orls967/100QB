CC       := cc
CXX      := c++
CFLAGS   := -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off

BUILD    := build

.PHONY: all build-c build-cpp build-rust build-go test run clean

all: build-c build-cpp build-rust build-go

$(BUILD):
	mkdir -p $(BUILD)

build-c: c/walk.c | $(BUILD)
	$(CC) $(CFLAGS) c/walk.c -o $(BUILD)/walk_c -lm

build-cpp: cpp/walk.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) cpp/walk.cpp -o $(BUILD)/walk_cpp

build-rust: rust/src/main.rs rust/Cargo.toml
	cd rust && cargo build --release --offline

build-go: go/main.go go/go.mod | $(BUILD)
	cd go && go vet ./... && go build -o ../$(BUILD)/walk_go .

test: all
	./$(BUILD)/walk_c s01
	./$(BUILD)/walk_cpp s01
	./rust/target/release/walk_rust s01
	./$(BUILD)/walk_go s01
	ruby -wc ruby/walk.rb && ruby ruby/walk.rb s01
	../../../.venv/bin/python py/walk.py s01
	node --experimental-strip-types ts/walk.ts s01

run: all
	bash run_all.sh

clean:
	rm -rf $(BUILD) rust/target
