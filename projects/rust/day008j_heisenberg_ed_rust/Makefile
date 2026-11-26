CARGO ?= cargo
LINES = $(shell wc -l Cargo.toml $(sort $(wildcard src/*.rs)) | tail -n 1 | awk '{print $$1}')
RUSTC_VER = $(shell rustc --version)
ARCH = $(shell uname -m)
CORES = $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

.PHONY: all test run lines clean

all:
	$(CARGO) build --release --offline

test: all
	./target/release/heis selftest

run: all
	./target/release/heis all --lines $(LINES) --rustc "$(RUSTC_VER)" --arch "$(ARCH)" --cores $(CORES)

lines:
	@wc -l Cargo.toml $(sort $(wildcard src/*.rs))
	@test $(LINES) -le 499

clean:
	$(CARGO) clean
