CC ?= cc
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
SRCS = $(sort $(wildcard src/*.c))
HDRS = $(sort $(wildcard src/*.h))
LINES = $(shell wc -l Makefile $(SRCS) $(HDRS) | tail -n 1 | awk '{print $$1}')
CC_VER = $(shell $(CC) --version | head -n 1)
ARCH = $(shell uname -m)
CORES = $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

.PHONY: all test run lines clean

all: build/sandpile

build/sandpile: $(SRCS) $(HDRS) Makefile
	@mkdir -p build out ../../../outputs/day008i ../../../outputs/figures
	$(CC) $(CFLAGS) $(SRCS) -o build/sandpile -lm

test: build/sandpile
	./build/sandpile selftest

run: build/sandpile
	./build/sandpile all --lines $(LINES) --cc "$(CC_VER)" --arch "$(ARCH)" --cores $(CORES)

lines:
	@wc -l Makefile $(SRCS) $(HDRS)
	@test $(LINES) -le 999

clean:
	rm -rf build
