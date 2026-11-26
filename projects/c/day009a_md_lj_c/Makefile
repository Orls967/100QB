CC ?= cc
CFLAGS = -std=c11 -O2 -Wall -Wextra -Werror -ffp-contract=off -fno-fast-math -Iinclude
LDFLAGS = -lm -lpthread
SCALE ?= full

SRCS = $(sort $(wildcard src/*.c))
HDRS = $(sort $(wildcard include/*.h))
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
LINES = $(shell wc -l $(HDRS) $(SRCS) | tail -n 1 | awk '{print $$1}')
CC_VER = $(shell $(CC) --version | head -n 1)
ARCH = $(shell uname -m)
CORES = $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

.PHONY: all test run lines clean

all: build/md

build/%.o: src/%.c $(HDRS) Makefile
	@mkdir -p build out ../../../outputs/day009a ../../../outputs/figures
	$(CC) $(CFLAGS) -c $< -o $@

build/md: $(OBJS)
	@mkdir -p build out ../../../outputs/day009a ../../../outputs/figures
	$(CC) $(CFLAGS) $(OBJS) -o build/md $(LDFLAGS)

test: build/md
	./build/md selftest

run: build/md
	./build/md all --scale $(SCALE) --lines $(LINES) --cc "$(CC_VER)" --arch "$(ARCH)" --cores $(CORES)

lines:
	@wc -l $(HDRS) $(SRCS)
	@echo "Total C lines (include/*.h + src/*.c): $(LINES)"
	@test $(LINES) -ge 6000 && test $(LINES) -le 9000

clean:
	rm -rf build
