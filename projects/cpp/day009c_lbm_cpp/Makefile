CXX      := c++
CXXFLAGS := -std=c++20 -O2 -ffp-contract=off -Wall -Wextra -Werror
LDFLAGS  := -lpthread

BUILD    := build
BIN      := $(BUILD)/lbm
SRCS     := $(sort $(wildcard src/*.cpp))
HDRS     := $(sort $(wildcard include/lbm/*.hpp))
OBJS     := $(patsubst src/%.cpp,$(BUILD)/%.o,$(SRCS))

SCALE   ?= full
THREADS ?= 4

.PHONY: all test run clean lines

all: $(BIN)

$(BUILD):
	mkdir -p $(BUILD) out ../../../outputs/day009c ../../../outputs/figures

$(BUILD)/%.o: src/%.cpp $(HDRS) | $(BUILD)
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@

$(BIN): $(OBJS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

test: $(BIN)
	./$(BIN) selftest

run: $(BIN)
	./$(BIN) all --scale $(SCALE) --threads $(THREADS) --lines $$(wc -l $(HDRS) $(SRCS) | tail -n 1 | awk '{print $$1}') --compiler "$$($(CXX) --version | head -n 1)"

lines:
	@wc -l $(HDRS) $(SRCS)
	@total=$$(wc -l $(HDRS) $(SRCS) | tail -n 1 | awk '{print $$1}'); \
	echo "Total C++ lines (include/lbm/*.hpp + src/*.cpp): $$total"; \
	test "$$total" -ge 6000 -a "$$total" -le 9000 || (echo "Line count $$total outside [6000, 9000]" && exit 1)

clean:
	rm -rf $(BUILD)
