CXX      := c++
CXXFLAGS := -O3 -std=c++20 -Wall -Wextra -Werror -fno-fast-math -ffp-contract=off
LDFLAGS  := -lpthread

BUILD    := build
BIN      := $(BUILD)/wator
TEST_BIN := $(BUILD)/test_wator

.PHONY: all test run clean

all: $(BIN) $(TEST_BIN)

$(BUILD):
	mkdir -p $(BUILD) data

$(BIN): src/main.cpp src/wator.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) src/main.cpp -o $@ $(LDFLAGS)

$(TEST_BIN): tests/test_wator.cpp src/wator.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) tests/test_wator.cpp -o $@ $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

run: $(BIN)
	./$(BIN) --out-dir data --fig-dir ../../../outputs/figures

clean:
	rm -rf $(BUILD)
