CC      := cc
CXX     := c++
WARN    := -Wall -Wextra -Werror
OPT     := -O3 -fno-fast-math -ffp-contract=off
CFLAGS  := -std=c11 $(WARN) $(OPT)
CXXFLAGS:= -std=c++20 $(WARN) $(OPT)
LDFLAGS := -lm -lpthread

BUILD   := build
C_BIN   := $(BUILD)/ising_c
CPP_BIN := $(BUILD)/ising_cpp
C_TEST  := $(BUILD)/test_c
CPP_TEST:= $(BUILD)/test_cpp
C_ASAN  := $(BUILD)/test_c_asan
CPP_ASAN:= $(BUILD)/test_cpp_asan

.PHONY: all test asan clean parity

all: $(C_BIN) $(CPP_BIN) $(C_TEST) $(CPP_TEST)

$(BUILD):
	mkdir -p $(BUILD) data

$(C_BIN): c/main.c c/ising.c c/ising.h c/rng.h | $(BUILD)
	$(CC) $(CFLAGS) c/main.c c/ising.c -o $@ $(LDFLAGS)

$(CPP_BIN): cpp/main.cpp cpp/ising.cpp cpp/ising.hpp cpp/rng.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) cpp/main.cpp cpp/ising.cpp -o $@ $(LDFLAGS)

$(C_TEST): c/test_c.c c/ising.c c/ising.h c/rng.h | $(BUILD)
	$(CC) $(CFLAGS) c/test_c.c c/ising.c -o $@ $(LDFLAGS)

$(CPP_TEST): cpp/test_cpp.cpp cpp/ising.cpp cpp/ising.hpp cpp/rng.hpp | $(BUILD)
	$(CXX) $(CXXFLAGS) cpp/test_cpp.cpp cpp/ising.cpp -o $@ $(LDFLAGS)

test: $(C_TEST) $(CPP_TEST)
	./$(C_TEST)
	./$(CPP_TEST)

asan: | $(BUILD)
	$(CC) -std=c11 $(WARN) -O1 -g -fsanitize=address,undefined c/test_c.c c/ising.c -o $(C_ASAN) $(LDFLAGS)
	$(CXX) -std=c++20 $(WARN) -O1 -g -fsanitize=address,undefined cpp/test_cpp.cpp cpp/ising.cpp -o $(CPP_ASAN) $(LDFLAGS)
	./$(C_ASAN)
	./$(CPP_ASAN)

parity: $(C_BIN) $(CPP_BIN)
	./$(C_BIN) --mode parity --out-dir data
	./$(CPP_BIN) --mode parity --out-dir data

clean:
	rm -rf $(BUILD)
