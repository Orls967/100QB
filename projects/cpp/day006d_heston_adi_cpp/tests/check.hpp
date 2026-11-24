#pragma once
#include <cstdio>
#include <cstdlib>
#include <vector>

struct TestReg {
    struct Entry { const char* name; void (*fn)(); };
    static std::vector<Entry>& all() { static std::vector<Entry> v; return v; }
    static int& fails() { static int f = 0; return f; }
};

struct RegAdder {
    RegAdder(const char* name, void (*fn)()) { TestReg::all().push_back({name, fn}); }
};

#define TEST_CASE(name) static void name(); static RegAdder reg_##name(#name, name); static void name()
#define CHECK(cond) do { if (!(cond)) { \
    std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
    ++TestReg::fails(); \
} } while (0)
