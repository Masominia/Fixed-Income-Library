// A very small test harness: CHECK records a failure and carries on, each test file registers one function.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

struct Suite {
    int checks = 0;
    int failures = 0;
    void record(bool ok, const char* expr, const char* file, int line) {
        ++checks;
        if (!ok) {
            ++failures;
            std::printf("FAILED %s:%d  %s\n", file, line, expr);
        }
    }
};

#define CHECK(cond) suite.record((cond), #cond, __FILE__, __LINE__)
#define CHECK_NEAR(x, y, tol) suite.record(std::fabs((x) - (y)) <= (tol), #x " ~ " #y, __FILE__, __LINE__)
#define CHECK_THROWS(expr)                                 \
    do {                                                   \
        bool threw_ = false;                               \
        try {                                              \
            (void)(expr);                                  \
        } catch (const std::exception&) {                  \
            threw_ = true;                                 \
        }                                                  \
        suite.record(threw_, #expr " throws", __FILE__, __LINE__); \
    } while (0)

void test_root(Suite& suite);
void test_linsolve(Suite& suite);
void test_interp(Suite& suite);
void test_curve(Suite& suite);
void test_swap(Suite& suite);
void test_bootstrap(Suite& suite);
void test_bond(Suite& suite);
void test_nelson_siegel(Suite& suite);
void test_black(Suite& suite);
void test_shortrate(Suite& suite);
void test_hullwhite(Suite& suite);
void test_hwtree(Suite& suite);
void test_calibration(Suite& suite);
void test_quotes(Suite& suite);
