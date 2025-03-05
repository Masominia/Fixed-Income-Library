#include "check.hpp"
#include "curvekit/linsolve.hpp"

using namespace ck;

void test_linsolve(Suite& suite) {
    // needs a row exchange: the first pivot is zero
    const Vector x = solve_linear({0, 2, 1, 1, 1, 0, 2, 0, 3}, {7, 3, 11});
    CHECK_NEAR(x[0], 1.0, 1e-14);
    CHECK_NEAR(x[1], 2.0, 1e-14);
    CHECK_NEAR(x[2], 3.0, 1e-14);
    CHECK_THROWS(solve_linear({1, 2, 2, 4}, {1, 1}));
    CHECK_THROWS(solve_linear({1, 2, 3}, {1, 1}));

    Vector t(20), one(20, 1.0), y(20);
    for (int i = 0; i < 20; ++i) {
        t[i] = 0.1 * i;
        y[i] = 2.0 - 3.0 * t[i];
    }
    const Vector beta = least_squares({one, t}, y);
    CHECK_NEAR(beta[0], 2.0, 1e-12);
    CHECK_NEAR(beta[1], -3.0, 1e-12);
}
