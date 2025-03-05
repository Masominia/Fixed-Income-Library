#include "check.hpp"
#include "curvekit/root.hpp"

using namespace ck;

void test_root(Suite& suite) {
    CHECK_NEAR(brent([](double x) { return x * x * x - 2.0 * x - 5.0; }, 2.0, 3.0), 2.0945514815423265, 1e-13);
    CHECK_NEAR(brent([](double x) { return std::cos(x) - x; }, 0.0, 1.0), 0.7390851332151607, 1e-13);
    CHECK_NEAR(brent([](double x) { return std::exp(x) - 1e-8; }, -40.0, 5.0), std::log(1e-8), 1e-10);
    CHECK(brent([](double x) { return x - 1.0; }, 1.0, 2.0) == 1.0);
    CHECK_THROWS(brent([](double x) { return x * x + 1.0; }, -1.0, 1.0));
    int calls = 0;
    brent([&](double x) { return ++calls, std::pow(x, 9); }, -1.0, 1.5, 1e-12);
    CHECK(calls < 150);  // a flat root must not send the method into a long crawl
    CHECK_NEAR(golden_minimum([](double x) { return (x - 0.3) * (x - 0.3) + 2.0; }, -1.0, 2.0), 0.3, 1e-7);
}
