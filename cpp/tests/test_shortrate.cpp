#include "check.hpp"
#include "curvekit/montecarlo.hpp"
#include "curvekit/shortrate.hpp"

using namespace ck;

void test_shortrate(Suite& suite) {
    const Vasicek v{0.3, 0.04, 0.012, 0.02};
    const CoxIngersollRoss c{0.3, 0.04, 0.08, 0.02};
    CHECK_NEAR(v.bond(0.0), 1.0, 1e-15);
    CHECK_NEAR(c.bond(0.0), 1.0, 1e-15);
    // short maturities yield the short rate; long ones a constant below the mean (convexity)
    CHECK_NEAR(v.zero(1e-6), 0.02, 1e-6);
    CHECK_NEAR(c.zero(1e-6), 0.02, 1e-6);
    CHECK_NEAR(v.zero(2000.0), 0.04 - 0.012 * 0.012 / (2.0 * 0.09), 1e-4);
    CHECK(c.zero(500.0) < 0.04 && c.zero(500.0) > 0.035 && c.feller());
    const CoxIngersollRoss wild{0.1, 0.01, 0.2, 0.02};
    CHECK(!wild.feller());
    // without volatility both reduce to the deterministic path of the rate
    const double t = 5.0, integral = 0.04 * t + (0.02 - 0.04) * (1.0 - std::exp(-0.3 * t)) / 0.3;
    CHECK_NEAR((Vasicek{0.3, 0.04, 1e-9, 0.02}.bond(t)), std::exp(-integral), 1e-12);
    CHECK_NEAR((CoxIngersollRoss{0.3, 0.04, 1e-5, 0.02}.bond(t)), std::exp(-integral), 1e-9);

    for (double m : {1.0, 5.0}) {
        const Estimate ev = mc_bond(v, m, static_cast<int>(50 * m), 20000, 3), ec = mc_bond(c, m, static_cast<int>(50 * m), 20000, 4);
        CHECK(std::fabs(ev.value - v.bond(m)) < 4.0 * ev.std_error + 2e-6);
        CHECK(std::fabs(ec.value - c.bond(m)) < 4.0 * ec.std_error + 1e-4);  // the Euler scheme has a small bias of its own
        CHECK(ev.std_error > 0.0 && ev.std_error < 1e-4);
    }
    // bond option: put-call parity and a positive time value
    const double k = v.bond(5.0) / v.bond(2.0);
    const double call = v.bond_option(2.0, 5.0, k, true), put = v.bond_option(2.0, 5.0, k, false);
    CHECK_NEAR(call - put, v.bond(5.0) - k * v.bond(2.0), 1e-15);
    CHECK(call > 0.0 && call < 0.05);
}
