#include "check.hpp"
#include "curvekit/schedule.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

void test_swap(Suite& suite) {
    const Vector s = schedule(0.0, 2.0, 2);
    CHECK(s.size() == 4);
    CHECK_NEAR(s[0], 0.5, 1e-15);
    CHECK_NEAR(s[3], 2.0, 1e-15);
    const Vector stub = schedule(0.0, 1.25, 2);  // short first period
    CHECK(stub.size() == 3);
    CHECK_NEAR(stub[0], 0.25, 1e-15);
    CHECK(schedule(1.0, 3.0, 1).size() == 2);
    CHECK_THROWS(schedule(1.0, 1.0, 2));

    const Curve c({1.0, 5.0, 10.0}, {0.02, 0.03, 0.035});
    const double rate = par_rate(c, 0.0, 5.0, 2);
    CHECK_NEAR(swap_value(c, {0.0, 5.0, rate, 2, true, 1.0}), 0.0, 1e-15);
    // a payer and a receiver with the same terms cancel
    CHECK_NEAR(swap_value(c, {1.0, 6.0, 0.03, 2, true, 100.0}) + swap_value(c, {1.0, 6.0, 0.03, 2, false, 100.0}), 0.0, 1e-13);
    // on a flat curve with annual payments the par rate is the annually compounded rate
    CHECK_NEAR(par_rate(Curve::flat(0.04), 0.0, 7.0, 1), std::exp(0.04) - 1.0, 1e-14);
    CHECK_NEAR(par_rate(Curve::flat(0.04), 3.0, 7.0, 1), std::exp(0.04) - 1.0, 1e-14);
    CHECK_NEAR(annuity(Curve::flat(0.0), 2.0, 7.0, 4), 5.0, 1e-13);
    // receiving floating gains when rates rise; the size is close to annuity times one basis point
    const Swap payer{0.0, 5.0, rate, 2, true, 1.0};
    CHECK(swap_dv01(c, payer) > 0.0);
    CHECK_NEAR(swap_dv01(c, payer) / (1e-4 * annuity(c, 0.0, 5.0, 2)), 1.0, 0.05);
}
