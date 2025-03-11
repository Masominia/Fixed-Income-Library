#include "check.hpp"
#include "curvekit/bond.hpp"

using namespace ck;

void test_bond(Suite& suite) {
    const Bond bond{10.0, 0.05, 2, 100.0};
    CHECK(cashflows(bond).size() == 20);
    CHECK_NEAR(cashflows(bond).back().amount, 102.5, 1e-13);
    // at a yield equal to the coupon a bond on a coupon date is worth par
    CHECK_NEAR(price_from_yield(bond, 0.05), 100.0, 1e-11);
    CHECK_NEAR(yield_to_maturity(bond, 100.0), 0.05, 1e-12);
    CHECK_NEAR(yield_to_maturity(bond, price_from_yield(bond, 0.0321)), 0.0321, 1e-12);
    CHECK(price_from_yield(bond, 0.04) > 100.0 && price_from_yield(bond, 0.06) < 100.0);

    // a zero-coupon bond: Macaulay duration is its maturity
    const Bond zero{7.0, 0.0, 1, 100.0};
    const YieldAnalytics z = yield_analytics(zero, 0.03);
    CHECK_NEAR(z.macaulay, 7.0, 1e-13);
    CHECK_NEAR(z.modified, 7.0 / 1.03, 1e-13);
    CHECK_NEAR(z.convexity, 7.0 * 8.0 / (1.03 * 1.03), 1e-12);

    // analytic derivatives against finite differences of the price
    const YieldAnalytics a = yield_analytics(bond, 0.045);
    const double h = 1e-5, p = price_from_yield(bond, 0.045), up = price_from_yield(bond, 0.045 + h), dn = price_from_yield(bond, 0.045 - h);
    CHECK_NEAR(a.modified, -(up - dn) / (2.0 * h * p), 1e-6);
    CHECK_NEAR(a.convexity, (up - 2.0 * p + dn) / (h * h * p), 1e-3);
    CHECK_NEAR(a.dv01, p - price_from_yield(bond, 0.0451), 1e-3);
    CHECK(a.macaulay > 7.0 && a.macaulay < 10.0);

    // on a flat curve the price is the sum of discounted flows and the key-rate durations add up to the total
    const Curve flat = Curve::flat(0.04);
    double direct = 0.0;
    for (const Cashflow& c : cashflows(bond)) direct += c.amount * std::exp(-0.04 * c.time);
    CHECK_NEAR(price(bond, flat), direct, 1e-12);
    const Curve c({1.0, 2.0, 5.0, 10.0, 20.0}, {0.02, 0.025, 0.03, 0.034, 0.036}, Interpolation::LinearZero);
    const Vector krd = key_rate_durations(bond, c);
    double sum = 0.0;
    for (double k : krd) sum += k;
    CHECK_NEAR(sum, effective_duration(bond, c), 1e-6);
    CHECK(krd[3] > krd[2] && krd[2] > krd[1] && krd[4] == 0.0);  // the 10-year pillar carries most of a 10-year bond
}
