#include "check.hpp"
#include "curvekit/black.hpp"
#include "curvekit/hullwhite.hpp"
#include "curvekit/shortrate.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

void test_hullwhite(Suite& suite) {
    const Curve c({0.5, 1.0, 2.0, 5.0, 10.0, 20.0}, {0.030, 0.032, 0.034, 0.037, 0.040, 0.041});
    const HullWhite hw(c, 0.05, 0.01);
    // at time zero, with the short rate equal to the instantaneous forward, the model returns the curve
    for (double t : {0.5, 3.0, 12.0}) CHECK_NEAR(hw.bond(0.0, t, c.instantaneous_forward(0.0)), c.discount(t), 1e-14);
    CHECK(hw.bond(2.0, 7.0, 0.05) < hw.bond(2.0, 7.0, 0.03));

    // fitted to a Vasicek curve, Hull-White is that Vasicek model
    const Vasicek v{0.2, 0.05, 0.015, 0.03};
    Vector t, z;
    for (double u = 0.25; u <= 12.0; u += 0.25) {
        t.push_back(u);
        z.push_back(v.zero(u));
    }
    const HullWhite same(Curve(t, z), 0.2, 0.015);
    CHECK_NEAR(same.bond_option(2.0, 6.0, 0.85, true), v.bond_option(2.0, 6.0, 0.85, true), 1e-12);
    // bond prices at a future time agree too; this goes through the interpolated forward curve
    const Vasicek later{0.2, 0.05, 0.015, 0.045};
    CHECK_NEAR(same.bond(2.0, 6.0, 0.045) / later.bond(4.0), 1.0, 1e-4);

    // parities
    const double k = 0.85;
    CHECK_NEAR(hw.bond_option(2.0, 7.0, k, true) - hw.bond_option(2.0, 7.0, k, false), c.discount(7.0) - k * c.discount(2.0), 1e-15);
    for (double strike : {0.02, 0.04, 0.06}) {
        const double gap = hw.swaption(2.0, 5.0, 2, strike, true) - hw.swaption(2.0, 5.0, 2, strike, false);
        CHECK_NEAR(gap, swap_value(c, {2.0, 7.0, strike, 2, true, 1.0}), 1e-13);
    }
    // a swaption on a one-period swap is a caplet
    const double f = c.forward(2.0, 2.5);
    CHECK_NEAR(hw.swaption(2.0, 0.5, 2, f, true), hw.caplet(2.0, 2.5, f), 1e-14);
    CHECK_NEAR(hw.cap(1.5, 2, 0.035), hw.caplet(0.5, 1.0, 0.035) + hw.caplet(1.0, 1.5, 0.035), 1e-15);
    // more volatility, more value; more mean reversion, less
    const double atm = par_rate(c, 2.0, 7.0, 2);
    CHECK(HullWhite(c, 0.05, 0.012).swaption(2.0, 5.0, 2, atm, true) > hw.swaption(2.0, 5.0, 2, atm, true));
    CHECK(HullWhite(c, 0.20, 0.010).swaption(2.0, 5.0, 2, atm, true) < hw.swaption(2.0, 5.0, 2, atm, true));
    CHECK_THROWS(HullWhite(c, 0.0, 0.01));
}
