#include "check.hpp"
#include "curvekit/hullwhite.hpp"
#include "curvekit/hwtree.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

void test_hwtree(Suite& suite) {
    const Curve c({0.5, 1.0, 2.0, 5.0, 10.0}, {0.030, 0.032, 0.034, 0.037, 0.040});
    const HullWhite hw(c, 0.1, 0.01);
    const HullWhiteTree tree(c, 0.1, 0.01, 7.0, 32);
    CHECK(tree.steps() == 224 && tree.width() == 59);
    // the fit to the curve is exact by construction, also beyond the step where the tree stops widening
    for (double t : {0.5, 1.0, 2.5, 5.0, 7.0}) CHECK_NEAR(tree.discount_bond(t) / c.discount(t), 1.0, 1e-13);
    CHECK_THROWS(tree.discount_bond(0.51));
    CHECK_THROWS(tree.discount_bond(8.0));

    const double k = c.discount(7.0) / c.discount(2.0);
    CHECK_NEAR(tree.bond_option(2.0, 7.0, k, true) / hw.bond_option(2.0, 7.0, k, true), 1.0, 3e-3);
    // far out of the money the relative error is larger
    CHECK_NEAR(tree.bond_option(2.0, 7.0, 0.9 * k, false) / hw.bond_option(2.0, 7.0, 0.9 * k, false), 1.0, 5e-2);
    // put-call parity holds on the tree to rounding, whatever the discretisation error of each option
    CHECK_NEAR(tree.bond_option(2.0, 7.0, k, true) - tree.bond_option(2.0, 7.0, k, false), c.discount(7.0) - k * c.discount(2.0), 1e-13);

    const double atm = par_rate(c, 2.0, 7.0, 2);
    const double european = tree.swaption(2.0, 5.0, 2, atm, true, false);
    CHECK_NEAR(european / hw.swaption(2.0, 5.0, 2, atm, true), 1.0, 3e-3);
    CHECK_NEAR(european - tree.swaption(2.0, 5.0, 2, atm, false, false), 0.0, 1e-13);
    // the right to exercise later as well is worth something, for both sides
    CHECK(tree.swaption(2.0, 5.0, 2, atm, true, true) > european * 1.05);
    CHECK(tree.swaption(2.0, 5.0, 2, atm, false, true) > european * 1.02);
    // a Bermudan is worth at least each European option it contains
    const double bermudan = tree.swaption(2.0, 5.0, 2, atm, true, true);
    for (double later : {2.5, 4.0, 6.0}) CHECK(bermudan >= tree.swaption(later, 7.0 - later, 2, atm, true, false) - 1e-12);
    // with a single exercise date left, Bermudan and European coincide
    CHECK_NEAR(tree.swaption(6.5, 0.5, 2, atm, true, true), tree.swaption(6.5, 0.5, 2, atm, true, false), 1e-15);
    // finer trees get closer to the closed form
    const double coarse = HullWhiteTree(c, 0.1, 0.01, 7.0, 4).swaption(2.0, 5.0, 2, atm, true, false);
    CHECK(std::fabs(european - hw.swaption(2.0, 5.0, 2, atm, true)) < std::fabs(coarse - hw.swaption(2.0, 5.0, 2, atm, true)));
}
