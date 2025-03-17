#include "check.hpp"
#include "curvekit/black.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

void test_black(Suite& suite) {
    // at the money: F (2 N(sd / 2) - 1)
    CHECK_NEAR(black(0.03, 0.03, 0.2, 4.0, true), 0.03 * (2.0 * 0.5792597094391030 - 1.0), 1e-15);
    // put-call parity, both models
    for (double k : {0.01, 0.03, 0.06}) {
        CHECK_NEAR(black(0.03, k, 0.3, 2.0, true) - black(0.03, k, 0.3, 2.0, false), 0.03 - k, 1e-15);
        CHECK_NEAR(bachelier(0.03, k, 0.008, 2.0, true) - bachelier(0.03, k, 0.008, 2.0, false), 0.03 - k, 1e-15);
    }
    CHECK_NEAR(bachelier(0.03, 0.03, 0.01, 1.0, true), 0.01 * 0.3989422804014327, 1e-15);
    CHECK_NEAR(black(0.03, 0.02, 0.0, 1.0, true), 0.01, 1e-16);
    CHECK(black(0.03, 0.02, 0.5, 0.0, false) == 0.0);
    // at the money and for small volatility the two models agree when normal vol = lognormal vol times forward
    CHECK_NEAR(black(0.03, 0.03, 0.1, 1.0, true) / bachelier(0.03, 0.03, 0.003, 1.0, true), 1.0, 1e-3);

    for (double vol : {0.05, 0.2, 0.8})
        for (double k : {0.015, 0.03, 0.05}) CHECK_NEAR(implied_black_vol(black(0.03, k, vol, 3.0, k > 0.03), 0.03, k, 3.0, k > 0.03), vol, 1e-9);
    CHECK_THROWS(implied_black_vol(0.05, 0.03, 0.03, 1.0, true));
    CHECK_THROWS(implied_black_vol(0.0, 0.03, 0.03, 1.0, true));

    const Curve c({1.0, 5.0, 10.0}, {0.02, 0.03, 0.035});
    // a cap is the sum of its caplets, and worth more with more volatility
    const double cap = black_cap(c, 2.0, 2, 0.03, 0.25);
    CHECK_NEAR(cap, black_caplet(c, 0.5, 1.0, 0.03, 0.25) + black_caplet(c, 1.0, 1.5, 0.03, 0.25) + black_caplet(c, 1.5, 2.0, 0.03, 0.25), 1e-15);
    CHECK(black_cap(c, 2.0, 2, 0.03, 0.35) > cap);
    // payer minus receiver swaption is the forward swap
    const double k = 0.03;
    const double gap = black_swaption(c, 2.0, 5.0, 2, k, 0.2, true) - black_swaption(c, 2.0, 5.0, 2, k, 0.2, false);
    CHECK_NEAR(gap, swap_value(c, {2.0, 7.0, k, 2, true, 1.0}), 1e-15);
}
