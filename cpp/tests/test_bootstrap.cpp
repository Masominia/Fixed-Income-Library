#include "check.hpp"
#include "curvekit/bootstrap.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

static std::vector<Quote> market() {
    using K = Quote::Kind;
    return {{K::Deposit, 0.0, 0.5, 0.030, 1}, {K::Fra, 0.5, 1.0, 0.034, 1}, {K::Swap, 0.0, 2.0, 0.035, 2}, {K::Swap, 0.0, 5.0, 0.039, 2},
            {K::Swap, 0.0, 10.0, 0.041, 2}, {K::Swap, 0.0, 30.0, 0.040, 2}};
}

void test_bootstrap(Suite& suite) {
    for (Interpolation mode : {Interpolation::LinearZero, Interpolation::LogLinear, Interpolation::MonotoneCubic}) {
        const Bootstrap b = bootstrap(market(), mode);
        CHECK(b.max_error < 1e-13);
        CHECK(b.passes <= (mode == Interpolation::MonotoneCubic ? 20 : 2));
        CHECK(b.curve.times().size() == 6);
        for (const Quote& q : market()) CHECK_NEAR(implied_rate(q, b.curve), q.rate, 1e-13);
        CHECK_NEAR(b.curve.discount(0.5), 1.0 / 1.015, 1e-14);
        CHECK_NEAR(b.curve.discount(1.0), 1.0 / (1.015 * 1.017), 1e-14);
        // forwards stay positive on this market
        bool positive = true;
        for (double t = 0.0; t <= 30.0; t += 0.1) positive = positive && b.curve.instantaneous_forward(t) > 0.0;
        CHECK(positive);
    }
    // quotes generated from a known curve give that curve back at the pillars
    const Curve truth({1.0, 2.0, 3.0, 5.0}, {0.020, 0.025, 0.028, 0.031}, Interpolation::LogLinear);
    std::vector<Quote> q;
    for (double t : {1.0, 2.0, 3.0, 5.0}) q.push_back({Quote::Kind::Swap, 0.0, t, par_rate(truth, 0.0, t, 1), 1});
    const Bootstrap back = bootstrap(q, Interpolation::LogLinear);
    for (std::size_t i = 0; i < 4; ++i) CHECK_NEAR(back.curve.zero_rates()[i], truth.zero_rates()[i], 1e-13);
    // the order of the quotes does not matter
    std::reverse(q.begin(), q.end());
    CHECK_NEAR(bootstrap(q, Interpolation::LogLinear).curve.zero(3.0), 0.028, 1e-13);
    q.push_back(q.front());
    CHECK_THROWS(bootstrap(q, Interpolation::LogLinear));
    CHECK_THROWS(bootstrap({}, Interpolation::LogLinear));
}
