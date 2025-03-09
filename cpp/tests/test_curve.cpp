#include "check.hpp"
#include "curvekit/curve.hpp"

using namespace ck;

void test_curve(Suite& suite) {
    const Curve flat = Curve::flat(0.03);
    CHECK_NEAR(flat.discount(2.0), std::exp(-0.06), 1e-15);
    CHECK_NEAR(flat.zero(7.5), 0.03, 1e-15);
    CHECK_NEAR(flat.instantaneous_forward(40.0), 0.03, 1e-15);
    CHECK_NEAR(flat.forward(1.0, 1.5), (std::exp(0.015) - 1.0) / 0.5, 1e-14);
    CHECK(flat.discount(0.0) == 1.0);

    const Vector t = {0.5, 1.0, 2.0, 5.0, 10.0}, z = {0.020, 0.024, 0.029, 0.034, 0.036};
    for (Interpolation mode : {Interpolation::LinearZero, Interpolation::LogLinear, Interpolation::MonotoneCubic}) {
        const Curve c(t, z, mode);
        for (std::size_t i = 0; i < t.size(); ++i) CHECK_NEAR(c.zero(t[i]), z[i], 1e-14);
        // the instantaneous forward is the derivative of -ln discount
        for (double u : {0.2, 0.75, 1.6, 3.3, 8.0, 12.0}) {
            const double h = 1e-5;
            CHECK_NEAR(c.instantaneous_forward(u), (c.log_discount(u + h) - c.log_discount(u - h)) / (2.0 * h), 1e-8);
        }
        // discount factors fall with time
        bool decreasing = true;
        for (double u = 0.0; u < 14.0; u += 0.05) decreasing = decreasing && c.discount(u + 0.05) < c.discount(u);
        CHECK(decreasing);
        CHECK_NEAR(c.shifted(0.001).zero(2.0), 0.030, 1e-14);
        CHECK_NEAR(c.bumped(2, 0.001).zero(2.0), 0.030, 1e-14);
        CHECK_NEAR(c.bumped(2, 0.001).zero(5.0), 0.034, 1e-14);
        CHECK(std::string(interpolation_name(mode)) == interpolation_name(parse_interpolation(interpolation_name(mode))));
    }
    // piecewise-constant forwards jump at the pillars, the cubic does not
    const Curve steps(t, z, Interpolation::LogLinear), smooth(t, z, Interpolation::MonotoneCubic);
    CHECK(std::fabs(steps.instantaneous_forward(2.0 + 1e-9) - steps.instantaneous_forward(2.0 - 1e-9)) > 1e-3);
    CHECK(std::fabs(smooth.instantaneous_forward(2.0 + 1e-9) - smooth.instantaneous_forward(2.0 - 1e-9)) < 1e-8);

    CHECK_THROWS(Curve({1.0, 1.0}, {0.02, 0.03}));
    CHECK_THROWS(Curve({1.0}, {0.02, 0.03}));
    CHECK_THROWS(flat.discount(-1.0));
    CHECK_THROWS(parse_interpolation("spline"));
}
