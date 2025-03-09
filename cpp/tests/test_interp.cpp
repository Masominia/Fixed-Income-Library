#include "check.hpp"
#include "curvekit/interp.hpp"

using namespace ck;

void test_interp(Suite& suite) {
    const Vector x = {0.0, 1.0, 2.0, 4.0};
    CHECK(locate(x, -1.0) == 0 && locate(x, 0.5) == 0 && locate(x, 1.0) == 1 && locate(x, 3.9) == 2 && locate(x, 9.0) == 2);

    // a straight line is reproduced exactly
    const Vector line = {1.0, 3.0, 5.0, 9.0};
    const Vector m = monotone_slopes(x, line);
    for (double s : m) CHECK_NEAR(s, 2.0, 1e-14);
    double d = 0.0;
    CHECK_NEAR(hermite(2.0, 4.0, 5.0, 9.0, m[2], m[3], 3.3, &d), 7.6, 1e-14);
    CHECK_NEAR(d, 2.0, 1e-14);

    // monotone data with a sharp bend: no overshoot anywhere
    const Vector bend = {0.0, 0.0, 0.1, 5.0};
    const Vector mb = monotone_slopes(x, bend);
    bool monotone = true;
    double previous = 0.0;
    for (double v = 0.0; v <= 4.0; v += 0.01) {
        const std::size_t i = locate(x, v);
        const double y = hermite(x[i], x[i + 1], bend[i], bend[i + 1], mb[i], mb[i + 1], v);
        if (y < previous - 1e-14 || y > 5.0 + 1e-14) monotone = false;
        previous = y;
    }
    CHECK(monotone);
    // a local extremum in the data gets a zero slope
    CHECK(monotone_slopes({0.0, 1.0, 2.0}, {0.0, 1.0, 0.0})[1] == 0.0);

    // the derivative returned by the segment agrees with a finite difference
    const double h = 1e-6;
    const double up = hermite(0.0, 2.0, 1.0, 4.0, 0.5, 3.0, 0.7 + h), dn = hermite(0.0, 2.0, 1.0, 4.0, 0.5, 3.0, 0.7 - h);
    hermite(0.0, 2.0, 1.0, 4.0, 0.5, 3.0, 0.7, &d);
    CHECK_NEAR(d, (up - dn) / (2.0 * h), 1e-8);
}
