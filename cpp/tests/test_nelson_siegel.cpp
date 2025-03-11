#include "check.hpp"
#include "curvekit/nelson_siegel.hpp"

using namespace ck;

void test_nelson_siegel(Suite& suite) {
    const NelsonSiegel truth{0.045, -0.02, 0.015, 2.5};
    CHECK_NEAR(truth.zero(1e-12), 0.045 - 0.02, 1e-9);       // short end: level + slope
    CHECK_NEAR(truth.zero(1e6), 0.045, 1e-6);                // long end: level
    Vector t, z;
    for (double u : {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0}) {
        t.push_back(u);
        z.push_back(truth.zero(u));
    }
    const NelsonSiegelFit fit = fit_nelson_siegel(t, z);
    CHECK(fit.rmse < 1e-9);
    CHECK_NEAR(fit.model.tau, 2.5, 1e-4);
    CHECK_NEAR(fit.model.level, 0.045, 1e-7);
    CHECK_NEAR(fit.model.slope, -0.02, 1e-7);
    CHECK_NEAR(fit.model.curvature, 0.015, 1e-6);
    // a curve outside the family is fitted approximately, and the reported errors are consistent
    z[4] += 0.002;
    const NelsonSiegelFit rough = fit_nelson_siegel(t, z);
    CHECK(rough.rmse > 1e-4 && rough.rmse < 0.002 && rough.max_error >= rough.rmse);
}
