#include "check.hpp"
#include "curvekit/calibration.hpp"

using namespace ck;

void test_calibration(Suite& suite) {
    const Curve c({0.5, 1.0, 2.0, 5.0, 10.0, 20.0}, {0.030, 0.032, 0.034, 0.037, 0.040, 0.041});
    std::vector<SwaptionQuote> quotes;
    for (double e : {1.0, 3.0, 7.0})
        for (double t : {1.0, 5.0, 10.0}) quotes.push_back({e, t, model_black_vol(c, 0.08, 0.012, e, t, 2)});
    // volatilities fall with expiry and tenor when rates mean-revert
    CHECK(quotes[0].black_vol > quotes[2].black_vol && quotes[0].black_vol > quotes[6].black_vol);
    const Calibration fit = calibrate_hull_white(c, quotes, 2);
    CHECK_NEAR(fit.mean_reversion, 0.08, 1e-5);
    CHECK_NEAR(fit.sigma, 0.012, 1e-7);
    CHECK(fit.rmse_vol < 1e-7 && fit.evaluations > 100 && fit.evaluations < 5000);
    // quotes that no parameter pair can match are fitted approximately
    quotes[4].black_vol *= 1.10;
    const Calibration rough = calibrate_hull_white(c, quotes, 2);
    CHECK(rough.rmse_vol > 1e-3 && rough.rmse_vol < 0.02);
    CHECK_NEAR(rough.sigma, 0.012, 0.002);
}
