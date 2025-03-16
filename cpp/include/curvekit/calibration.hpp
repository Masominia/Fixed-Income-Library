// Calibration of the two Hull-White parameters to at-the-money swaption volatilities.
#pragma once

#include <vector>

#include "curvekit/curve.hpp"

namespace ck {

struct SwaptionQuote {
    double expiry = 1.0;
    double tenor = 1.0;
    double black_vol = 0.2;
};

struct Calibration {
    double mean_reversion = 0.0;
    double sigma = 0.0;
    double rmse_vol = 0.0;  // root mean square difference between model and quoted Black volatilities
    int evaluations = 0;
};

// Black volatility implied by the model price of an at-the-money payer swaption.
double model_black_vol(const Curve& curve, double mean_reversion, double sigma, double expiry, double tenor, int frequency);

// Minimises the relative price errors over (mean reversion, sigma) with a scan followed by a pattern search in logarithms.
Calibration calibrate_hull_white(const Curve& curve, const std::vector<SwaptionQuote>& quotes, int frequency);

}  // namespace ck
