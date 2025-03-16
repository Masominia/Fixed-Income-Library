#include "curvekit/calibration.hpp"

#include <cmath>
#include <limits>

#include "curvekit/black.hpp"
#include "curvekit/hullwhite.hpp"
#include "curvekit/swap.hpp"

namespace ck {

double model_black_vol(const Curve& curve, double a, double sigma, double expiry, double tenor, int frequency) {
    const double rate = par_rate(curve, expiry, expiry + tenor, frequency), ann = annuity(curve, expiry, expiry + tenor, frequency);
    const double price = HullWhite(curve, a, sigma).swaption(expiry, tenor, frequency, rate, true);
    return implied_black_vol(price / ann, rate, rate, expiry, true);
}

Calibration calibrate_hull_white(const Curve& curve, const std::vector<SwaptionQuote>& quotes, int frequency) {
    Vector strikes, targets;
    for (const SwaptionQuote& q : quotes) {
        strikes.push_back(par_rate(curve, q.expiry, q.expiry + q.tenor, frequency));
        targets.push_back(black_swaption(curve, q.expiry, q.tenor, frequency, strikes.back(), q.black_vol, true));
    }
    Calibration out;
    auto objective = [&](double log_a, double log_s) {
        ++out.evaluations;
        const HullWhite model(curve, std::exp(log_a), std::exp(log_s));
        double sum = 0.0;
        for (std::size_t i = 0; i < quotes.size(); ++i) {
            const double e = model.swaption(quotes[i].expiry, quotes[i].tenor, frequency, strikes[i], true) / targets[i] - 1.0;
            sum += e * e;
        }
        return sum;
    };
    double best_a = 0.0, best_s = 0.0, best = std::numeric_limits<double>::infinity();
    for (double la = std::log(0.005); la <= std::log(0.5) + 1e-9; la += std::log(100.0) / 8.0)
        for (double ls = std::log(0.001); ls <= std::log(0.05) + 1e-9; ls += std::log(50.0) / 12.0) {
            const double v = objective(la, ls);
            if (v < best) best = v, best_a = la, best_s = ls;
        }
    double step = 0.25;
    while (step > 1e-7) {
        bool improved = false;
        const double moves[4][2] = {{step, 0.0}, {-step, 0.0}, {0.0, step}, {0.0, -step}};
        for (const auto& m : moves) {
            const double v = objective(best_a + m[0], best_s + m[1]);
            if (v < best) best = v, best_a += m[0], best_s += m[1], improved = true;
        }
        if (!improved) step *= 0.5;
    }
    out.mean_reversion = std::exp(best_a);
    out.sigma = std::exp(best_s);
    double sum = 0.0;
    for (const SwaptionQuote& q : quotes) {
        const double e = model_black_vol(curve, out.mean_reversion, out.sigma, q.expiry, q.tenor, frequency) - q.black_vol;
        sum += e * e;
    }
    out.rmse_vol = std::sqrt(sum / static_cast<double>(quotes.size()));
    return out;
}

}  // namespace ck
