#include "curvekit/black.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "curvekit/normal.hpp"
#include "curvekit/root.hpp"
#include "curvekit/schedule.hpp"
#include "curvekit/swap.hpp"

namespace ck {

double black(double f, double k, double vol, double expiry, bool call) {
    const double sd = vol * std::sqrt(expiry), w = call ? 1.0 : -1.0;
    if (sd <= 0.0) return std::max(w * (f - k), 0.0);
    const double d1 = std::log(f / k) / sd + 0.5 * sd;
    return w * (f * normal_cdf(w * d1) - k * normal_cdf(w * (d1 - sd)));
}

double bachelier(double f, double k, double vol, double expiry, bool call) {
    const double sd = vol * std::sqrt(expiry), w = call ? 1.0 : -1.0;
    if (sd <= 0.0) return std::max(w * (f - k), 0.0);
    const double d = (f - k) / sd;
    return w * (f - k) * normal_cdf(w * d) + sd * normal_pdf(d);
}

double implied_black_vol(double price, double f, double k, double expiry, bool call) {
    const double intrinsic = std::max(call ? f - k : k - f, 0.0), upper = call ? f : k;
    if (!(price > intrinsic && price < upper)) throw std::invalid_argument("implied_black_vol: price outside the no-arbitrage bounds");
    return brent([&](double v) { return black(f, k, v, expiry, call) - price; }, 1e-8, 20.0, 1e-14);
}

double black_caplet(const Curve& curve, double start, double end, double strike, double vol) {
    return curve.discount(end) * (end - start) * black(curve.forward(start, end), strike, vol, start, true);
}

double black_cap(const Curve& curve, double maturity, int frequency, double strike, double vol) {
    const Vector times = schedule(0.0, maturity, frequency);
    double value = 0.0;
    for (std::size_t i = 1; i < times.size(); ++i) value += black_caplet(curve, times[i - 1], times[i], strike, vol);
    return value;
}

double black_swaption(const Curve& curve, double expiry, double tenor, int frequency, double strike, double vol, bool payer) {
    const double a = annuity(curve, expiry, expiry + tenor, frequency);
    return a * black(par_rate(curve, expiry, expiry + tenor, frequency), strike, vol, expiry, payer);
}

}  // namespace ck
