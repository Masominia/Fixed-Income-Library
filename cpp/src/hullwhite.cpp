#include "curvekit/hullwhite.hpp"

#include <cmath>
#include <stdexcept>

#include "curvekit/normal.hpp"
#include "curvekit/root.hpp"
#include "curvekit/schedule.hpp"

namespace ck {

HullWhite::HullWhite(const Curve& curve, double mean_reversion, double sigma) : curve_(curve), a_(mean_reversion), sigma_(sigma) {
    if (!(a_ > 0.0) || !(sigma_ > 0.0)) throw std::invalid_argument("HullWhite: mean reversion and sigma must be positive");
}

double HullWhite::loading(double t) const { return (1.0 - std::exp(-a_ * t)) / a_; }

double HullWhite::bond(double t, double maturity, double r) const {
    const double b = loading(maturity - t);
    const double log_a = curve_.log_discount(t) - curve_.log_discount(maturity) + b * curve_.instantaneous_forward(t) -
                         sigma_ * sigma_ / (4.0 * a_) * (1.0 - std::exp(-2.0 * a_ * t)) * b * b;
    return std::exp(log_a - b * r);
}

double HullWhite::bond_option(double expiry, double maturity, double strike, bool call) const {
    const double sp = sigma_ * std::sqrt((1.0 - std::exp(-2.0 * a_ * expiry)) / (2.0 * a_)) * loading(maturity - expiry);
    const double long_bond = curve_.discount(maturity), short_bond = curve_.discount(expiry);
    const double h = std::log(long_bond / (strike * short_bond)) / sp + 0.5 * sp;
    const double w = call ? 1.0 : -1.0;
    return w * (long_bond * normal_cdf(w * h) - strike * short_bond * normal_cdf(w * (h - sp)));
}

double HullWhite::caplet(double start, double end, double strike) const {
    const double gross = 1.0 + strike * (end - start);
    return gross * bond_option(start, end, 1.0 / gross, false);
}

double HullWhite::cap(double maturity, int frequency, double strike) const {
    const Vector times = schedule(0.0, maturity, frequency);
    double value = 0.0;
    for (std::size_t i = 1; i < times.size(); ++i) value += caplet(times[i - 1], times[i], strike);
    return value;
}

double HullWhite::swaption(double expiry, double tenor, int frequency, double strike, bool payer) const {
    const Vector times = schedule(expiry, expiry + tenor, frequency);
    Vector amounts(times.size());
    double previous = expiry;
    for (std::size_t i = 0; i < times.size(); ++i) {
        amounts[i] = strike * (times[i] - previous);
        previous = times[i];
    }
    amounts.back() += 1.0;
    // the short rate at which the fixed leg is worth par; every bond price falls as the rate rises
    const double critical = brent(
        [&](double r) {
            double v = -1.0;
            for (std::size_t i = 0; i < times.size(); ++i) v += amounts[i] * bond(expiry, times[i], r);
            return v;
        },
        -1.0, 1.0, 1e-15);
    double value = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) value += amounts[i] * bond_option(expiry, times[i], bond(expiry, times[i], critical), !payer);
    return value;
}

}  // namespace ck
