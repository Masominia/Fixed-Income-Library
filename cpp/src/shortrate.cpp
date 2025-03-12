#include "curvekit/shortrate.hpp"

#include <cmath>

#include "curvekit/normal.hpp"

namespace ck {

static double loading(double speed, double t) { return speed < 1e-10 ? t : (1.0 - std::exp(-speed * t)) / speed; }

double Vasicek::bond(double t) const {
    const double b = loading(speed, t);
    const double log_a = (b - t) * (speed * speed * mean - 0.5 * sigma * sigma) / (speed * speed) - sigma * sigma * b * b / (4.0 * speed);
    return std::exp(log_a - b * r0);
}

double Vasicek::zero(double t) const { return -std::log(bond(t)) / t; }

double Vasicek::bond_option(double expiry, double maturity, double strike, bool call) const {
    const double sp = sigma * std::sqrt((1.0 - std::exp(-2.0 * speed * expiry)) / (2.0 * speed)) * loading(speed, maturity - expiry);
    const double long_bond = bond(maturity), short_bond = bond(expiry);
    const double h = std::log(long_bond / (strike * short_bond)) / sp + 0.5 * sp;
    const double w = call ? 1.0 : -1.0;
    return w * (long_bond * normal_cdf(w * h) - strike * short_bond * normal_cdf(w * (h - sp)));
}

double CoxIngersollRoss::bond(double t) const {
    const double g = std::sqrt(speed * speed + 2.0 * sigma * sigma);
    const double e = std::expm1(g * t);
    const double denom = (g + speed) * e + 2.0 * g;
    const double b = 2.0 * e / denom;
    const double log_a = (2.0 * speed * mean / (sigma * sigma)) * (std::log(2.0 * g) + 0.5 * (speed + g) * t - std::log(denom));
    return std::exp(log_a - b * r0);
}

double CoxIngersollRoss::zero(double t) const { return -std::log(bond(t)) / t; }

}  // namespace ck
