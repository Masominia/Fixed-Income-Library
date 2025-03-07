#include "curvekit/bond.hpp"

#include <cmath>

#include "curvekit/root.hpp"
#include "curvekit/schedule.hpp"

namespace ck {

std::vector<Cashflow> cashflows(const Bond& b) {
    std::vector<Cashflow> flows;
    const double coupon = b.face * b.coupon / b.frequency;
    for (double t : schedule(0.0, b.maturity, b.frequency)) flows.push_back({t, coupon});
    flows.back().amount += b.face;
    return flows;
}

double price(const Bond& b, const Curve& curve) {
    double p = 0.0;
    for (const Cashflow& c : cashflows(b)) p += c.amount * curve.discount(c.time);
    return p;
}

double price_from_yield(const Bond& b, double y) {
    double p = 0.0;
    for (const Cashflow& c : cashflows(b)) p += c.amount * std::pow(1.0 + y / b.frequency, -c.time * b.frequency);
    return p;
}

double yield_to_maturity(const Bond& b, double target) {
    return brent([&](double y) { return price_from_yield(b, y) - target; }, -0.9 * b.frequency, 10.0, 1e-15);
}

YieldAnalytics yield_analytics(const Bond& b, double y) {
    const double f = b.frequency, base = 1.0 + y / f;
    double p = 0.0, first = 0.0, second = 0.0;
    for (const Cashflow& c : cashflows(b)) {
        const double pv = c.amount * std::pow(base, -c.time * f);
        p += pv;
        first += c.time * pv;
        second += c.time * (c.time + 1.0 / f) * pv;
    }
    YieldAnalytics a;
    a.macaulay = first / p;
    a.modified = a.macaulay / base;
    a.convexity = second / (p * base * base);
    a.dv01 = a.modified * p * 1e-4;
    return a;
}

double effective_duration(const Bond& b, const Curve& curve, double bump) {
    return (price(b, curve.shifted(-bump)) - price(b, curve.shifted(bump))) / (2.0 * bump * price(b, curve));
}

Vector key_rate_durations(const Bond& b, const Curve& curve, double bump) {
    const double p = price(b, curve);
    Vector k(curve.times().size());
    for (std::size_t i = 0; i < k.size(); ++i) k[i] = (price(b, curve.bumped(i, -bump)) - price(b, curve.bumped(i, bump))) / (2.0 * bump * p);
    return k;
}

}  // namespace ck
