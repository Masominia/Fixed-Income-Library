// Fixed-coupon bullet bonds: pricing on a curve and yield-based analytics.
#pragma once

#include "curvekit/curve.hpp"

namespace ck {

struct Bond {
    double maturity = 0.0;  // years from today
    double coupon = 0.0;    // annual rate
    int frequency = 2;      // coupons per year
    double face = 100.0;
};

struct Cashflow {
    double time = 0.0;
    double amount = 0.0;
};

// Remaining cash flows; the next coupon is paid in full (prices are dirty prices).
std::vector<Cashflow> cashflows(const Bond& bond);

double price(const Bond& bond, const Curve& curve);

// Price for a yield compounded `frequency` times per year.
double price_from_yield(const Bond& bond, double yield);
double yield_to_maturity(const Bond& bond, double price);

struct YieldAnalytics {
    double macaulay = 0.0;   // years
    double modified = 0.0;   // -dP/dy / P
    double convexity = 0.0;  // d2P/dy2 / P
    double dv01 = 0.0;       // price fall for a one basis point rise in yield
};

YieldAnalytics yield_analytics(const Bond& bond, double yield);

// -dP/dz / P for a parallel move of the curve's zero rates, by central difference.
double effective_duration(const Bond& bond, const Curve& curve, double bump = 1e-4);

// Sensitivity to each pillar's zero rate separately, same units as the effective duration.
Vector key_rate_durations(const Bond& bond, const Curve& curve, double bump = 1e-4);

}  // namespace ck
