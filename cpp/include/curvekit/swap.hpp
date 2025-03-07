// Fixed-for-floating interest rate swaps on a single curve.
#pragma once

#include "curvekit/curve.hpp"

namespace ck {

struct Swap {
    double start = 0.0;
    double maturity = 0.0;
    double fixed_rate = 0.0;
    int frequency = 1;    // fixed payments per year
    bool payer = true;    // pays fixed, receives floating
    double notional = 1.0;
};

// Sum of accrual fraction times discount factor over the fixed payments.
double annuity(const Curve& curve, double start, double maturity, int frequency);

// Fixed rate that makes the swap worth zero.
double par_rate(const Curve& curve, double start, double maturity, int frequency);

double swap_value(const Curve& curve, const Swap& swap);

// Value change for a one basis point rise of every zero rate.
double swap_dv01(const Curve& curve, const Swap& swap);

}  // namespace ck
