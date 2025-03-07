#include "curvekit/swap.hpp"

#include "curvekit/schedule.hpp"

namespace ck {

double annuity(const Curve& curve, double start, double maturity, int frequency) {
    double a = 0.0, previous = start;
    for (double t : schedule(start, maturity, frequency)) {
        a += (t - previous) * curve.discount(t);
        previous = t;
    }
    return a;
}

double par_rate(const Curve& curve, double start, double maturity, int frequency) {
    return (curve.discount(start) - curve.discount(maturity)) / annuity(curve, start, maturity, frequency);
}

double swap_value(const Curve& curve, const Swap& s) {
    const double floating = curve.discount(s.start) - curve.discount(s.maturity);
    const double fixed = s.fixed_rate * annuity(curve, s.start, s.maturity, s.frequency);
    return s.notional * (s.payer ? floating - fixed : fixed - floating);
}

double swap_dv01(const Curve& curve, const Swap& s) { return swap_value(curve.shifted(1e-4), s) - swap_value(curve, s); }

}  // namespace ck
