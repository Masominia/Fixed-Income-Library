// Hull-White trinomial tree fitted to the initial curve, for products without a closed form.
#pragma once

#include "curvekit/curve.hpp"

namespace ck {

class HullWhiteTree {
  public:
    // The tree covers [0, horizon] with `steps_per_year` equal steps per year.
    HullWhiteTree(const Curve& curve, double mean_reversion, double sigma, double horizon, int steps_per_year);

    int steps() const { return steps_; }
    int width() const { return jmax_; }

    // Price today of one unit paid at `maturity`; equals the curve's discount factor up to rounding.
    double discount_bond(double maturity) const;

    double bond_option(double expiry, double maturity, double strike, bool call) const;

    // Swaption with exercise at `expiry` only, or, if `bermudan`, also at every later fixed payment time before maturity.
    double swaption(double expiry, double tenor, int frequency, double strike, bool payer, bool bermudan) const;

  private:
    int index(double time) const;
    int reach(int step) const { return step < jmax_ ? step : jmax_; }
    // Replaces the values at step + 1 by the discounted expectations at `step`.
    void rollback(Vector& values, int step) const;

    double dt_, dx_;
    int steps_, jmax_;
    Vector alpha_;               // shift of the short rate at each step
    Vector up_, mid_, down_;     // branching probabilities per node index j + jmax
};

}  // namespace ck
