// The Hull-White model, dr = (theta(t) - a r) dt + sigma dW, fitted to an initial curve: closed-form prices.
#pragma once

#include "curvekit/curve.hpp"

namespace ck {

class HullWhite {
  public:
    HullWhite(const Curve& curve, double mean_reversion, double sigma);

    double mean_reversion() const { return a_; }
    double sigma() const { return sigma_; }
    const Curve& curve() const { return curve_; }

    // Price at time t of the bond maturing at T when the short rate is r.
    double bond(double t, double maturity, double r) const;

    // European option expiring at `expiry` on the zero-coupon bond maturing at `maturity`.
    double bond_option(double expiry, double maturity, double strike, bool call) const;

    // Caplet on the simply compounded rate over [start, end]: a put on a zero-coupon bond.
    double caplet(double start, double end, double strike) const;
    double cap(double maturity, int frequency, double strike) const;

    // European swaption by Jamshidian's decomposition into options on zero-coupon bonds.
    double swaption(double expiry, double tenor, int frequency, double strike, bool payer) const;

  private:
    double loading(double t) const;
    Curve curve_;
    double a_, sigma_;
};

}  // namespace ck
