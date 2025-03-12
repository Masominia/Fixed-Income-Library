// Two classical one-factor short-rate models with closed-form bond prices.
#pragma once

namespace ck {

// dr = speed (mean - r) dt + sigma dW
struct Vasicek {
    double speed = 0.1, mean = 0.03, sigma = 0.01, r0 = 0.03;
    double bond(double maturity) const;
    double zero(double maturity) const;
    // European option expiring at `expiry` on the bond maturing at `maturity`.
    double bond_option(double expiry, double maturity, double strike, bool call) const;
};

// dr = speed (mean - r) dt + sigma sqrt(r) dW
struct CoxIngersollRoss {
    double speed = 0.1, mean = 0.03, sigma = 0.05, r0 = 0.03;
    double bond(double maturity) const;
    double zero(double maturity) const;
    // True if 2 speed mean >= sigma^2, in which case the rate never reaches zero.
    bool feller() const { return 2.0 * speed * mean >= sigma * sigma; }
};

}  // namespace ck
