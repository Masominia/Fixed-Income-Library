// Black and Bachelier formulas for options on forward rates, and the market conventions built on them.
#pragma once

#include "curvekit/curve.hpp"

namespace ck {

// Undiscounted option on a lognormal forward. `call` pays max(F - K, 0).
double black(double forward, double strike, double vol, double expiry, bool call);

// Undiscounted option on a normally distributed forward; `vol` is in absolute rate units.
double bachelier(double forward, double strike, double vol, double expiry, bool call);

// Lognormal volatility that reproduces an undiscounted price. Throws if the price is outside the no-arbitrage bounds.
double implied_black_vol(double price, double forward, double strike, double expiry, bool call);

// Caplet on the simply compounded rate over [start, end], fixed at `start`, paid at `end`.
double black_caplet(const Curve& curve, double start, double end, double strike, double vol);

// Cap from today to `maturity` with the given payment frequency and one flat volatility; the first period is skipped.
double black_cap(const Curve& curve, double maturity, int frequency, double strike, double vol);

// Option expiring at `expiry` to enter a swap of length `tenor`.
double black_swaption(const Curve& curve, double expiry, double tenor, int frequency, double strike, double vol, bool payer);

}  // namespace ck
