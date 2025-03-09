// Builds a curve that reprices a set of market quotes exactly.
#pragma once

#include <vector>

#include "curvekit/curve.hpp"

namespace ck {

struct Quote {
    enum class Kind { Deposit, Fra, Swap };
    Kind kind = Kind::Deposit;
    double start = 0.0;     // zero for deposits and swaps
    double maturity = 0.0;  // defines the pillar
    double rate = 0.0;
    int frequency = 1;      // fixed payments per year, swaps only
};

// Rate implied by the curve for the instrument of the quote.
double implied_rate(const Quote& quote, const Curve& curve);

struct Bootstrap {
    Curve curve;
    int passes = 0;          // sweeps over the instruments until nothing moved
    double max_error = 0.0;  // largest |implied rate - quoted rate|
};

// One pillar per quote, at its maturity. Each pillar is solved in turn with the others fixed, and the sweep is repeated
// until the zero rates stop changing; one sweep is enough for local interpolation, the cubic needs a few.
Bootstrap bootstrap(std::vector<Quote> quotes, Interpolation mode);

}  // namespace ck
