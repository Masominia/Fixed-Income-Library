// One-dimensional root finding and minimisation.
#pragma once

#include <functional>

namespace ck {

// Brent's method on a bracket [lo, hi] with f(lo) f(hi) <= 0. Throws std::invalid_argument if the root is not bracketed.
double brent(const std::function<double(double)>& f, double lo, double hi, double tolerance = 1e-14);

// Golden-section search for a minimum of a unimodal function on [lo, hi].
double golden_minimum(const std::function<double(double)>& f, double lo, double hi, double tolerance = 1e-10);

}  // namespace ck
