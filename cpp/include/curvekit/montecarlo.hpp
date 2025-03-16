// Monte Carlo estimates of discount factors under the short-rate models, as an independent check of the formulas.
#pragma once

#include <cstdint>

#include "curvekit/shortrate.hpp"

namespace ck {

struct Estimate {
    double value = 0.0;
    double std_error = 0.0;
};

// E[exp(-integral of r)] with the exact transition of the rate and the trapezoid rule for the integral.
Estimate mc_bond(const Vasicek& model, double maturity, int steps, int paths, std::uint64_t seed);

// Same with an Euler scheme in which a negative rate is replaced by zero inside the drift and the diffusion.
Estimate mc_bond(const CoxIngersollRoss& model, double maturity, int steps, int paths, std::uint64_t seed);

}  // namespace ck
