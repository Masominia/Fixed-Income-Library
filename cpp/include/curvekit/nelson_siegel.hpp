// The Nelson-Siegel curve: level, slope and curvature with one decay time.
#pragma once

#include "curvekit/types.hpp"

namespace ck {

struct NelsonSiegel {
    double level = 0.0, slope = 0.0, curvature = 0.0, tau = 1.0;
    double zero(double t) const;  // continuously compounded zero rate
};

struct NelsonSiegelFit {
    NelsonSiegel model;
    double rmse = 0.0;       // root mean square error of the zero rates
    double max_error = 0.0;
};

// Least squares on zero rates: for a given decay time the three coefficients are a linear problem,
// the decay time itself is found by a scan followed by a golden-section search.
NelsonSiegelFit fit_nelson_siegel(const Vector& times, const Vector& zero_rates);

}  // namespace ck
