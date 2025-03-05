// Payment times of a regular schedule, in years.
#pragma once

#include "curvekit/types.hpp"

namespace ck {

// Payment times in (start, end], spaced 1/frequency and generated backwards from `end`, so a short first period is possible.
Vector schedule(double start, double end, int frequency);

}  // namespace ck
