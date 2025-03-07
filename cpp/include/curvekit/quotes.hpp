// Text format for market quotes: one instrument per line, `kind start maturity rate frequency`.
#pragma once

#include <string>
#include <vector>

#include "curvekit/bootstrap.hpp"

namespace ck {

// Lines starting with # and empty lines are skipped. Kinds: deposit, fra, swap.
std::vector<Quote> read_quotes(const std::string& path);

}  // namespace ck
