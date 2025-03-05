#include "curvekit/schedule.hpp"

#include <algorithm>
#include <stdexcept>

namespace ck {

Vector schedule(double start, double end, int frequency) {
    if (frequency <= 0 || !(end > start)) throw std::invalid_argument("schedule: need a positive frequency and end after start");
    Vector times;
    for (int k = 0;; ++k) {
        const double t = end - static_cast<double>(k) / frequency;
        if (t <= start + 1e-9) break;
        times.push_back(t);
    }
    std::reverse(times.begin(), times.end());
    return times;
}

}  // namespace ck
