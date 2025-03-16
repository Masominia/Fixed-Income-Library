#include "curvekit/montecarlo.hpp"

#include <algorithm>
#include <cmath>

#include "curvekit/rng.hpp"

namespace ck {

template <class Step>
static Estimate average_discount(double r0, double maturity, int steps, int paths, std::uint64_t seed, Step step) {
    Rng rng(seed);
    const double dt = maturity / steps;
    double sum = 0.0, sum2 = 0.0;
    for (int p = 0; p < paths; ++p) {
        // antithetic pair
        double ra = r0, rb = r0, ia = 0.0, ib = 0.0;
        for (int i = 0; i < steps; ++i) {
            const double z = rng.normal();
            const double na = step(ra, z), nb = step(rb, -z);
            ia += 0.5 * (ra + na) * dt;
            ib += 0.5 * (rb + nb) * dt;
            ra = na;
            rb = nb;
        }
        const double v = 0.5 * (std::exp(-ia) + std::exp(-ib));
        sum += v;
        sum2 += v * v;
    }
    Estimate e;
    e.value = sum / paths;
    e.std_error = std::sqrt(std::max(sum2 / paths - e.value * e.value, 0.0) / paths);
    return e;
}

Estimate mc_bond(const Vasicek& m, double maturity, int steps, int paths, std::uint64_t seed) {
    const double dt = maturity / steps, decay = std::exp(-m.speed * dt);
    const double sd = m.sigma * std::sqrt((1.0 - decay * decay) / (2.0 * m.speed));
    return average_discount(m.r0, maturity, steps, paths, seed, [&](double r, double z) { return m.mean + (r - m.mean) * decay + sd * z; });
}

Estimate mc_bond(const CoxIngersollRoss& m, double maturity, int steps, int paths, std::uint64_t seed) {
    const double dt = maturity / steps, root = std::sqrt(dt);
    return average_discount(m.r0, maturity, steps, paths, seed, [&](double r, double z) {
        const double positive = std::max(r, 0.0);
        return r + m.speed * (m.mean - positive) * dt + m.sigma * std::sqrt(positive) * root * z;
    });
}

}  // namespace ck
