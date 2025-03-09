#include "curvekit/bootstrap.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "curvekit/root.hpp"
#include "curvekit/swap.hpp"

namespace ck {

double implied_rate(const Quote& q, const Curve& curve) {
    switch (q.kind) {
        case Quote::Kind::Deposit: return (1.0 / curve.discount(q.maturity) - 1.0) / q.maturity;
        case Quote::Kind::Fra: return curve.forward(q.start, q.maturity);
        default: return par_rate(curve, 0.0, q.maturity, q.frequency);
    }
}

Bootstrap bootstrap(std::vector<Quote> quotes, Interpolation mode) {
    if (quotes.empty()) throw std::invalid_argument("bootstrap: no quotes");
    std::sort(quotes.begin(), quotes.end(), [](const Quote& a, const Quote& b) { return a.maturity < b.maturity; });
    Vector times, zeros;
    for (const Quote& q : quotes) {
        if (!times.empty() && q.maturity <= times.back() + 1e-9) throw std::invalid_argument("bootstrap: two quotes share a maturity");
        times.push_back(q.maturity);
        zeros.push_back(q.rate);
    }
    Bootstrap out{Curve(times, zeros, mode), 0, 0.0};
    for (out.passes = 1; out.passes <= 50; ++out.passes) {
        double moved = 0.0;
        for (std::size_t k = 0; k < quotes.size(); ++k) {
            const double before = out.curve.zero_rates()[k];
            const double solved = brent(
                [&](double z) {
                    out.curve.set_zero(k, z);
                    return implied_rate(quotes[k], out.curve) - quotes[k].rate;
                },
                -0.5, 1.0, 1e-15);
            out.curve.set_zero(k, solved);
            moved = std::max(moved, std::fabs(solved - before));
        }
        if (moved < 1e-14) break;
    }
    for (const Quote& q : quotes) out.max_error = std::max(out.max_error, std::fabs(implied_rate(q, out.curve) - q.rate));
    return out;
}

}  // namespace ck
