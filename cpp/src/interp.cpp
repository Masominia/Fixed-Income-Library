#include "curvekit/interp.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ck {

std::size_t locate(const Vector& x, double v) {
    if (x.size() < 2) throw std::invalid_argument("locate: need at least two nodes");
    const auto it = std::upper_bound(x.begin(), x.end(), v);
    const std::ptrdiff_t i = (it - x.begin()) - 1;
    return static_cast<std::size_t>(std::clamp<std::ptrdiff_t>(i, 0, static_cast<std::ptrdiff_t>(x.size()) - 2));
}

Vector monotone_slopes(const Vector& x, const Vector& y) {
    const std::size_t n = x.size();
    if (n < 2 || y.size() != n) throw std::invalid_argument("monotone_slopes: need matching arrays of at least two nodes");
    Vector secant(n - 1), m(n);
    for (std::size_t i = 0; i + 1 < n; ++i) secant[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    if (n == 2) return {secant[0], secant[0]};
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if (secant[i - 1] * secant[i] <= 0.0) {
            m[i] = 0.0;
        } else {
            // weighted harmonic mean of the neighbouring secants
            const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
            const double w0 = 2.0 * h1 + h0, w1 = h1 + 2.0 * h0;
            m[i] = (w0 + w1) / (w0 / secant[i - 1] + w1 / secant[i]);
        }
    }
    // one-sided three-point end slopes, limited so that the end segments stay monotone
    auto end_slope = [](double h0, double h1, double d0, double d1) {
        double s = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (s * d0 <= 0.0) return 0.0;
        if (d0 * d1 <= 0.0 && std::fabs(s) > 3.0 * std::fabs(d0)) return 3.0 * d0;
        return s;
    };
    m[0] = end_slope(x[1] - x[0], x[2] - x[1], secant[0], secant[1]);
    m[n - 1] = end_slope(x[n - 1] - x[n - 2], x[n - 2] - x[n - 3], secant[n - 2], secant[n - 3]);
    return m;
}

double hermite(double x0, double x1, double y0, double y1, double m0, double m1, double x, double* derivative) {
    const double h = x1 - x0, s = (x - x0) / h;
    const double h00 = (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s), h10 = s * (1.0 - s) * (1.0 - s);
    const double h01 = s * s * (3.0 - 2.0 * s), h11 = s * s * (s - 1.0);
    if (derivative) {
        const double d00 = 6.0 * s * (s - 1.0), d10 = (1.0 - s) * (1.0 - 3.0 * s), d01 = -d00, d11 = s * (3.0 * s - 2.0);
        *derivative = (d00 * y0 + d01 * y1) / h + d10 * m0 + d11 * m1;
    }
    return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1;
}

}  // namespace ck
