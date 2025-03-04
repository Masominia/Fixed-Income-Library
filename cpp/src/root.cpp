#include "curvekit/root.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace ck {

double brent(const std::function<double(double)>& f, double a, double b, double tolerance) {
    double fa = f(a), fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if ((fa > 0.0) == (fb > 0.0)) throw std::invalid_argument("brent: root not bracketed");
    if (std::fabs(fa) < std::fabs(fb)) std::swap(a, b), std::swap(fa, fb);
    double c = a, fc = fa, d = b - a;
    bool bisected = true;
    for (int it = 0; it < 200; ++it) {
        if (fb == 0.0 || std::fabs(b - a) <= tolerance * (1.0 + std::fabs(b))) return b;
        double s;
        if (fa != fc && fb != fc) {
            s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }
        const double mid = 0.5 * (a + b);
        const bool outside = (s - mid) * (s - b) > 0.0 || !(std::fabs(s - b) < std::fabs((3.0 * a + b) / 4.0 - b) + std::fabs(b - a));
        const bool slow = bisected ? std::fabs(s - b) >= 0.5 * std::fabs(b - c) : std::fabs(s - b) >= 0.5 * std::fabs(c - d);
        if (outside || slow || !std::isfinite(s)) {
            s = mid;
            bisected = true;
        } else {
            bisected = false;
        }
        const double fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if ((fa > 0.0) != (fs > 0.0)) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if (std::fabs(fa) < std::fabs(fb)) std::swap(a, b), std::swap(fa, fb);
    }
    return b;
}

double golden_minimum(const std::function<double(double)>& f, double lo, double hi, double tolerance) {
    const double ratio = 0.6180339887498949;
    double x1 = hi - ratio * (hi - lo), x2 = lo + ratio * (hi - lo);
    double f1 = f(x1), f2 = f(x2);
    while (hi - lo > tolerance * (1.0 + std::fabs(lo) + std::fabs(hi))) {
        if (f1 < f2) {
            hi = x2, x2 = x1, f2 = f1;
            x1 = hi - ratio * (hi - lo);
            f1 = f(x1);
        } else {
            lo = x1, x1 = x2, f1 = f2;
            x2 = lo + ratio * (hi - lo);
            f2 = f(x2);
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace ck
