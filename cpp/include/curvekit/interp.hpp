// Interpolation helpers used by the curve.
#pragma once

#include "curvekit/types.hpp"

namespace ck {

// Index i of the interval [x[i], x[i+1]] that contains v, clamped to the first or last interval.
std::size_t locate(const Vector& x, double v);

// Node slopes of the shape-preserving cubic of Fritsch and Carlson: monotone data give a monotone interpolant.
Vector monotone_slopes(const Vector& x, const Vector& y);

// Cubic Hermite segment through (x0, y0) and (x1, y1) with slopes m0 and m1. If `derivative` is given it receives dy/dx.
double hermite(double x0, double x1, double y0, double y1, double m0, double m1, double x, double* derivative = nullptr);

}  // namespace ck
