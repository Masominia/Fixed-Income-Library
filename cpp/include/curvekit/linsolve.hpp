// Small dense linear systems and least squares.
#pragma once

#include "curvekit/types.hpp"

namespace ck {

// Solves a x = b for a square matrix stored row by row, by Gaussian elimination with partial pivoting.
// Throws std::runtime_error if the matrix is singular.
Vector solve_linear(Vector a, Vector b);

// Least-squares coefficients for columns `basis[k]` against `target` (normal equations; meant for a handful of columns).
Vector least_squares(const std::vector<Vector>& basis, const Vector& target);

}  // namespace ck
