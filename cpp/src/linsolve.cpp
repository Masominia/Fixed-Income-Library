#include "curvekit/linsolve.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace ck {

Vector solve_linear(Vector a, Vector b) {
    const std::size_t n = b.size();
    if (a.size() != n * n) throw std::invalid_argument("solve_linear: matrix and right-hand side do not match");
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a[r * n + col]) > std::fabs(a[pivot * n + col])) pivot = r;
        if (std::fabs(a[pivot * n + col]) < 1e-300) throw std::runtime_error("solve_linear: singular matrix");
        if (pivot != col) {
            for (std::size_t k = 0; k < n; ++k) std::swap(a[col * n + k], a[pivot * n + k]);
            std::swap(b[col], b[pivot]);
        }
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = a[r * n + col] / a[col * n + col];
            for (std::size_t k = col; k < n; ++k) a[r * n + k] -= f * a[col * n + k];
            b[r] -= f * b[col];
        }
    }
    Vector x(n);
    for (std::size_t i = n; i-- > 0;) {
        double s = b[i];
        for (std::size_t k = i + 1; k < n; ++k) s -= a[i * n + k] * x[k];
        x[i] = s / a[i * n + i];
    }
    return x;
}

Vector least_squares(const std::vector<Vector>& basis, const Vector& target) {
    const std::size_t m = basis.size();
    Vector a(m * m, 0.0), b(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        if (basis[i].size() != target.size()) throw std::invalid_argument("least_squares: column length differs from target");
        for (std::size_t j = 0; j < m; ++j)
            for (std::size_t k = 0; k < target.size(); ++k) a[i * m + j] += basis[i][k] * basis[j][k];
        for (std::size_t k = 0; k < target.size(); ++k) b[i] += basis[i][k] * target[k];
    }
    return solve_linear(a, b);
}

}  // namespace ck
