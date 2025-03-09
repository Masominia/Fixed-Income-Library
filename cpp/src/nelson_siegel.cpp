#include "curvekit/nelson_siegel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "curvekit/linsolve.hpp"
#include "curvekit/root.hpp"

namespace ck {

static void loadings(double t, double tau, double& slope, double& curvature) {
    const double x = t / tau;
    slope = x < 1e-8 ? 1.0 - 0.5 * x : (1.0 - std::exp(-x)) / x;
    curvature = slope - std::exp(-x);
}

double NelsonSiegel::zero(double t) const {
    double s, c;
    loadings(t, tau, s, c);
    return level + slope * s + curvature * c;
}

static NelsonSiegelFit fit_for_tau(const Vector& times, const Vector& zeros, double tau) {
    std::vector<Vector> basis(3, Vector(times.size(), 1.0));
    for (std::size_t i = 0; i < times.size(); ++i) loadings(times[i], tau, basis[1][i], basis[2][i]);
    const Vector beta = least_squares(basis, zeros);
    NelsonSiegelFit fit;
    fit.model = {beta[0], beta[1], beta[2], tau};
    double sse = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
        const double e = fit.model.zero(times[i]) - zeros[i];
        sse += e * e;
        fit.max_error = std::max(fit.max_error, std::fabs(e));
    }
    fit.rmse = std::sqrt(sse / static_cast<double>(times.size()));
    return fit;
}

NelsonSiegelFit fit_nelson_siegel(const Vector& times, const Vector& zeros) {
    double best_tau = 1.0, best = std::numeric_limits<double>::infinity();
    for (double log_tau = std::log(0.1); log_tau <= std::log(20.0); log_tau += 0.1) {
        const double rmse = fit_for_tau(times, zeros, std::exp(log_tau)).rmse;
        if (rmse < best) best = rmse, best_tau = std::exp(log_tau);
    }
    const double refined = golden_minimum([&](double lt) { return fit_for_tau(times, zeros, std::exp(lt)).rmse; }, std::log(best_tau) - 0.1, std::log(best_tau) + 0.1);
    return fit_for_tau(times, zeros, std::exp(refined));
}

}  // namespace ck
