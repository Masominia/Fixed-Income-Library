#include "curvekit/curve.hpp"

#include <cmath>
#include <cstring>
#include <stdexcept>

#include "curvekit/interp.hpp"

namespace ck {

Interpolation parse_interpolation(const char* name) {
    if (!std::strcmp(name, "linear_zero")) return Interpolation::LinearZero;
    if (!std::strcmp(name, "log_linear")) return Interpolation::LogLinear;
    if (!std::strcmp(name, "monotone_cubic")) return Interpolation::MonotoneCubic;
    throw std::invalid_argument("interpolation must be linear_zero, log_linear or monotone_cubic");
}

const char* interpolation_name(Interpolation mode) {
    switch (mode) {
        case Interpolation::LinearZero: return "linear_zero";
        case Interpolation::LogLinear: return "log_linear";
        default: return "monotone_cubic";
    }
}

Curve::Curve(Vector times, Vector zero_rates, Interpolation mode) : t_(std::move(times)), z_(std::move(zero_rates)), mode_(mode) {
    if (t_.empty() || t_.size() != z_.size()) throw std::invalid_argument("Curve: need one zero rate per pillar");
    for (std::size_t i = 0; i < t_.size(); ++i)
        if (!(t_[i] > (i == 0 ? 0.0 : t_[i - 1]))) throw std::invalid_argument("Curve: pillar times must be positive and increasing");
    rebuild();
}

Curve Curve::flat(double rate) { return Curve({1.0}, {rate}, Interpolation::LogLinear); }

void Curve::rebuild() {
    x_.assign(1, 0.0);
    y_.assign(1, 0.0);
    for (std::size_t i = 0; i < t_.size(); ++i) {
        x_.push_back(t_[i]);
        y_.push_back(z_[i] * t_[i]);
    }
    m_ = mode_ == Interpolation::MonotoneCubic ? monotone_slopes(x_, y_) : Vector();
}

double Curve::evaluate(double t, double* derivative) const {
    if (t < 0.0) throw std::invalid_argument("Curve: negative time");
    const std::size_t n = t_.size();
    if (mode_ == Interpolation::LinearZero) {
        double z = z_.back(), dz = 0.0;
        if (t <= t_.front()) {
            z = z_.front();
        } else if (t < t_.back()) {
            const std::size_t i = locate(t_, t);
            dz = (z_[i + 1] - z_[i]) / (t_[i + 1] - t_[i]);
            z = z_[i] + dz * (t - t_[i]);
        }
        if (derivative) *derivative = z + t * dz;
        return z * t;
    }
    if (t >= t_.back()) {  // flat forward beyond the last pillar
        const double slope = mode_ == Interpolation::MonotoneCubic ? m_.back() : (y_[n] - y_[n - 1]) / (x_[n] - x_[n - 1]);
        if (derivative) *derivative = slope;
        return y_.back() + slope * (t - t_.back());
    }
    const std::size_t i = locate(x_, t);
    if (mode_ == Interpolation::LogLinear) {
        const double slope = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
        if (derivative) *derivative = slope;
        return y_[i] + slope * (t - x_[i]);
    }
    return hermite(x_[i], x_[i + 1], y_[i], y_[i + 1], m_[i], m_[i + 1], t, derivative);
}

double Curve::log_discount(double t) const { return evaluate(t, nullptr); }

double Curve::discount(double t) const { return std::exp(-evaluate(t, nullptr)); }

double Curve::zero(double t) const { return t > 1e-12 ? evaluate(t, nullptr) / t : instantaneous_forward(0.0); }

double Curve::forward(double start, double end) const {
    if (!(end > start)) throw std::invalid_argument("Curve::forward: end must be after start");
    return (discount(start) / discount(end) - 1.0) / (end - start);
}

double Curve::instantaneous_forward(double t) const {
    double d = 0.0;
    evaluate(t, &d);
    return d;
}

void Curve::set_zero(std::size_t pillar, double rate) {
    z_.at(pillar) = rate;
    rebuild();
}

Curve Curve::shifted(double dz) const {
    Vector z = z_;
    for (double& v : z) v += dz;
    return Curve(t_, z, mode_);
}

Curve Curve::bumped(std::size_t pillar, double dz) const {
    Vector z = z_;
    z.at(pillar) += dz;
    return Curve(t_, z, mode_);
}

}  // namespace ck
