// A discount curve defined by continuously compounded zero rates at pillar times (in years).
#pragma once

#include "curvekit/types.hpp"

namespace ck {

enum class Interpolation {
    LinearZero,     // zero rates linear in time
    LogLinear,      // log discount factors linear in time: piecewise-constant forward rates
    MonotoneCubic,  // shape-preserving cubic on log discount factors: continuous forward rates
};

Interpolation parse_interpolation(const char* name);
const char* interpolation_name(Interpolation mode);

class Curve {
  public:
    Curve(Vector times, Vector zero_rates, Interpolation mode = Interpolation::MonotoneCubic);
    static Curve flat(double rate);

    double discount(double t) const;
    double zero(double t) const;                    // continuously compounded
    double forward(double start, double end) const; // simply compounded over [start, end]
    double instantaneous_forward(double t) const;
    double log_discount(double t) const;            // -ln discount(t)

    const Vector& times() const { return t_; }
    const Vector& zero_rates() const { return z_; }
    Interpolation mode() const { return mode_; }

    void set_zero(std::size_t pillar, double rate);
    Curve shifted(double dz) const;                     // every zero rate moved by dz
    Curve bumped(std::size_t pillar, double dz) const;  // one pillar's zero rate moved by dz

  private:
    double evaluate(double t, double* derivative) const;
    void rebuild();

    Vector t_, z_;
    Vector x_, y_, m_;  // nodes of -ln discount, including the origin, and their slopes
    Interpolation mode_;
};

}  // namespace ck
