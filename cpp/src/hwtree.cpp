#include "curvekit/hwtree.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "curvekit/schedule.hpp"

namespace ck {

HullWhiteTree::HullWhiteTree(const Curve& curve, double a, double sigma, double horizon, int steps_per_year) {
    if (!(a > 0.0) || !(sigma > 0.0) || steps_per_year <= 0 || !(horizon > 0.0)) throw std::invalid_argument("HullWhiteTree: invalid arguments");
    dt_ = 1.0 / steps_per_year;
    steps_ = static_cast<int>(std::lround(horizon * steps_per_year));
    dx_ = sigma * std::sqrt(3.0 * dt_);
    jmax_ = std::min(static_cast<int>(std::ceil(0.184 / (a * dt_))), steps_);
    jmax_ = std::max(jmax_, 1);
    const int n = 2 * jmax_ + 1;
    up_.resize(n), mid_.resize(n), down_.resize(n);
    for (int j = -jmax_; j <= jmax_; ++j) {
        const double m = a * j * dt_, m2 = m * m;
        const int k = j + jmax_;
        if (j == jmax_) {  // branches to j, j - 1, j - 2
            up_[k] = 7.0 / 6.0 + 0.5 * (m2 - 3.0 * m), mid_[k] = -1.0 / 3.0 - m2 + 2.0 * m, down_[k] = 1.0 / 6.0 + 0.5 * (m2 - m);
        } else if (j == -jmax_) {  // branches to j + 2, j + 1, j
            up_[k] = 1.0 / 6.0 + 0.5 * (m2 + m), mid_[k] = -1.0 / 3.0 - m2 - 2.0 * m, down_[k] = 7.0 / 6.0 + 0.5 * (m2 + 3.0 * m);
        } else {  // branches to j + 1, j, j - 1
            up_[k] = 1.0 / 6.0 + 0.5 * (m2 - m), mid_[k] = 2.0 / 3.0 - m2, down_[k] = 1.0 / 6.0 + 0.5 * (m2 + m);
        }
    }
    // forward induction on state prices: choose each shift so that the next discount factor is matched
    alpha_.resize(steps_);
    Vector q(n, 0.0), next(n);
    q[jmax_] = 1.0;
    for (int i = 0; i < steps_; ++i) {
        const int r = reach(i);
        double sum = 0.0;
        for (int j = -r; j <= r; ++j) sum += q[j + jmax_] * std::exp(-j * dx_ * dt_);
        alpha_[i] = (std::log(sum) + curve.log_discount((i + 1) * dt_)) / dt_;
        std::fill(next.begin(), next.end(), 0.0);
        for (int j = -r; j <= r; ++j) {
            const int k = j + jmax_;
            const double w = q[k] * std::exp(-(alpha_[i] + j * dx_) * dt_);
            const int centre = (j == jmax_ && i >= jmax_) ? k - 1 : ((j == -jmax_ && i >= jmax_) ? k + 1 : k);
            next[centre + 1] += w * up_[k];
            next[centre] += w * mid_[k];
            next[centre - 1] += w * down_[k];
        }
        q.swap(next);
    }
}

int HullWhiteTree::index(double time) const {
    const int i = static_cast<int>(std::lround(time / dt_));
    if (std::fabs(i * dt_ - time) > 1e-9 || i < 0 || i > steps_) throw std::invalid_argument("HullWhiteTree: time is not on the grid");
    return i;
}

void HullWhiteTree::rollback(Vector& v, int step) const {
    const int r = reach(step);
    Vector out(v.size(), 0.0);
    for (int j = -r; j <= r; ++j) {
        const int k = j + jmax_;
        const int centre = (j == jmax_ && step >= jmax_) ? k - 1 : ((j == -jmax_ && step >= jmax_) ? k + 1 : k);
        out[k] = std::exp(-(alpha_[step] + j * dx_) * dt_) * (up_[k] * v[centre + 1] + mid_[k] * v[centre] + down_[k] * v[centre - 1]);
    }
    v.swap(out);
}

double HullWhiteTree::discount_bond(double maturity) const {
    Vector v(2 * jmax_ + 1, 1.0);
    for (int i = index(maturity) - 1; i >= 0; --i) rollback(v, i);
    return v[jmax_];
}

double HullWhiteTree::bond_option(double expiry, double maturity, double strike, bool call) const {
    const int e = index(expiry);
    Vector v(2 * jmax_ + 1, 1.0);
    for (int i = index(maturity) - 1; i >= e; --i) rollback(v, i);
    for (double& x : v) x = std::max(call ? x - strike : strike - x, 0.0);
    for (int i = e - 1; i >= 0; --i) rollback(v, i);
    return v[jmax_];
}

double HullWhiteTree::swaption(double expiry, double tenor, int frequency, double strike, bool payer, bool bermudan) const {
    const Vector times = schedule(expiry, expiry + tenor, frequency);
    const int e = index(expiry), last = index(times.back());
    Vector coupon(last + 1, 0.0);
    std::vector<char> exercise(last + 1, 0);
    double previous = expiry;
    for (double t : times) {
        coupon[index(t)] = strike * (t - previous);
        if (bermudan && index(t) < last) exercise[index(t)] = 1;
        previous = t;
    }
    exercise[e] = 1;
    const int n = 2 * jmax_ + 1;
    Vector fixed(n, 1.0 + coupon[last]), option(n, 0.0);  // value of the remaining fixed leg including notional, and of the option
    for (int i = last - 1; i >= 0; --i) {
        if (i >= e) rollback(fixed, i);
        rollback(option, i);
        if (i >= e && exercise[i])
            for (int k = 0; k < n; ++k) option[k] = std::max(option[k], payer ? 1.0 - fixed[k] : fixed[k] - 1.0);
        if (i >= e)
            for (int k = 0; k < n; ++k) fixed[k] += coupon[i];
    }
    return option[jmax_];
}

}  // namespace ck
