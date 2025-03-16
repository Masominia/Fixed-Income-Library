// Command line front end: curvekit <command> <quotes.txt> key=value ...
// Output: `key value` lines and `row ...` lines, read by the Python package.
//
//   curve      bootstrap and print pillars and forward rates
//   bond       price and risk measures of a fixed-coupon bond
//   nelson     Nelson-Siegel fit to the bootstrapped zero rates
//   hw         Hull-White closed forms against the trinomial tree
//   vols       Black volatilities implied by Hull-White swaption prices
//   calibrate  recover Hull-White parameters from swaption volatilities
//   shortrate  Vasicek and Cox-Ingersoll-Ross bond prices against Monte Carlo

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "curvekit/black.hpp"
#include "curvekit/bond.hpp"
#include "curvekit/bootstrap.hpp"
#include "curvekit/calibration.hpp"
#include "curvekit/hullwhite.hpp"
#include "curvekit/hwtree.hpp"
#include "curvekit/montecarlo.hpp"
#include "curvekit/nelson_siegel.hpp"
#include "curvekit/quotes.hpp"
#include "curvekit/swap.hpp"

using namespace ck;

struct Args {
    std::map<std::string, std::string> kv;
    double num(const std::string& key, double fallback) const {
        const auto it = kv.find(key);
        return it == kv.end() ? fallback : std::stod(it->second);
    }
    int count(const std::string& key, int fallback) const { return static_cast<int>(num(key, fallback)); }
    std::string text(const std::string& key, const std::string& fallback) const {
        const auto it = kv.find(key);
        return it == kv.end() ? fallback : it->second;
    }
    Vector list(const std::string& key, const std::string& fallback) const {
        Vector out;
        std::stringstream ss(text(key, fallback));
        std::string cell;
        while (std::getline(ss, cell, ',')) out.push_back(std::stod(cell));
        return out;
    }
};

static int run(const std::string& command, const std::string& path, const Args& a) {
    if (command == "shortrate") {
        const int paths = a.count("paths", 100000), steps_per_year = a.count("steps_per_year", 50);
        const Vasicek v{a.num("speed", 0.3), a.num("mean", 0.04), a.num("sigma", 0.012), a.num("r0", 0.02)};
        const CoxIngersollRoss c{a.num("speed", 0.3), a.num("mean", 0.04), a.num("cir_sigma", 0.08), a.num("r0", 0.02)};
        std::printf("feller %d\n", c.feller() ? 1 : 0);
        for (double t : a.list("maturities", "1,2,5,10")) {
            const int steps = static_cast<int>(std::lround(t * steps_per_year));
            const Estimate ev = mc_bond(v, t, steps, paths, 11), ec = mc_bond(c, t, steps, paths, 12);
            std::printf("row vasicek %.6g %.12g %.12g %.6g\n", t, v.bond(t), ev.value, ev.std_error);
            std::printf("row cir %.6g %.12g %.12g %.6g\n", t, c.bond(t), ec.value, ec.std_error);
        }
        return 0;
    }

    const std::vector<Quote> quotes = read_quotes(path);
    const Bootstrap boot = bootstrap(quotes, parse_interpolation(a.text("mode", "monotone_cubic").c_str()));
    const Curve& curve = boot.curve;
    const double a_hw = a.num("a", 0.05), sigma = a.num("sigma", 0.01);
    const int freq = a.count("frequency", 2);

    if (command == "curve") {
        std::printf("mode %s\npasses %d\nmax_error %.6g\n", interpolation_name(curve.mode()), boot.passes, boot.max_error);
        for (std::size_t i = 0; i < curve.times().size(); ++i)
            std::printf("row pillar %.10g %.12g %.12g\n", curve.times()[i], curve.zero_rates()[i], curve.discount(curve.times()[i]));
        const double step = a.num("step", 0.05), last = curve.times().back();
        for (double t = 0.0; t <= last + 1e-9; t += step) std::printf("row forward %.10g %.12g %.12g\n", t, curve.instantaneous_forward(t), curve.zero(t));
        return 0;
    }

    if (command == "bond") {
        const Bond bond{a.num("maturity", 10.0), a.num("coupon", 0.04), freq, 100.0};
        const double p = price(bond, curve), y = yield_to_maturity(bond, p);
        const YieldAnalytics ya = yield_analytics(bond, y);
        std::printf("price %.12g\nyield %.12g\nmacaulay %.12g\nmodified %.12g\nconvexity %.12g\ndv01 %.12g\n", p, y, ya.macaulay, ya.modified, ya.convexity, ya.dv01);
        std::printf("effective_duration %.12g\n", effective_duration(bond, curve));
        const Vector krd = key_rate_durations(bond, curve);
        for (std::size_t i = 0; i < krd.size(); ++i) std::printf("row krd %.10g %.12g\n", curve.times()[i], krd[i]);
        for (double bp : a.list("shifts", "-200,-100,-50,-10,10,50,100,200")) {
            const double dy = bp * 1e-4, exact = price_from_yield(bond, y + dy);
            const double first = p * (1.0 - ya.modified * dy), second = first + 0.5 * p * ya.convexity * dy * dy;
            std::printf("row shift %.6g %.12g %.12g %.12g\n", bp, exact, first, second);
        }
        return 0;
    }

    if (command == "nelson") {
        const NelsonSiegelFit fit = fit_nelson_siegel(curve.times(), curve.zero_rates());
        std::printf("level %.10g\nslope %.10g\ncurvature %.10g\ntau %.10g\nrmse %.6g\nmax_error %.6g\n", fit.model.level, fit.model.slope, fit.model.curvature, fit.model.tau, fit.rmse,
                    fit.max_error);
        for (std::size_t i = 0; i < curve.times().size(); ++i) std::printf("row fit %.10g %.12g %.12g\n", curve.times()[i], curve.zero_rates()[i], fit.model.zero(curve.times()[i]));
        return 0;
    }

    if (command == "hw") {
        const HullWhite model(curve, a_hw, sigma);
        const double expiry = a.num("expiry", 2.0), tenor = a.num("tenor", 5.0), horizon = expiry + tenor;
        const double atm = par_rate(curve, expiry, horizon, freq), strike = atm + a.num("offset", 0.0);
        const double bond_strike = curve.discount(horizon) / curve.discount(expiry);
        const double cf_bond = model.bond_option(expiry, horizon, bond_strike, true);
        const double cf_caplet = model.caplet(expiry, expiry + 0.5, atm);
        const double cf_payer = model.swaption(expiry, tenor, freq, strike, true), cf_receiver = model.swaption(expiry, tenor, freq, strike, false);
        const Swap forward_swap{expiry, horizon, strike, freq, true, 1.0};
        std::printf("atm %.12g\nstrike %.12g\n", atm, strike);
        std::printf("bond_option %.12g\ncaplet %.12g\npayer %.12g\nreceiver %.12g\n", cf_bond, cf_caplet, cf_payer, cf_receiver);
        std::printf("parity_gap %.6g\n", cf_payer - cf_receiver - swap_value(curve, forward_swap));
        std::printf("black_caplet_vol %.12g\n", implied_black_vol(cf_caplet / (0.5 * curve.discount(expiry + 0.5)), curve.forward(expiry, expiry + 0.5), atm, expiry, true));
        for (double n : a.list("steps", "4,8,16,32,64")) {
            const HullWhiteTree tree(curve, a_hw, sigma, horizon, static_cast<int>(n));
            double fit = 0.0;
            for (double t = 0.5; t <= horizon + 1e-9; t += 0.5) fit = std::max(fit, std::fabs(tree.discount_bond(t) / curve.discount(t) - 1.0));
            const double gross = 1.0 + atm * 0.5;
            std::printf("row tree %d %d %.6g %.12g %.12g %.12g %.12g %.12g\n", static_cast<int>(n), tree.width(), fit, tree.bond_option(expiry, horizon, bond_strike, true),
                        gross * tree.bond_option(expiry, expiry + 0.5, 1.0 / gross, false), tree.swaption(expiry, tenor, freq, strike, true, false),
                        tree.swaption(expiry, tenor, freq, strike, false, false), tree.swaption(expiry, tenor, freq, strike, true, true));
        }
        return 0;
    }

    if (command == "vols") {
        for (double e : a.list("expiries", "1,2,5,10"))
            for (double t : a.list("tenors", "1,2,5,10")) std::printf("row vol %.6g %.6g %.12g %.12g\n", e, t, model_black_vol(curve, a_hw, sigma, e, t, freq), par_rate(curve, e, e + t, freq));
        return 0;
    }

    if (command == "calibrate") {
        std::vector<SwaptionQuote> market;
        const double noise = a.num("noise", 0.0);
        int k = 0;
        for (double e : a.list("expiries", "1,2,5,10"))
            for (double t : a.list("tenors", "1,2,5,10")) {
                // deterministic pseudo-noise so that the run is reproducible
                const double wobble = noise * std::sin(12.9898 * ++k);
                market.push_back({e, t, model_black_vol(curve, a_hw, sigma, e, t, freq) * (1.0 + wobble)});
            }
        const Calibration c = calibrate_hull_white(curve, market, freq);
        std::printf("true_a %.10g\ntrue_sigma %.10g\nfitted_a %.10g\nfitted_sigma %.10g\nrmse_vol %.6g\nevaluations %d\nquotes %zu\n", a_hw, sigma, c.mean_reversion, c.sigma, c.rmse_vol,
                    c.evaluations, market.size());
        return 0;
    }

    throw std::invalid_argument("unknown command " + command);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: curvekit <curve|bond|nelson|hw|vols|calibrate|shortrate> <quotes.txt> key=value ...\n");
        return 1;
    }
    Args a;
    for (int i = 3; i < argc; ++i) {
        const std::string s = argv[i];
        const auto eq = s.find('=');
        if (eq != std::string::npos) a.kv[s.substr(0, eq)] = s.substr(eq + 1);
    }
    try {
        return run(argv[1], argv[2], a);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "curvekit: %s\n", e.what());
        return 1;
    }
}
